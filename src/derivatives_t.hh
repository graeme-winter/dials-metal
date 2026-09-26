// The analytical derivatives, written once and compiled at either precision.
//
// Same discipline as `target.h`, and for the same reason: a separate float
// implementation disagreeing with the double one could be the precision or a
// transcription error, and no measurement can tell those apart. Compiled with
// `double` this must reproduce `crystal_derivatives`, `detector_derivatives`
// and `beam_derivatives` exactly; only then does the `float` answer mean
// anything.
//
// The question it settles: a finite difference of the target does not survive
// float32 -- measured at 100 per cent relative error at the step refinement
// uses -- because the residual is a difference of positions of order two
// thousand pixels and its absolute error is epsilon times the position. An
// analytical derivative never forms that difference. Whether it therefore
// keeps its digits is an expectation until measured, and the last expectation
// of this kind was wrong by the whole result.

#pragma once

#include <cmath>

#include "target.hh"

namespace mxi {

//: Everything about one predicted reflection the derivatives need, at
//: precision T. The double-precision twin of this is `PredictionState`.
template <typename T>
struct TargetState {
  bool valid = false;
  T r0[3], r_phi[3], s1[3], v[3];
  T D[9];      // inverse of (fast | slow | origin) as columns
  T R[9];      // goniometer rotation at the diffracting angle
  T phi = 0;
  T volume = 0;  // (axis x r_phi) . s0, the denominator of eqn (40)
};

namespace target_detail {

template <typename T>
inline bool invert(const T (&m)[9], T (&out)[9]) {
  const T det = m[0] * (m[4] * m[8] - m[5] * m[7]) -
                m[1] * (m[3] * m[8] - m[5] * m[6]) +
                m[2] * (m[3] * m[7] - m[4] * m[6]);
  if (!(std::abs(det) > T(0))) return false;
  const T i = T(1) / det;
  out[0] = (m[4] * m[8] - m[5] * m[7]) * i;
  out[1] = (m[2] * m[7] - m[1] * m[8]) * i;
  out[2] = (m[1] * m[5] - m[2] * m[4]) * i;
  out[3] = (m[5] * m[6] - m[3] * m[8]) * i;
  out[4] = (m[0] * m[8] - m[2] * m[6]) * i;
  out[5] = (m[2] * m[3] - m[0] * m[5]) * i;
  out[6] = (m[3] * m[7] - m[4] * m[6]) * i;
  out[7] = (m[1] * m[6] - m[0] * m[7]) * i;
  out[8] = (m[0] * m[4] - m[1] * m[3]) * i;
  return true;
}

}  // namespace target_detail

template <typename T>
TargetState<T> target_state(const TargetModel<T> &m, int h, int k, int l, T z) {
  using target_detail::cross;
  using target_detail::dot;
  using target_detail::invert;
  using target_detail::mat_mul;
  using target_detail::mat_vec;
  using target_detail::rodrigues;
  TargetState<T> s;

  const T hkl[3] = {static_cast<T>(h), static_cast<T>(k), static_cast<T>(l)};
  mat_vec(m.A, hkl, s.r0);

  // The Ewald solution, exactly as in evaluate_target.
  T u[3];
  mat_vec(m.fixed, s.r0, u);
  const T r_squared = dot(u, u);
  if (!(r_squared > T(0))) return s;

  const T setting_t[9] = {m.setting[0], m.setting[3], m.setting[6],
                          m.setting[1], m.setting[4], m.setting[7],
                          m.setting[2], m.setting[5], m.setting[8]};
  T s0p[3];
  mat_vec(setting_t, m.s0, s0p);
  const T a0 = dot(u, m.axis) * dot(s0p, m.axis);
  const T b = dot(s0p, u) - a0;
  T axis_cross_u[3];
  cross(m.axis, u, axis_cross_u);
  const T c = dot(s0p, axis_cross_u);
  const T d = -T(0.5) * r_squared - a0;
  const T amplitude = std::sqrt(b * b + c * c);
  if (!(amplitude > T(0))) return s;
  const T cosine = d / amplitude;
  if (cosine < T(-1) || cosine > T(1)) return s;

  const T two_pi = static_cast<T>(6.283185307179586476925286766559);
  const T degrees = static_cast<T>(0.017453292519943295769236907685);
  const T phi_obs = (m.osc_start + (z - m.z_offset) * m.osc_width) * degrees;
  const T centre = std::atan2(c, b);
  const T offset = std::acos(cosine);

  s.phi = phi_obs;
  T best = static_cast<T>(1e30);
  for (int i = 0; i < 2; ++i) {
    const T root = i == 0 ? centre - offset : centre + offset;
    T gap = std::fmod(root - phi_obs, two_pi);
    if (gap > T(0.5) * two_pi) gap -= two_pi;
    if (gap < T(-0.5) * two_pi) gap += two_pi;
    if (std::abs(gap) < best) {
      best = std::abs(gap);
      s.phi = phi_obs + gap;
    }
  }

  T rot[9], tmp[9];
  rodrigues(m.axis, s.phi, rot);
  mat_mul(m.setting, rot, tmp);
  mat_mul(tmp, m.fixed, s.R);
  mat_vec(s.R, s.r0, s.r_phi);
  for (int i = 0; i < 3; ++i) s.s1[i] = m.s0[i] + s.r_phi[i];

  // d with the panel basis as COLUMNS, so d (X, Y, 1)^T = alpha s1.
  const T dmat[9] = {m.fast[0], m.slow[0], m.origin[0],
                     m.fast[1], m.slow[1], m.origin[1],
                     m.fast[2], m.slow[2], m.origin[2]};
  if (!invert(dmat, s.D)) return s;
  mat_vec(s.D, s.s1, s.v);

  T axis_cross_r[3];
  cross(m.axis, s.r_phi, axis_cross_r);
  s.volume = dot(axis_cross_r, m.s0);
  s.valid = std::abs(s.v[2]) > T(0);
  return s;
}

//: dX, dY in millimetres on the panel, dphi in radians.
template <typename T>
struct Derivative3 {
  T dX = 0, dY = 0, dphi = 0;
};

namespace target_detail {

// The shared tail of every derivative: given d(r_phi)/dp, produce dX, dY.
template <typename T>
inline void finish(const TargetState<T> &s, const T (&dr_phi)[3],
                   Derivative3<T> *out) {
  T dv[3];
  mat_vec(s.D, dr_phi, dv);
  const T w = s.v[2];
  out->dX = (w * dv[0] - s.v[0] * dv[2]) / (w * w);
  out->dY = (w * dv[1] - s.v[1] * dv[2]) / (w * w);
}

}  // namespace target_detail

template <typename T>
void crystal_derivatives_t(const TargetModel<T> &m, const TargetState<T> &s,
                           int h, int k, int l, Derivative3<T> (&out)[9]) {
  using target_detail::cross;
  using target_detail::dot;
  using target_detail::finish;
  using target_detail::mat_vec;
  if (!s.valid || !(std::abs(s.volume) > T(0))) return;
  const T hkl[3] = {static_cast<T>(h), static_cast<T>(k), static_cast<T>(l)};
  T e_cross_r[3];
  cross(m.axis, s.r_phi, e_cross_r);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      T dr0[3] = {T(0), T(0), T(0)};
      dr0[i] = hkl[j];
      T R_dr0[3];
      mat_vec(s.R, dr0, R_dr0);
      const T dphi = -dot(R_dr0, s.s1) / s.volume;
      T dr_phi[3];
      for (int c = 0; c < 3; ++c) dr_phi[c] = e_cross_r[c] * dphi + R_dr0[c];
      out[i * 3 + j].dphi = dphi;
      finish(s, dr_phi, &out[i * 3 + j]);
    }
  }
}

template <typename T>
void detector_derivatives_t(const TargetModel<T> &m, const TargetState<T> &s,
                            Derivative3<T> (&out)[6]) {
  using target_detail::cross;
  using target_detail::mat_vec;
  if (!s.valid) return;
  const T w = s.v[2];

  for (int k = 0; k < 6; ++k) {
    T dfast[3] = {T(0), T(0), T(0)};
    T dslow[3] = {T(0), T(0), T(0)};
    T dorigin[3] = {T(0), T(0), T(0)};
    if (k < 3) {
      dorigin[k] = T(1);
    } else {
      T axis[3] = {T(0), T(0), T(0)};
      axis[k - 3] = T(1);
      cross(axis, m.fast, dfast);
      cross(axis, m.slow, dslow);
      const T offset[3] = {m.origin[0] - m.centre[0], m.origin[1] - m.centre[1],
                           m.origin[2] - m.centre[2]};
      cross(axis, offset, dorigin);
    }
    // -D (dd/dp) v, with dd/dp having those three as its columns.
    const T dd[9] = {dfast[0], dslow[0], dorigin[0], dfast[1], dslow[1],
                     dorigin[1], dfast[2], dslow[2], dorigin[2]};
    T ddv[3], dv[3];
    mat_vec(dd, s.v, ddv);
    mat_vec(s.D, ddv, dv);
    out[k].dphi = T(0);
    out[k].dX = -(w * dv[0] - s.v[0] * dv[2]) / (w * w);
    out[k].dY = -(w * dv[1] - s.v[1] * dv[2]) / (w * w);
  }
}

template <typename T>
void beam_derivatives_t(const TargetModel<T> &m, const TargetState<T> &s,
                        T wavelength, Derivative3<T> (&out)[2]) {
  using target_detail::cross;
  using target_detail::dot;
  using target_detail::finish;
  if (!s.valid || !(std::abs(s.volume) > T(0))) return;

  // s0 points along propagation, so the direction it was built from is its
  // negation; the perpendicular basis has to be built from the same vector the
  // double version uses or the two parameters mean different things.
  const T length = std::sqrt(dot(m.s0, m.s0));
  const T d[3] = {-m.s0[0] / length, -m.s0[1] / length, -m.s0[2] / length};
  const T z_axis[3] = {T(0), T(0), T(1)};
  T u[3];
  cross(z_axis, d, u);
  if (!(std::sqrt(dot(u, u)) > T(1e-6))) {
    const T x_axis[3] = {T(1), T(0), T(0)};
    cross(x_axis, d, u);
  }
  const T u_length = std::sqrt(dot(u, u));
  for (int i = 0; i < 3; ++i) u[i] /= u_length;
  T v[3];
  cross(d, u, v);

  T e_cross_r[3];
  cross(m.axis, s.r_phi, e_cross_r);
  const T *basis[2] = {u, v};
  for (int k = 0; k < 2; ++k) {
    T ds0[3];
    for (int i = 0; i < 3; ++i) ds0[i] = -basis[k][i] / wavelength;
    const T dphi = -dot(s.r_phi, ds0) / s.volume;
    T dr_phi[3];
    for (int i = 0; i < 3; ++i) dr_phi[i] = e_cross_r[i] * dphi + ds0[i];
    out[k].dphi = dphi;
    finish(s, dr_phi, &out[k]);
  }
}

}  // namespace mxi
