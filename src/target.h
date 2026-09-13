// The refinement target, written once and compiled at whatever precision.
//
// The question this exists to answer: can the target be computed in float32,
// as it would have to be on an Apple GPU, which has no double precision at
// all?
//
// The trap in answering it is that a separate float implementation and the
// existing double one could disagree for two quite different reasons -- a
// transcription error, or the precision itself -- and the measurement cannot
// tell them apart. So this is a single implementation, templated on the scalar
// type. Compiled with `double` it must reproduce `centroid_residual` exactly;
// only once that holds does the `float` result mean anything.
//
// What is NOT in here, deliberately: the scan-varying spline. Its evaluation
// is a convex combination of four matrices with weights summing to one, which
// is about the most benign arithmetic in the pipeline, and templating it would
// add code without adding a question. The caller passes the setting matrix it
// wants used. A device implementation would evaluate the spline in float too,
// and that is worth measuring separately rather than folding in here.

#pragma once

#include <cmath>

#include "geometry.h"

namespace mxi {

//: The model a target evaluation needs, at the precision it will be evaluated
//: in. Built from an Experiment by `narrow`.
template <typename T>
struct TargetModel {
  T A[9];         // setting matrix, row-major
  T s0[3];        // incident beam, already negated
  T axis[3];      // rotation axis, laboratory frame
  T fixed[9];     // goniometer fixed rotation
  T setting[9];   // goniometer setting rotation
  T fast[3], slow[3], origin[3];
  T pixel_size[2];
  T mu = 0, thickness = 0;
  bool parallax = false;
  T osc_start = 0, osc_width = 0, z_offset = 0;
};

//: Observed minus calculated, in pixels and images, at precision T.
template <typename T>
struct TargetResidual {
  bool valid = false;
  T dx = 0, dy = 0, dz = 0;
};

// Build the model at precision T from the double-precision experiment. The
// setting matrix is passed separately so the caller decides whether the spline
// was evaluated in double or not.
template <typename T>
TargetModel<T> narrow(const Experiment &e, std::size_t panel, const Mat3 &A) {
  TargetModel<T> m;
  for (std::size_t i = 0; i < 9; ++i) {
    m.A[i] = static_cast<T>(A.m[i]);
    m.fixed[i] = static_cast<T>(e.goniometer.fixed.m[i]);
    m.setting[i] = static_cast<T>(e.goniometer.setting.m[i]);
  }
  const Vec3 s0 = e.beam.s0();
  const Vec3 axis = e.goniometer.axis.normalized();
  const Panel &p = e.detector[panel];
  for (std::size_t i = 0; i < 3; ++i) {
    m.s0[i] = static_cast<T>(s0[i]);
    m.axis[i] = static_cast<T>(axis[i]);
    m.fast[i] = static_cast<T>(p.fast[i]);
    m.slow[i] = static_cast<T>(p.slow[i]);
    m.origin[i] = static_cast<T>(p.origin[i]);
  }
  m.pixel_size[0] = static_cast<T>(p.pixel_size[0]);
  m.pixel_size[1] = static_cast<T>(p.pixel_size[1]);
  m.mu = static_cast<T>(p.mu);
  m.thickness = static_cast<T>(p.thickness);
  m.parallax = p.parallax;
  m.osc_start = static_cast<T>(e.scan.osc_start);
  m.osc_width = static_cast<T>(e.scan.osc_width);
  m.z_offset = static_cast<T>(e.scan.z_offset);
  return m;
}

namespace target_detail {

template <typename T>
inline void mat_vec(const T (&m)[9], const T (&v)[3], T (&out)[3]) {
  out[0] = m[0] * v[0] + m[1] * v[1] + m[2] * v[2];
  out[1] = m[3] * v[0] + m[4] * v[1] + m[5] * v[2];
  out[2] = m[6] * v[0] + m[7] * v[1] + m[8] * v[2];
}

template <typename T>
inline T dot(const T (&a)[3], const T (&b)[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T>
inline void cross(const T (&a)[3], const T (&b)[3], T (&out)[3]) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

template <typename T>
inline void rodrigues(const T (&u)[3], T angle, T (&out)[9]) {
  const T c = std::cos(angle), s = std::sin(angle), t = T(1) - c;
  out[0] = t * u[0] * u[0] + c;
  out[1] = t * u[0] * u[1] - s * u[2];
  out[2] = t * u[0] * u[2] + s * u[1];
  out[3] = t * u[0] * u[1] + s * u[2];
  out[4] = t * u[1] * u[1] + c;
  out[5] = t * u[1] * u[2] - s * u[0];
  out[6] = t * u[0] * u[2] - s * u[1];
  out[7] = t * u[1] * u[2] + s * u[0];
  out[8] = t * u[2] * u[2] + c;
}

template <typename T>
inline void mat_mul(const T (&a)[9], const T (&b)[9], T (&out)[9]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      T sum = T(0);
      for (int k = 0; k < 3; ++k) sum += a[i * 3 + k] * b[k * 3 + j];
      out[i * 3 + j] = sum;
    }
  }
}

}  // namespace target_detail

// The whole chain: reciprocal lattice point, Ewald solution, rotation, panel
// intersection, parallax, residual. Mirrors `centroid_residual` step for step.
template <typename T>
TargetResidual<T> evaluate_target(const TargetModel<T> &m, int h, int k, int l,
                                  T px_fast, T px_slow, T z) {
  using namespace target_detail;
  TargetResidual<T> out;

  const T hkl[3] = {static_cast<T>(h), static_cast<T>(k), static_cast<T>(l)};
  T r0[3];
  mat_vec(m.A, hkl, r0);

  // Ewald: solve A0 + B cos(phi) + C sin(phi) = -|r|^2 / 2 in the frame the
  // rotation acts in.
  T u[3];
  mat_vec(m.fixed, r0, u);
  const T r_squared = dot(u, u);
  if (!(r_squared > T(0))) return out;

  T s0p[3];
  const T setting_t[9] = {m.setting[0], m.setting[3], m.setting[6],
                          m.setting[1], m.setting[4], m.setting[7],
                          m.setting[2], m.setting[5], m.setting[8]};
  mat_vec(setting_t, m.s0, s0p);

  const T u_parallel = dot(u, m.axis);
  const T s_parallel = dot(s0p, m.axis);
  const T a0 = u_parallel * s_parallel;
  const T b = dot(s0p, u) - a0;
  T axis_cross_u[3];
  cross(m.axis, u, axis_cross_u);
  const T c = dot(s0p, axis_cross_u);
  const T d = -T(0.5) * r_squared - a0;

  const T amplitude = std::sqrt(b * b + c * c);
  if (!(amplitude > T(0))) return out;
  const T cosine = d / amplitude;
  if (cosine < T(-1) || cosine > T(1)) return out;
  const T centre = std::atan2(c, b);
  const T offset = std::acos(cosine);

  const T two_pi = static_cast<T>(6.283185307179586476925286766559);
  const T degrees = static_cast<T>(0.017453292519943295769236907685);
  const T phi_obs = (m.osc_start + (z - m.z_offset) * m.osc_width) * degrees;

  // The root nearer the observation.
  T phi = phi_obs;
  T best = static_cast<T>(1e30);
  for (int i = 0; i < 2; ++i) {
    const T root = i == 0 ? centre - offset : centre + offset;
    T gap = std::fmod(root - phi_obs, two_pi);
    if (gap > T(0.5) * two_pi) gap -= two_pi;
    if (gap < T(-0.5) * two_pi) gap += two_pi;
    if (std::abs(gap) < best) {
      best = std::abs(gap);
      phi = phi_obs + gap;
    }
  }

  // s1 = s0 + R(phi) r0, with R = setting * rot(axis, phi) * fixed.
  T rot[9], tmp[9], R[9];
  rodrigues(m.axis, phi, rot);
  mat_mul(m.setting, rot, tmp);
  mat_mul(tmp, m.fixed, R);
  T r_phi[3];
  mat_vec(R, r0, r_phi);
  const T s1[3] = {m.s0[0] + r_phi[0], m.s0[1] + r_phi[1], m.s0[2] + r_phi[2]};

  // Panel intersection: solve origin + uu * fast + vv * slow = tt * s1_hat.
  const T s1_length = std::sqrt(dot(s1, s1));
  if (!(s1_length > T(0))) return out;
  const T dir[3] = {s1[0] / s1_length, s1[1] / s1_length, s1[2] / s1_length};

  // Columns fast, slow, -dir; solved by Cramer's rule, which is what a device
  // would do for a three by three rather than calling a factorisation.
  T col2[3] = {-dir[0], -dir[1], -dir[2]};
  T cs[3];
  cross(m.slow, col2, cs);
  const T det = dot(m.fast, cs);
  if (std::abs(det) <= T(0)) return out;
  const T rhs[3] = {-m.origin[0], -m.origin[1], -m.origin[2]};
  // Cramer: with columns (fast, slow, col2), the solution components are the
  // determinants with each column replaced by the right-hand side. Note the
  // ORDER of every cross product below -- reversing one flips the sign of that
  // component, and the first version of this had two of them backwards, which
  // showed up as every ray appearing to travel away from the detector.
  T c1[3], c2[3];
  cross(rhs, col2, c1);     // det(fast, rhs, col2) = fast . (rhs x col2)
  cross(m.slow, rhs, c2);   // det(fast, slow, rhs) = fast . (slow x rhs)
  const T mm_fast = dot(rhs, cs) / det;
  const T mm_slow = dot(m.fast, c1) / det;
  const T tt = dot(m.fast, c2) / det;
  if (!(tt > T(0))) return out;

  // Parallax, then millimetres to pixels.
  T off_fast = T(0), off_slow = T(0);
  if (m.parallax && m.mu > T(0) && m.thickness > T(0)) {
    const T lab[3] = {m.origin[0] + mm_fast * m.fast[0] + mm_slow * m.slow[0],
                      m.origin[1] + mm_fast * m.fast[1] + mm_slow * m.slow[1],
                      m.origin[2] + mm_fast * m.fast[2] + mm_slow * m.slow[2]};
    const T length = std::sqrt(dot(lab, lab));
    const T unit[3] = {lab[0] / length, lab[1] / length, lab[2] / length};
    T normal[3];
    cross(m.fast, m.slow, normal);
    const T n_length = std::sqrt(dot(normal, normal));
    const T n[3] = {normal[0] / n_length, normal[1] / n_length,
                    normal[2] / n_length};
    const T cos_theta = std::abs(dot(unit, n));
    if (cos_theta > T(0)) {
      const T attenuation = T(1) / m.mu;
      const T path = m.thickness / cos_theta;
      const T depth = attenuation - (path + attenuation) * std::exp(-m.mu * path);
      off_fast = depth * dot(unit, m.fast);
      off_slow = depth * dot(unit, m.slow);
    }
  }

  const T cal_fast = (mm_fast + off_fast) / m.pixel_size[0];
  const T cal_slow = (mm_slow + off_slow) / m.pixel_size[1];
  const T cal_z = m.osc_width != T(0)
                      ? (phi / degrees - m.osc_start) / m.osc_width + m.z_offset
                      : m.z_offset;

  out.valid = true;
  out.dx = px_fast - cal_fast;
  out.dy = px_slow - cal_slow;
  out.dz = z - cal_z;
  return out;
}

}  // namespace mxi
