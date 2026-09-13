#include "derivatives.h"

#include <algorithm>
#include <cmath>

#include "predict.h"

namespace mxi {

PredictionState prediction_state(const Experiment &e, std::size_t panel, int h,
                                 int k, int l, double z) {
  PredictionState s;
  if (!e.crystal || panel >= e.detector.size()) return s;

  const Vec3 hkl{static_cast<double>(h), static_cast<double>(k),
                 static_cast<double>(l)};
  s.r0 = e.setting_at(z) * hkl;

  const Intersections cross = ewald_intersections(e, s.r0);
  if (!cross.any) return s;

  // The root nearer the observation, matching what the refinement target does.
  const double phi_obs = e.scan.phi_from_z(z);
  const double two_pi = 2.0 * 3.14159265358979323846;
  double best_gap = 1e30;
  for (int i = 0; i < 2; ++i) {
    double gap = std::fmod(cross.phi[i] - phi_obs, two_pi);
    if (gap > 0.5 * two_pi) gap -= two_pi;
    if (gap < -0.5 * two_pi) gap += two_pi;
    if (std::abs(gap) < best_gap) {
      best_gap = std::abs(gap);
      s.phi = phi_obs + gap;
    }
  }

  s.R_phi = e.goniometer.rotation_at(s.phi);
  s.r_phi = s.R_phi * s.r0;
  s.s1 = e.beam.s0() + s.r_phi;
  s.axis = e.goniometer.lab_axis();

  const Panel &p = e.detector[panel];
  bool ok = false;
  // d has the panel basis and origin as its COLUMNS, so that
  // d (X, Y, 1)^T = alpha s1 reproduces eqn (1). Building it from rows would
  // transpose the detector, which is the obvious way to get this wrong.
  s.D = Mat3::from_columns(p.fast, p.slow, p.origin).inverse(&ok);
  if (!ok) return s;
  s.v = s.D * s.s1;

  s.volume = s.axis.cross(s.r_phi).dot(e.beam.s0());
  s.valid = std::abs(s.v.z) > 0.0;
  return s;
}

std::array<CentroidDerivative, 9> crystal_derivatives(const PredictionState &s,
                                                      int h, int k, int l) {
  std::array<CentroidDerivative, 9> out{};
  if (!s.valid || s.volume == 0.0) return out;

  const double hkl[3] = {static_cast<double>(h), static_cast<double>(k),
                         static_cast<double>(l)};
  const Vec3 e_cross_r = s.axis.cross(s.r_phi);
  const double w = s.v.z;
  const double w2 = w * w;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      // r0 = A h, so the derivative with respect to element A(i, j) is h_j
      // sitting in row i and nothing elsewhere. This is where refining A
      // directly, rather than through U and B, pays for itself.
      Vec3 dr0{0.0, 0.0, 0.0};
      dr0[i] = hkl[j];

      // eqn (40). ds0/dp is zero: the beam does not depend on a crystal
      // parameter.
      const Vec3 R_dr0 = s.R_phi * dr0;
      const double dphi = -R_dr0.dot(s.s1) / s.volume;

      // eqn (46). The first term is the constraint that the only motion of
      // r_phi the diffraction condition permits is the one that keeps it on
      // the Ewald sphere.
      const Vec3 dr_phi = e_cross_r * dphi + R_dr0;

      // eqn (45), with the detector fixed, then eqn (43).
      const Vec3 dv = s.D * dr_phi;
      CentroidDerivative &d = out[i * 3 + j];
      d.dphi = dphi;
      d.dX = (w * dv.x - s.v.x * dv.z) / w2;
      d.dY = (w * dv.y - s.v.y * dv.z) / w2;
    }
  }
  return out;
}

namespace {

// The panel centre, in the laboratory frame. Rotations act about this point.
Vec3 panel_centre(const Panel &p) {
  return p.lab_coord_mm(0.5 * static_cast<double>(p.image_size[0]) * p.pixel_size[0],
                        0.5 * static_cast<double>(p.image_size[1]) * p.pixel_size[1]);
}

}  // namespace

Panel perturb_panel(const Panel &p, const double shift[6]) {
  Panel out = p;
  const Vec3 turn{shift[3], shift[4], shift[5]};
  const double angle = turn.norm();
  if (angle > 0.0) {
    const Vec3 centre = panel_centre(p);
    const Mat3 r = rotation(turn / angle, angle);
    out.fast = r * p.fast;
    out.slow = r * p.slow;
    out.origin = centre + r * (p.origin - centre);
  }
  out.origin += Vec3{shift[0], shift[1], shift[2]};
  return out;
}

Beam perturb_beam(const Beam &b, const double shift[2]) {
  Beam out = b;
  const Vec3 d = b.direction.normalized();
  Vec3 u = Vec3{0.0, 0.0, 1.0}.cross(d);
  if (u.norm() < 1e-6) u = Vec3{1.0, 0.0, 0.0}.cross(d);
  u = u.normalized();
  const Vec3 v = d.cross(u);
  out.direction = (d + u * shift[0] + v * shift[1]).normalized();
  return out;
}

std::array<CentroidDerivative, 6> detector_derivatives(const PredictionState &s,
                                                       const Panel &p) {
  std::array<CentroidDerivative, 6> out{};
  if (!s.valid) return out;

  const Vec3 centre = panel_centre(p);
  const double w = s.v.z;
  const double w2 = w * w;

  for (std::size_t k = 0; k < 6; ++k) {
    // Columns of dd/dp: the derivative of (fast | slow | origin).
    Vec3 dfast{0.0, 0.0, 0.0}, dslow{0.0, 0.0, 0.0}, dorigin{0.0, 0.0, 0.0};
    if (k < 3) {
      dorigin[k] = 1.0;  // a pure translation moves only the origin
    } else {
      // The derivative of a rotation at the identity is the cross product with
      // the axis, so an infinitesimal rotation about laboratory axis j takes
      // any vector x to e_j x x.
      Vec3 axis{0.0, 0.0, 0.0};
      axis[k - 3] = 1.0;
      dfast = axis.cross(p.fast);
      dslow = axis.cross(p.slow);
      dorigin = axis.cross(p.origin - centre);
    }
    const Mat3 dd = Mat3::from_columns(dfast, dslow, dorigin);

    // eqn (47) with dr0/dp, ds0/dp and dphi/dp all zero: only the detector
    // term survives.
    const Vec3 dv = -(s.D * (dd * s.v));
    out[k].dphi = 0.0;
    out[k].dX = (w * dv.x - s.v.x * dv.z) / w2;
    out[k].dY = (w * dv.y - s.v.y * dv.z) / w2;
  }
  return out;
}

std::array<CentroidDerivative, 2> beam_derivatives(const PredictionState &s,
                                                   const Beam &b) {
  std::array<CentroidDerivative, 2> out{};
  if (!s.valid || s.volume == 0.0) return out;

  const Vec3 d = b.direction.normalized();
  Vec3 u = Vec3{0.0, 0.0, 1.0}.cross(d);
  if (u.norm() < 1e-6) u = Vec3{1.0, 0.0, 0.0}.cross(d);
  u = u.normalized();
  const Vec3 basis[2] = {u, d.cross(u)};

  const Vec3 e_cross_r = s.axis.cross(s.r_phi);
  const double w = s.v.z;
  const double w2 = w * w;

  for (std::size_t k = 0; k < 2; ++k) {
    // s0 = -direction / wavelength, and a tilt perpendicular to the direction
    // has no first-order effect on its length, so the normalisation drops out.
    const Vec3 ds0 = -basis[k] / b.wavelength;

    // eqn (40) with dr0/dp zero.
    const double dphi = -s.r_phi.dot(ds0) / s.volume;
    // eqn (47) with dd/dp and dr0/dp zero.
    const Vec3 ds1 = e_cross_r * dphi + ds0;
    const Vec3 dv = s.D * ds1;
    out[k].dphi = dphi;
    out[k].dX = (w * dv.x - s.v.x * dv.z) / w2;
    out[k].dY = (w * dv.y - s.v.y * dv.z) / w2;
  }
  return out;
}

SplineWeights spline_weights(const Experiment &e, double z) {
  SplineWeights out;
  if (!e.crystal || !e.crystal->scan_varying()) {
    out.count = 1;
    out.index[0] = 0;
    out.weight[0] = 1.0;
    return out;
  }
  const auto n = static_cast<long>(e.crystal->A_points.size());
  const double images = static_cast<double>(e.scan.num_images());
  const double t = images > 0.0 ? std::fmax(0.0, std::fmin(1.0, z / images)) : 0.0;
  const double segments = static_cast<double>(n + 1);
  const double u = t * segments;
  const auto i = static_cast<long>(std::fmin(std::floor(u), segments - 1.0));
  const double f = u - static_cast<double>(i);
  const double f2 = f * f, f3 = f2 * f;

  // The same basis as Crystal::A_at, and the same padding: segment i draws on
  // padded entries Q(i) .. Q(i+3), where Q(m) = P[clamp(m - 2, 0, n - 1)].
  const double b[4] = {(1.0 - 3.0 * f + 3.0 * f2 - f3) / 6.0,
                       (4.0 - 6.0 * f2 + 3.0 * f3) / 6.0,
                       (1.0 + 3.0 * f + 3.0 * f2 - 3.0 * f3) / 6.0, f3 / 6.0};
  out.count = 4;
  for (long m = 0; m < 4; ++m) {
    const long j = std::max(0L, std::min(i + m - 2, n - 1));
    out.index[m] = static_cast<std::size_t>(j);
    out.weight[m] = b[m];
  }
  return out;
}

}  // namespace mxi
