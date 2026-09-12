#include "geometry.h"

#include <cmath>
#include <tuple>

namespace mxi {

std::pair<double, double> Panel::parallax_offset(double mm_fast,
                                                 double mm_slow) const {
  if (!parallax || mu <= 0.0 || thickness <= 0.0) return {0.0, 0.0};
  const Vec3 u = lab_coord_mm(mm_fast, mm_slow).normalized();
  const double cos_theta = std::abs(u.dot(normal()));
  if (cos_theta <= 0.0) return {0.0, 0.0};
  const double attenuation_length = 1.0 / mu;
  // Mean depth of interaction along the ray, truncated by the sensor: a photon
  // that gets through the full thickness is not recorded at all, which is the
  // second term.
  const double path = thickness / cos_theta;
  const double transmitted = std::exp(-mu * path);
  double depth = attenuation_length - (path + attenuation_length) * transmitted;
  if (parallax_conditional) {
    const double absorbed = 1.0 - transmitted;
    if (absorbed > 0.0) depth /= absorbed;
  }
  return {depth * u.dot(fast), depth * u.dot(slow)};
}

std::pair<double, double> Panel::px_to_mm(double px_fast,
                                          double px_slow) const {
  const double raw_fast = px_fast * pixel_size[0];
  const double raw_slow = px_slow * pixel_size[1];
  const auto offset = parallax_offset(raw_fast, raw_slow);
  return {raw_fast - offset.first, raw_slow - offset.second};
}

std::pair<double, double> Panel::mm_to_px(double mm_fast,
                                          double mm_slow) const {
  const auto offset = parallax_offset(mm_fast, mm_slow);
  return {(mm_fast + offset.first) / pixel_size[0],
          (mm_slow + offset.second) / pixel_size[1]};
}

Scan Scan::from_oscillation(const std::vector<double> &oscillation_deg,
                           std::int64_t first, std::int64_t last) {
  Scan s;
  s.first_image = first;
  s.last_image = last;
  if (oscillation_deg.empty()) return s;
  s.osc_start = oscillation_deg.front();
  const std::size_t n = oscillation_deg.size();
  if (n < 2) return s;

  s.osc_width = (oscillation_deg.back() - oscillation_deg.front()) /
                static_cast<double>(n - 1);
  if (s.osc_width == 0.0) return s;
  for (std::size_t i = 1; i < n; ++i) {
    const double width = oscillation_deg[i] - oscillation_deg[i - 1];
    s.max_width_deviation =
        std::fmax(s.max_width_deviation,
                  std::abs(width - s.osc_width) / std::abs(s.osc_width));
  }
  return s;
}

Goniometer Goniometer::from_axes(const std::vector<Vec3> &axes,
                                 const std::vector<double> &angles_deg,
                                 std::size_t scan_axis) {
  Goniometer g;
  if (scan_axis >= axes.size()) return g;
  g.axis = axes[scan_axis].normalized();

  // The axes are numbered from the SAMPLE outwards towards the laboratory.
  // In the l-cysteine goniometer the sample is attached to GON_PHI (axes[0]),
  // which sits on GON_OMEGA (axes[1]), which is bolted to the floor.
  //
  // So a vector fixed to the sample is carried first by axes[0], then by
  // axes[1], and so on: each axis further out applies LATER and therefore
  // multiplies on the LEFT. Accumulating the other way round composes the
  // stack inside out.
  //
  // With only one axis on either side of the scan axis the two orders are the
  // same matrix, which is why no dataset to hand can tell them apart -- the
  // l-cysteine goniometer has two axes, so at most one is ever below the scan
  // axis. This order is from the physical arrangement, not from a measurement.
  Mat3 fixed = Mat3::identity();
  for (std::size_t i = 0; i < scan_axis; ++i) {
    fixed = rotation(axes[i], Scan::radians(angles_deg[i])) * fixed;
  }
  Mat3 setting = Mat3::identity();
  for (std::size_t i = scan_axis + 1; i < axes.size(); ++i) {
    setting = rotation(axes[i], Scan::radians(angles_deg[i])) * setting;
  }
  g.fixed = fixed;
  g.setting = setting;
  return g;
}

std::optional<std::pair<double, double>> Panel::intersect(const Vec3 &s1) const {
  // Solve origin + u*fast_mm + v*slow_mm = t*s1_hat for (u, v, t).
  const Vec3 d = s1.normalized();
  const Mat3 basis = Mat3::from_columns(fast, slow, -d);
  bool ok = false;
  const Mat3 inv = basis.inverse(&ok);
  if (!ok) return std::nullopt;  // ray parallel to the panel

  const Vec3 uvt = inv * (-origin);
  const double t = uvt.z;
  if (!(t > 0.0)) return std::nullopt;  // behind the sample

  // uvt is where the ray meets the panel face, in millimetres; the pixel that
  // fires is displaced from it by the parallax offset.
  const auto px = mm_to_px(uvt.x, uvt.y);
  const double px_fast = px.first;
  const double px_slow = px.second;

  // A ray aimed at the exact corner of the panel comes back at -4e-13 pixels
  // rather than zero, because the intersection is solved in millimetres and
  // divided back. A bare `< 0` test then rejects it. That matters in two
  // places: a reflection predicted exactly on an edge is silently lost, and on
  // a tiled detector a ray striking the seam between two panels can belong to
  // neither of them.
  //
  // A nanopixel is nine orders of magnitude below anything physically
  // meaningful and three above the round-off, so it separates the two cases
  // cleanly without widening the panel in any sense that matters.
  constexpr double kEdgeTolerance = 1e-9;
  if (px_fast < -kEdgeTolerance ||
      px_fast >= static_cast<double>(image_size[0]) + kEdgeTolerance) {
    return std::nullopt;
  }
  if (px_slow < -kEdgeTolerance ||
      px_slow >= static_cast<double>(image_size[1]) + kEdgeTolerance) {
    return std::nullopt;
  }
  return std::make_pair(px_fast, px_slow);
}

std::optional<std::tuple<std::size_t, double, double>> Detector::intersect(
    const Vec3 &s1) const {
  for (std::size_t i = 0; i < panels.size(); ++i) {
    if (auto hit = panels[i].intersect(s1)) {
      return std::make_tuple(i, hit->first, hit->second);
    }
  }
  return std::nullopt;
}

double UnitCell::volume() const {
  const double ca = std::cos(Scan::radians(alpha));
  const double cb = std::cos(Scan::radians(beta));
  const double cg = std::cos(Scan::radians(gamma));
  const double t = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg;
  return t > 0.0 ? a * b * c * std::sqrt(t) : 0.0;
}

Crystal Crystal::from_real_space(const Vec3 &a, const Vec3 &b, const Vec3 &c) {
  Crystal crystal;
  // Real-space vectors are the ROWS of A inverse, so A is the plain inverse of
  // the matrix built from them -- not the inverse transpose. Both are diagonal
  // for an orthogonal cell, so the wrong one survives any cubic test case.
  crystal.A = Mat3::from_rows(a, b, c).inverse();
  return crystal;
}

Vec3 Crystal::real_a() const { return A.inverse().row(0); }
Vec3 Crystal::real_b() const { return A.inverse().row(1); }
Vec3 Crystal::real_c() const { return A.inverse().row(2); }

UnitCell Crystal::cell() const {
  const Mat3 real = A.inverse();
  const Vec3 a = real.row(0), b = real.row(1), c = real.row(2);
  UnitCell u;
  u.a = a.norm();
  u.b = b.norm();
  u.c = c.norm();
  auto angle = [](const Vec3 &p, const Vec3 &q) {
    const double d = p.dot(q) / (p.norm() * q.norm());
    return Scan::degrees(std::acos(std::fmax(-1.0, std::fmin(1.0, d))));
  };
  u.alpha = angle(b, c);
  u.beta = angle(a, c);
  u.gamma = angle(a, b);
  return u;
}

Mat3 Crystal::A_at(double t) const {
  if (A_points.size() < 2) return A;
  const double span = static_cast<double>(A_points.size() - 1);
  // Clamped rather than extrapolated. A reflection predicted a little outside
  // the scan should be modelled by the nearest end of it, not by a linear
  // continuation of whatever the crystal was doing when the scan stopped.
  const double u = std::fmax(0.0, std::fmin(1.0, t)) * span;
  const auto i = static_cast<std::size_t>(std::fmin(std::floor(u), span - 1.0));
  const double f = u - static_cast<double>(i);
  const Mat3 &a = A_points[i];
  const Mat3 &b = A_points[i + 1];
  Mat3 out;
  for (std::size_t k = 0; k < 9; ++k) out.m[k] = a.m[k] * (1.0 - f) + b.m[k] * f;
  return out;
}

Mat3 Experiment::setting_at(double z) const {
  if (!crystal) return Mat3::identity();
  if (!crystal->scan_varying()) return crystal->A;
  const double n = static_cast<double>(scan.num_images());
  return crystal->A_at(n > 0.0 ? z / n : 0.0);
}

double Crystal::d_spacing(int h, int k, int l) const {
  const Vec3 q = A * Vec3{static_cast<double>(h), static_cast<double>(k),
                          static_cast<double>(l)};
  const double n = q.norm();
  return n > 0.0 ? 1.0 / n : 0.0;
}

Vec3 lab_scattering_vector(const Experiment &e, std::size_t panel,
                           double px_fast, double px_slow) {
  const Vec3 lab = e.detector[panel].lab_coord_px(px_fast, px_slow);
  const Vec3 s1 = lab.normalized() / e.beam.wavelength;
  return s1 - e.beam.s0();
}

Vec3 reciprocal_lattice_point(const Experiment &e, std::size_t panel,
                              double px_fast, double px_slow, double z) {
  const Vec3 q = lab_scattering_vector(e, panel, px_fast, px_slow);
  const double phi = e.scan.phi_from_z(z);
  // Undo the goniometer rotation. Inverting the rotation by negating the angle
  // rather than inverting the matrix: it is exact for the rotation part and
  // avoids a determinant on a matrix that is orthogonal by construction.
  const Mat3 inverse_rotation = e.goniometer.fixed.transpose() *
                                rotation(e.goniometer.axis, -phi) *
                                e.goniometer.setting.transpose();
  return inverse_rotation * q;
}

std::vector<Vec3> reciprocal_lattice_points(const Experiment &e,
                                            const std::vector<Observation> &obs) {
  std::vector<Vec3> out;
  out.reserve(obs.size());
  for (const Observation &o : obs) {
    out.push_back(
        reciprocal_lattice_point(e, o.panel, o.px_fast, o.px_slow, o.z));
  }
  return out;
}

}  // namespace mxi
