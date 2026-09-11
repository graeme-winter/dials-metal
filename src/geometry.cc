#include "geometry.h"

#include <cmath>
#include <tuple>

namespace mxi {

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

  const double px_fast = uvt.x / pixel_size[0];
  const double px_slow = uvt.y / pixel_size[1];

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

double Crystal::d_spacing(int h, int k, int l) const {
  const Vec3 q = A * Vec3{static_cast<double>(h), static_cast<double>(k),
                          static_cast<double>(l)};
  const double n = q.norm();
  return n > 0.0 ? 1.0 / n : 0.0;
}

Vec3 lab_scattering_vector(const Experiment &e, std::size_t panel,
                           double px_fast, double px_slow) {
  const Vec3 lab = e.detector[panel].lab_coord(px_fast, px_slow);
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
