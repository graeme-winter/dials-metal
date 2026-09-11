#include "predict.h"

#include <algorithm>
#include <cmath>

namespace mxi {

namespace {

constexpr double kPi = 3.14159265358979323846;

// Put an angle into the half-open interval starting at `from`, so that a
// predicted phi can be compared against the scan range without either end
// needing a special case.
double wrap_from(double phi, double from) {
  const double two_pi = 2.0 * kPi;
  double d = std::fmod(phi - from, two_pi);
  if (d < 0.0) d += two_pi;
  return from + d;
}

}  // namespace

Intersections ewald_intersections(const Experiment &e, const Vec3 &r0) {
  Intersections out;

  const Vec3 m2 = e.goniometer.axis.normalized();
  // Work in the frame the rotation is about: move the incident beam through
  // the inverse setting rotation rather than rotating the axis, so the
  // decomposition below is about a fixed vector.
  const Vec3 s0p = e.goniometer.setting.transpose() * e.beam.s0();
  const Vec3 u = e.goniometer.fixed * r0;

  const double r_squared = u.norm_squared();
  if (r_squared <= 0.0) return out;
  // Beyond twice the Ewald radius nothing can reach the sphere at any angle.
  if (r_squared > 4.0 * e.beam.s0().norm_squared()) return out;

  const double u_parallel = u.dot(m2);
  const double s_parallel = s0p.dot(m2);
  const double a0 = u_parallel * s_parallel;
  const double b = s0p.dot(u) - a0;
  const double c = s0p.dot(m2.cross(u));
  const double d = -0.5 * r_squared - a0;

  const double amplitude = std::sqrt(b * b + c * c);
  if (amplitude <= 0.0) return out;
  const double cosine = d / amplitude;
  // A point whose component along the axis never brings it to the sphere.
  // The equality case is a tangential grazing, which is a measure-zero event
  // that acos handles without special-casing.
  if (cosine < -1.0 || cosine > 1.0) return out;

  const double centre = std::atan2(c, b);
  const double offset = std::acos(cosine);

  out.any = true;
  out.phi[0] = centre - offset;
  out.phi[1] = centre + offset;

  // Entering or exiting the sphere, from the sign of d/dphi of the distance to
  // the sphere centre. The derivative of s0 . R r at the crossing is
  // s0' . (m2 x R r), and the point is entering while that is negative.
  //
  // CONVENTION: this is a belief about DIALS' definition and is unvalidated;
  // if it is backwards, every `entering` flag is inverted and every keyed join
  // against DIALS output pairs the wrong two observations. See
  // docs/conventions.md.
  for (int i = 0; i < 2; ++i) {
    const Vec3 rotated = rotation(m2, out.phi[i]) * u;
    out.entering[i] = s0p.dot(m2.cross(rotated)) < 0.0;
  }
  return out;
}

std::array<int, 3> index_bounds(const Crystal &crystal, double d_min,
                                int max_index) {
  const Mat3 real = crystal.A.inverse();
  std::array<int, 3> bounds{};
  for (int i = 0; i < 3; ++i) {
    const double length = real.row(static_cast<std::size_t>(i)).norm();
    const double limit = d_min > 0.0 ? length / d_min : 0.0;
    bounds[static_cast<std::size_t>(i)] =
        std::min(max_index, static_cast<int>(std::ceil(limit)) + 1);
  }
  return bounds;
}

namespace {

// Turn one Ewald intersection into a Prediction, or reject it.
bool build(const Experiment &e, const PredictOptions &options, int h, int k,
           int l, const Vec3 &r0, double phi, bool entering,
           Prediction *out) {
  if (!options.allow_outside_scan) {
    const double start = e.scan.phi_start();
    const double end = e.scan.phi_end();
    const double lo = std::min(start, end);
    const double span = std::abs(end - start);
    // Wrapped into the scan's own interval, so a scan crossing 360 degrees and
    // a reflection predicted at -179 still meet.
    if (span < 2.0 * kPi) {
      const double wrapped = wrap_from(phi, lo);
      if (wrapped > lo + span) return false;
      phi = wrapped;
    }
  }

  const Mat3 r = e.goniometer.rotation_at(phi);
  const Vec3 s1 = e.beam.s0() + r * r0;

  auto hit = e.detector.intersect(s1);
  if (!hit) return false;

  out->h = h;
  out->k = k;
  out->l = l;
  out->entering = entering;
  out->phi = phi;
  out->panel = std::get<0>(*hit);
  out->px_fast = std::get<1>(*hit);
  out->px_slow = std::get<2>(*hit);
  out->z = e.scan.z_from_phi(phi);
  out->s1 = s1;
  return true;
}

}  // namespace

std::vector<Prediction> predict_indices(
    const Experiment &e, const std::vector<std::array<int, 3>> &indices,
    const PredictOptions &options) {
  std::vector<Prediction> out;
  if (!e.crystal) return out;
  const Mat3 &A = e.crystal->A;

  for (const std::array<int, 3> &hkl : indices) {
    if (hkl[0] == 0 && hkl[1] == 0 && hkl[2] == 0) continue;
    const Vec3 r0 = A * Vec3{static_cast<double>(hkl[0]),
                             static_cast<double>(hkl[1]),
                             static_cast<double>(hkl[2])};
    const Intersections cross = ewald_intersections(e, r0);
    if (!cross.any) continue;
    for (int i = 0; i < 2; ++i) {
      Prediction p;
      if (build(e, options, hkl[0], hkl[1], hkl[2], r0, cross.phi[i],
                cross.entering[i], &p)) {
        out.push_back(p);
      }
    }
  }
  return out;
}

std::vector<Prediction> predict(const Experiment &e,
                                const PredictOptions &options) {
  std::vector<Prediction> out;
  if (!e.crystal) return out;

  double d_min = options.d_min;
  if (d_min <= 0.0) {
    // The absolute limit: a reciprocal lattice point further out than the
    // diameter of the Ewald sphere can never diffract.
    d_min = 0.5 * e.beam.wavelength;
  }
  const double q_max = 1.0 / d_min;
  const std::array<int, 3> bounds =
      index_bounds(*e.crystal, d_min, options.max_index);
  const Mat3 &A = e.crystal->A;

  // One independent unit of work per (h, k, l). Nothing in the body touches
  // anything outside it, which is what makes this the natural first kernel;
  // the only shared state is the output, and on a device that becomes an
  // atomic append or a compaction pass.
  for (int h = -bounds[0]; h <= bounds[0]; ++h) {
    for (int k = -bounds[1]; k <= bounds[1]; ++k) {
      for (int l = -bounds[2]; l <= bounds[2]; ++l) {
        if (h == 0 && k == 0 && l == 0) continue;
        const Vec3 r0 = A * Vec3{static_cast<double>(h), static_cast<double>(k),
                                 static_cast<double>(l)};
        if (r0.norm() > q_max) continue;

        const Intersections cross = ewald_intersections(e, r0);
        if (!cross.any) continue;
        for (int i = 0; i < 2; ++i) {
          Prediction p;
          if (build(e, options, h, k, l, r0, cross.phi[i], cross.entering[i],
                    &p)) {
            out.push_back(p);
          }
        }
      }
    }
  }
  return out;
}

}  // namespace mxi
