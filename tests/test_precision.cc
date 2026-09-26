// Can the target be computed in single precision?
//
// The answer turns out to be: the residual yes, a finite-difference derivative
// no, and the second half of that is the opposite of what this project
// estimated before measuring it.
//
// The estimate went: a parameter step of 1e-6 relative changes the residual by
// 1.5e-3 of its own size, float32 has an epsilon of 1.2e-7, so about four
// digits survive. That reasoning is wrong, and the error is worth keeping
// written down. The residual is a DIFFERENCE of detector positions of order
// two thousand pixels, so its absolute error in float32 is not epsilon times
// the residual but epsilon times the position -- about 2.4e-4 px. The change
// being measured is 1.5e-3 x 0.32 px, about 5e-4 px. Signal and noise are the
// same size and nothing survives at all.
//
// Measuring beats estimating, and the way to measure is to compile one
// implementation at both precisions rather than write two.

#include <algorithm>
#include <cmath>
#include <vector>

#include "../src/derivatives.hh"
#include "../src/predict.hh"
#include "../src/refine.hh"
#include "../src/target.hh"
#include "check.hh"
#include "real_data.hh"
#include "../src/derivatives_t.hh"

namespace mxi {

namespace {

Vec3 vec(const double (&a)[3]) { return {a[0], a[1], a[2]}; }

Experiment insulin_experiment() {
  Experiment e;
  e.beam.direction = vec(real::kBeamDirection);
  e.beam.wavelength = real::kWavelength;
  Panel p;
  p.fast = vec(real::kRefinedFast);
  p.slow = vec(real::kRefinedSlow);
  p.origin = vec(real::kRefinedOrigin);
  p.pixel_size[0] = p.pixel_size[1] = real::kRefinedPixelSize[0];
  p.image_size[0] = real::kImageSize[0];
  p.image_size[1] = real::kImageSize[1];
  p.parallax = true;
  p.mu = real::kRefinedMu;
  p.thickness = real::kRefinedThickness;
  e.detector.panels.push_back(p);
  std::vector<Vec3> axes;
  std::vector<double> angles;
  for (int i = 0; i < 3; ++i) {
    axes.push_back(vec(real::kGonioAxes[i]));
    angles.push_back(real::kGonioAngles[i]);
  }
  e.goniometer = Goniometer::from_axes(axes, angles,
                                       static_cast<std::size_t>(real::kScanAxis));
  e.scan.first_image = real::kImageRange[0];
  e.scan.last_image = real::kImageRange[1];
  e.scan.osc_start = real::kOscStart;
  e.scan.osc_width = real::kOscWidth;
  e.crystal = Crystal::from_real_space(vec(real::kRealSpaceA), vec(real::kRealSpaceB),
                                       vec(real::kRealSpaceC));
  return e;
}

double percentile(std::vector<double> v, double p) {
  if (v.empty()) return -1.0;
  std::sort(v.begin(), v.end());
  return v[static_cast<std::size_t>(p * static_cast<double>(v.size() - 1))];
}

}  // namespace

TEST(the_templated_target_reproduces_the_real_one_in_double) {
  // Until this holds, nothing measured with the float version means anything:
  // a disagreement could be the precision or could be a transcription error,
  // and the measurement cannot tell them apart.
  const Experiment e = insulin_experiment();
  std::vector<double> difference;
  for (const real::Row &r : real::rows()) {
    const Residual reference =
        centroid_residual(e, 0, r.h, r.k, r.l, r.px_fast, r.px_slow, r.px_z);
    if (!reference.valid) continue;
    const auto model = narrow<double>(e, 0, e.setting_at(r.px_z));
    const auto got = evaluate_target<double>(model, r.h, r.k, r.l, r.px_fast,
                                             r.px_slow, r.px_z);
    check::is_true(got.valid, "the templated version must predict too");
    difference.push_back(std::hypot(got.dx - reference.dx, got.dy - reference.dy));
  }
  check::is_true(difference.size() > 30, "enough reflections");
  // Measured at 5e-13 px: the same arithmetic in a different order.
  check::is_true(percentile(difference, 0.99) < 1e-10,
                 "the transcription must be exact to rounding");
}

TEST(the_residual_itself_survives_single_precision) {
  // 2.2e-4 px median against residuals of 0.32 px -- under a tenth of a per
  // cent, random per reflection, averaging away over thirteen thousand of
  // them. This half of the question is fine.
  const Experiment e = insulin_experiment();
  std::vector<double> difference;
  for (const real::Row &r : real::rows()) {
    const auto wide = narrow<double>(e, 0, e.setting_at(r.px_z));
    const auto narrowed = narrow<float>(e, 0, e.setting_at(r.px_z));
    const auto a = evaluate_target<double>(wide, r.h, r.k, r.l, r.px_fast,
                                           r.px_slow, r.px_z);
    const auto b = evaluate_target<float>(
        narrowed, r.h, r.k, r.l, static_cast<float>(r.px_fast),
        static_cast<float>(r.px_slow), static_cast<float>(r.px_z));
    if (!a.valid || !b.valid) continue;
    difference.push_back(std::hypot(static_cast<double>(b.dx) - a.dx,
                                    static_cast<double>(b.dy) - a.dy));
  }
  check::is_true(difference.size() > 30, "enough reflections");
  check::is_true(percentile(difference, 0.5) < 1e-3, "median float error");
  check::is_true(percentile(difference, 0.99) < 5e-3, "and the tail");
}

TEST(a_finite_difference_derivative_does_not_survive_single_precision) {
  // The part that was estimated wrongly. This asserts the FAILURE, so that
  // nobody later assumes numerical differentiation would port as it stands --
  // the test would then fail and say so.
  //
  // Measured against the analytical derivative: in double at the step the
  // refinement uses, the median relative error is 4e-7. In float at the same
  // step it is 1.0, which is to say the derivative is entirely noise. At a
  // step near sqrt(epsilon), where a float finite difference is at its best,
  // it is 1.4e-2 -- under two digits, with a ninety-ninth percentile above
  // five, meaning some entries have the wrong sign.
  const Experiment e = insulin_experiment();
  double scale = 0.0;
  for (double m : e.crystal->A.m) scale = std::fmax(scale, std::abs(m));

  std::vector<double> in_double, in_float, in_float_big;
  for (const real::Row &r : real::rows()) {
    const PredictionState s = prediction_state(e, 0, r.h, r.k, r.l, r.px_z);
    if (!s.valid || std::abs(s.volume) < 0.05) continue;
    const auto exact = crystal_derivatives(s, r.h, r.k, r.l);
    double J[4];
    e.detector[0].mm_to_px_jacobian(s.v.x / s.v.z, s.v.y / s.v.z, J);
    const Mat3 A = e.setting_at(r.px_z);

    for (std::size_t p = 0; p < 9; ++p) {
      // d(residual_x)/dp, which is minus the derivative of the prediction.
      const double reference = -(J[0] * exact[p].dX + J[1] * exact[p].dY);
      if (std::abs(reference) < 1.0) continue;

      const auto difference = [&](auto tag, double step) -> double {
        using T = decltype(tag);
        Mat3 moved = A;
        moved.m[p] += step;
        const auto m0 = narrow<T>(e, 0, A);
        const auto m1 = narrow<T>(e, 0, moved);
        const auto a = evaluate_target<T>(m0, r.h, r.k, r.l,
                                          static_cast<T>(r.px_fast),
                                          static_cast<T>(r.px_slow),
                                          static_cast<T>(r.px_z));
        const auto b = evaluate_target<T>(m1, r.h, r.k, r.l,
                                          static_cast<T>(r.px_fast),
                                          static_cast<T>(r.px_slow),
                                          static_cast<T>(r.px_z));
        if (!a.valid || !b.valid) return -1.0;
        const double fd = (static_cast<double>(b.dx) - static_cast<double>(a.dx)) / step;
        return std::abs(fd - reference) / std::abs(reference);
      };

      const double d0 = difference(double(), 1e-6 * scale);
      const double f0 = difference(float(), 1e-6 * scale);
      const double f1 = difference(float(), 3e-4 * scale);
      if (d0 >= 0) in_double.push_back(d0);
      if (f0 >= 0) in_float.push_back(f0);
      if (f1 >= 0) in_float_big.push_back(f1);
    }
  }
  check::is_true(in_double.size() > 100, "enough comparisons");

  // Double is fine, which is what makes the float numbers meaningful rather
  // than a sign that the analytical reference is wrong.
  check::is_true(percentile(in_double, 0.5) < 1e-5, "double is accurate");
  // Float at the current step is worthless.
  check::is_true(percentile(in_float, 0.5) > 0.1,
                 "float at this step must be shown to fail, not assumed to work");
  // And even at its best step it does not reach three digits.
  check::is_true(percentile(in_float_big, 0.5) > 1e-3,
                 "float at sqrt(eps) is still not good enough");
}

// --------------------------------------------------------------------------
// the analytical derivative in single precision
// --------------------------------------------------------------------------


namespace {

// Every derivative for one reflection, flattened, at precision T.
template <typename T>
std::vector<double> templated_derivatives(const Experiment &e, int h, int k,
                                          int l, double z) {
  std::vector<double> out;
  const auto model = narrow<T>(e, 0, e.setting_at(z));
  const auto state = target_state<T>(model, h, k, l, static_cast<T>(z));
  if (!state.valid) return out;
  Derivative3<T> crystal[9], detector[6], beam[2];
  crystal_derivatives_t<T>(model, state, h, k, l, crystal);
  detector_derivatives_t<T>(model, state, detector);
  beam_derivatives_t<T>(model, state, static_cast<T>(e.beam.wavelength), beam);
  const auto push = [&out](const Derivative3<T> &d) {
    out.push_back(static_cast<double>(d.dX));
    out.push_back(static_cast<double>(d.dY));
    out.push_back(static_cast<double>(d.dphi));
  };
  for (const auto &d : crystal) push(d);
  for (const auto &d : detector) push(d);
  for (const auto &d : beam) push(d);
  return out;
}

std::vector<double> reference_derivatives(const Experiment &e,
                                          const PredictionState &s, int h,
                                          int k, int l) {
  std::vector<double> out;
  const auto crystal = crystal_derivatives(s, h, k, l);
  const auto detector = detector_derivatives(s, e.detector[0]);
  const auto beam = beam_derivatives(s, e.beam);
  const auto push = [&out](const CentroidDerivative &d) {
    out.push_back(d.dX);
    out.push_back(d.dY);
    out.push_back(d.dphi);
  };
  for (const auto &d : crystal) push(d);
  for (const auto &d : detector) push(d);
  for (const auto &d : beam) push(d);
  return out;
}

}  // namespace

TEST(the_templated_derivatives_reproduce_the_real_ones_in_double) {
  // Bit for bit, measured: the same operations in the same order. Anything
  // less and the float measurement below would be confounded by a
  // transcription difference.
  const Experiment e = insulin_experiment();
  double worst = 0.0;
  std::size_t compared = 0;
  for (const real::Row &r : real::rows()) {
    const PredictionState s = prediction_state(e, 0, r.h, r.k, r.l, r.px_z);
    if (!s.valid || std::abs(s.volume) < 0.05) continue;
    const std::vector<double> reference = reference_derivatives(e, s, r.h, r.k, r.l);
    const std::vector<double> got =
        templated_derivatives<double>(e, r.h, r.k, r.l, r.px_z);
    if (got.size() != reference.size()) continue;
    for (std::size_t i = 0; i < got.size(); ++i) {
      if (std::abs(reference[i]) < 1e-8) continue;
      worst = std::fmax(worst, std::abs(got[i] - reference[i]) / std::abs(reference[i]));
      ++compared;
    }
  }
  check::is_true(compared > 1000, "enough components compared");
  check::close(worst, 0.0, 0.0, "identical, not merely close");
}

TEST(the_analytical_derivative_does_survive_single_precision) {
  // The measurement this whole exercise was for, and the one that decides
  // whether a device port is possible at all.
  //
  //   finite difference, float, step 1e-6   median relative error 1.00
  //   analytical,        float              median relative error 1.3e-07
  //
  // Seven orders of magnitude, and the reason is structural rather than lucky:
  // an analytical derivative never forms the difference of two nearly equal
  // positions, so there is no cancellation to spend the significance on.
  const Experiment e = insulin_experiment();
  std::vector<double> relative;
  for (const real::Row &r : real::rows()) {
    const PredictionState s = prediction_state(e, 0, r.h, r.k, r.l, r.px_z);
    if (!s.valid || std::abs(s.volume) < 0.05) continue;
    const std::vector<double> wide =
        templated_derivatives<double>(e, r.h, r.k, r.l, r.px_z);
    const std::vector<double> narrowed =
        templated_derivatives<float>(e, r.h, r.k, r.l, r.px_z);
    if (wide.size() != narrowed.size() || wide.empty()) continue;
    for (std::size_t i = 0; i < wide.size(); ++i) {
      if (std::abs(wide[i]) < 1e-8) continue;
      relative.push_back(std::abs(narrowed[i] - wide[i]) / std::abs(wide[i]));
    }
  }
  check::is_true(relative.size() > 1000, "enough components compared");
  // Measured: median 1.3e-7, ninety-ninth percentile 1.2e-5.
  check::is_true(percentile(relative, 0.5) < 1e-6, "median, about seven digits");
  check::is_true(percentile(relative, 0.99) < 1e-3, "and the tail holds up");
}

}  // namespace mxi
