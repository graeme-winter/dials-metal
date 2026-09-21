// Summation integration and the quantum efficiency.

#include <cmath>
#include <vector>

#include "../src/integrate.h"
#include "check.h"

using namespace mxi;

namespace {

Panel silicon_panel() {
  Panel p;
  p.fast = {1.0, 0.0, 0.0};
  p.slow = {0.0, -1.0, 0.0};
  p.pixel_size[0] = p.pixel_size[1] = 0.075;
  p.image_size[0] = 2068;
  p.image_size[1] = 2162;
  p.origin = {-77.5, 81.1, -168.5};
  p.mu = 3.6631;
  p.thickness = 0.450;
  return p;
}

//: A box with a flat background and a block of signal in the middle.
Shoebox planted(double background, double signal_per_pixel, int signal_pixels) {
  Shoebox box;
  box.panel = 0;
  box.bbox[0] = 0;
  box.bbox[1] = 15;
  box.bbox[2] = 0;
  box.bbox[3] = 15;
  box.bbox[4] = 0;
  box.bbox[5] = 3;
  const std::size_t n = box.size();
  box.data.assign(n, static_cast<float>(background));
  box.background.assign(n, 0.0f);
  box.mask.assign(n, static_cast<std::uint8_t>(shoebox_mask::kValid |
                                               shoebox_mask::kBackground));
  int placed = 0;
  for (std::int32_t z = 1; z < 2 && placed < signal_pixels; ++z) {
    for (std::int32_t y = 5; y < 10 && placed < signal_pixels; ++y) {
      for (std::int32_t x = 5; x < 10 && placed < signal_pixels; ++x) {
        const std::size_t at = box.at(x, y, z);
        box.mask[at] = shoebox_mask::kValid | shoebox_mask::kForeground;
        box.data[at] = static_cast<float>(background + signal_per_pixel);
        ++placed;
      }
    }
  }
  return box;
}

}  // namespace

TEST(the_quantum_efficiency_is_the_fraction_the_sensor_stops) {
  // 1 - exp(-mu t / cos theta). Checked against the qe column of a DIALS
  // integrated.refl, where it is identical for every reflection to 2e-16; this
  // pins the two limits that fix the formula.
  const Panel p = silicon_panel();
  // Straight through the sensor: the shortest path, so the least absorbed.
  const Vec3 normal = p.fast.cross(p.slow).normalized();
  const double straight = quantum_efficiency(p, normal * 1.05);
  check::close(straight, 1.0 - std::exp(-p.mu * p.thickness), 1e-12,
               "normal incidence is 1 - exp(-mu t)");

  // Obliquely: a longer path through the same silicon, so more is stopped.
  const Vec3 oblique = (normal + p.fast * 0.5).normalized();
  const double slanted = quantum_efficiency(p, oblique);
  check::is_true(slanted > straight, "an oblique ray is absorbed more");

  const double cosine = std::fabs(oblique.dot(normal));
  check::close(slanted, 1.0 - std::exp(-p.mu * p.thickness / cosine), 1e-12,
               "and by exactly the extra path length");

  // A panel with no absorption model stops everything, rather than nothing:
  // returning zero would divide the intensity by zero downstream.
  Panel bare = p;
  bare.mu = 0.0;
  check::close(quantum_efficiency(bare, normal), 1.0, 0.0,
               "no absorption model means no correction");
}

TEST(a_planted_intensity_comes_back) {
  // Flat background, a known number of signal pixels each with a known excess.
  // The answer is arithmetic, not statistics -- so the robust estimator is
  // turned off with a large tuning constant, leaving the plain mean.
  //
  // With it on, a background of exactly 2.0 in every pixel comes back as
  // 2.068, and that is the estimator being right rather than wrong: a constant
  // is not a Poisson sample, and the Fisher consistency correction assumes the
  // data are drawn from the distribution it is correcting for. The GLM has its
  // own tests, on data that are.
  const double background = 2.0;
  const double excess = 37.0;
  const int pixels = 25;
  Shoebox box = planted(background, excess, pixels);
  IntegrateOptions options;
  options.background.tuning = 1e6;
  const IntegratedReflection r = integrate_shoebox(&box, options);
  check::is_true(r.valid, "integrated");
  check::equal(static_cast<long long>(r.n_foreground), pixels, "foreground count");
  check::close(r.background_mean, background, 1e-6, "background recovered");
  check::close(r.intensity, excess * pixels, 1e-4 * excess * pixels,
               "and the intensity is the excess over it");
}

TEST(the_variance_is_leslies_and_the_background_term_is_in_it) {
  // var = G [I + I_bg + (m/n) I_bg]. The third term is the uncertainty in the
  // background estimate and is what makes a large background region worth
  // having; a variance that omitted it would be too small for every weak
  // reflection in a dataset.
  const double background = 3.0;
  Shoebox box = planted(background, 10.0, 25);
  IntegrateOptions options;
  options.background.tuning = 1e6;  // the mean, so the arithmetic is exact
  const IntegratedReflection r = integrate_shoebox(&box, options);
  check::is_true(r.valid, "integrated");

  const double m = static_cast<double>(r.n_foreground);
  const double n = static_cast<double>(r.n_background);
  const double foreground_sum = r.intensity + r.background_sum;
  const double expected = foreground_sum + (m / n) * r.background_sum;
  check::close(r.variance, expected, 1e-6 * expected, "Leslie equation 11");
  check::is_true(r.background_sum_variance > 0.0,
                 "and the background term is not zero");
  check::is_true(r.variance > foreground_sum,
                 "so the variance exceeds the raw Poisson noise");
}

TEST(a_negative_intensity_still_has_a_positive_variance) {
  // A weak reflection whose foreground happens to fall below the background.
  // The intensity is negative and that is not an error -- throwing it away or
  // clamping it biases every merged intensity upwards. The variance must stay
  // positive, which is why it is written in terms of the counts rather than
  // the difference.
  Shoebox box = planted(5.0, 0.0, 25);
  for (std::size_t i = 0; i < box.size(); ++i) {
    if (box.mask[i] & shoebox_mask::kForeground) box.data[i] = 3.0f;
  }
  const IntegratedReflection r = integrate_shoebox(&box);
  check::is_true(r.valid, "integrated");
  check::is_true(r.intensity < 0.0, "the intensity is negative");
  check::is_true(r.variance > 0.0, "and the variance is not");
}

TEST(the_gain_multiplies_the_variance_and_not_the_intensity) {
  Shoebox box = planted(2.0, 20.0, 25);
  IntegrateOptions options;
  const IntegratedReflection one = integrate_shoebox(&box, options);
  Shoebox again = planted(2.0, 20.0, 25);
  options.gain = 3.0;
  const IntegratedReflection three = integrate_shoebox(&again, options);
  check::close(three.intensity, one.intensity, 1e-9, "same intensity");
  check::close(three.variance, 3.0 * one.variance, 1e-9 * three.variance,
               "three times the variance");
}

TEST(too_few_background_pixels_is_refused_rather_than_guessed) {
  Shoebox box = planted(2.0, 20.0, 25);
  for (std::size_t i = 0; i < box.size(); ++i) {
    if (box.mask[i] & shoebox_mask::kBackground) box.mask[i] = 0;
  }
  const IntegratedReflection r = integrate_shoebox(&box);
  check::is_true(!r.valid, "not integrated");
  check::is_true(r.too_few_background, "and it says why");
}

TEST(the_fitted_background_is_left_in_the_shoebox) {
  // So that a saved shoebox carries what was subtracted from it, and a picture
  // of one can be believed.
  Shoebox box = planted(4.0, 15.0, 25);
  IntegrateOptions options;
  options.background.tuning = 1e6;
  const IntegratedReflection r = integrate_shoebox(&box, options);
  check::is_true(r.valid, "integrated");
  for (float v : box.background) {
    check::close(v, r.background_mean, 1e-5, "background written back");
  }
}
