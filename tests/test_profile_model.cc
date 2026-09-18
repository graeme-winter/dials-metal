// The profile model, against spots whose angular spread is known in advance.
//
// Matching dials.integrate's number is not by itself evidence of a correct
// estimator -- there are enough knobs in the recipe that a wrong one can be
// tuned onto a right answer. So these tests plant a spread and check it comes
// back, which is a statement about the estimator rather than about agreement.

#include <cmath>
#include <vector>

#include "../src/profile_model.h"
#include "../src/shoebox.h"
#include "check.h"

using namespace mxi;

namespace {

// A flat panel square to the beam, so the geometry is checkable by hand.
Experiment simple_experiment(double distance) {
  Experiment e;
  e.beam.direction = {0.0, 0.0, 1.0};
  e.beam.wavelength = 1.0;
  Panel p;
  p.fast = {1.0, 0.0, 0.0};
  p.slow = {0.0, 1.0, 0.0};
  p.pixel_size[0] = p.pixel_size[1] = 0.1;
  p.image_size[0] = p.image_size[1] = 1000;
  p.origin = {-50.0, -50.0, distance};
  e.detector.panels.push_back(p);
  e.goniometer.axis = {1.0, 0.0, 0.0};
  e.scan.first_image = 1;
  e.scan.last_image = 10;
  e.scan.osc_start = 0.0;
  e.scan.osc_width = 0.1;
  return e;
}

// A shoebox holding a single pixel of signal, at the box centre.
Shoebox one_pixel(std::int32_t x, std::int32_t y, float count) {
  Shoebox box;
  box.panel = 0;
  box.bbox[0] = x;
  box.bbox[1] = x + 1;
  box.bbox[2] = y;
  box.bbox[3] = y + 1;
  box.bbox[4] = 0;
  box.bbox[5] = 1;
  box.data = {count};
  box.mask = {shoebox_mask::kValid | shoebox_mask::kForeground};
  box.background = {0.0f};
  return box;
}

Vec3 ray_through(const Experiment &e, double px_fast, double px_slow) {
  const Panel &p = e.detector[0];
  const Vec3 lab = p.lab_coord_mm(px_fast * p.pixel_size[0], px_slow * p.pixel_size[1]);
  return lab / lab.norm();
}

}  // namespace

TEST(a_spot_on_its_own_beam_has_no_angular_spread) {
  // One pixel, and the beam through its centre: the only angle in the sum is
  // zero. If this is not zero the pixel centre convention is wrong, and every
  // other number here would be wrong with it by a constant.
  const Experiment e = simple_experiment(200.0);
  const Shoebox box = one_pixel(500, 500, 100.0f);
  const Vec3 s1 = ray_through(e, 500.5, 500.5);
  double variance = -1.0;
  check::is_true(spot_angular_variance(e, box, s1, &variance), "has a variance");
  check::close(variance, 0.0, 1e-24, "and it is zero");
}

TEST(the_spread_of_two_pixels_is_the_angle_between_them) {
  // Two pixels of equal weight either side of the beam. The weighted mean of
  // the squared angles is the square of the half-separation, and the sample
  // variance divides by (w - 1), so the answer is known exactly.
  const Experiment e = simple_experiment(200.0);
  Shoebox box;
  box.panel = 0;
  box.bbox[0] = 500;
  box.bbox[1] = 502;
  box.bbox[2] = 500;
  box.bbox[3] = 501;
  box.bbox[4] = 0;
  box.bbox[5] = 1;
  box.data = {50.0f, 50.0f};
  box.mask = {static_cast<std::uint8_t>(shoebox_mask::kValid | shoebox_mask::kForeground),
              static_cast<std::uint8_t>(shoebox_mask::kValid | shoebox_mask::kForeground)};
  box.background = {0.0f, 0.0f};

  const Vec3 s1 = ray_through(e, 501.0, 500.5);  // midway between the two
  // The two angles are computed rather than assumed equal: a pixel half a
  // pixel to the left and one half a pixel to the right do not subtend the
  // same angle, because the angle is not linear in position. Assuming they did
  // put this test 2.5e-7 out, which is small enough to look like rounding and
  // is not.
  const Vec3 reference = s1 / s1.norm();
  const double left =
      std::acos(std::fmin(1.0, ray_through(e, 500.5, 500.5).dot(reference)));
  const double right =
      std::acos(std::fmin(1.0, ray_through(e, 501.5, 500.5).dot(reference)));

  double variance = 0.0;
  check::is_true(spot_angular_variance(e, box, s1, &variance), "has a variance");
  // 50 counts on each, over (100 - 1).
  check::close(variance, (50.0 * left * left + 50.0 * right * right) / 99.0, 1e-20,
               "known exactly");
}

TEST(the_variance_is_weighted_by_the_counts) {
  // Moving weight onto the pixel nearer the beam must lower the spread. A sum
  // that ignored the counts would give the same answer for both.
  const Experiment e = simple_experiment(200.0);
  const auto spread = [&](float near_count, float far_count) {
    Shoebox box;
    box.panel = 0;
    box.bbox[0] = 500;
    box.bbox[1] = 503;
    box.bbox[2] = 500;
    box.bbox[3] = 501;
    box.bbox[4] = 0;
    box.bbox[5] = 1;
    box.data = {near_count, 0.0f, far_count};
    const std::uint8_t on = shoebox_mask::kValid | shoebox_mask::kForeground;
    box.mask = {on, on, on};
    box.background = {0.0f, 0.0f, 0.0f};
    double variance = 0.0;
    spot_angular_variance(e, box, ray_through(e, 500.5, 500.5), &variance);
    return variance;
  };
  check::is_true(spread(90.0f, 10.0f) < spread(50.0f, 50.0f),
                 "weight near the beam narrows it");
  check::is_true(spread(10.0f, 90.0f) > spread(50.0f, 50.0f),
                 "and weight away from it widens it");
}

TEST(background_is_subtracted_before_the_spread_is_measured) {
  // Kabsch section 3.1 says to subtract the background first. A flat pedestal
  // under the peak would otherwise pull the spread towards the box edges,
  // which is where most of the pixels are.
  const Experiment e = simple_experiment(200.0);
  Shoebox box;
  box.panel = 0;
  box.bbox[0] = 499;
  box.bbox[1] = 502;
  box.bbox[2] = 500;
  box.bbox[3] = 501;
  box.bbox[4] = 0;
  box.bbox[5] = 1;
  const std::uint8_t on = shoebox_mask::kValid | shoebox_mask::kForeground;
  box.mask = {on, on, on};
  const Vec3 s1 = ray_through(e, 500.5, 500.5);

  box.data = {0.0f, 100.0f, 0.0f};
  box.background = {0.0f, 0.0f, 0.0f};
  double clean = 0.0;
  spot_angular_variance(e, box, s1, &clean);

  box.data = {20.0f, 120.0f, 20.0f};
  box.background = {20.0f, 20.0f, 20.0f};
  double pedestal = 0.0;
  spot_angular_variance(e, box, s1, &pedestal);
  check::close(pedestal, clean, 1e-20, "the pedestal must not widen the spot");
}

TEST(a_spot_with_one_count_has_no_variance_and_is_skipped) {
  // The sample variance of a single observation is undefined, not zero.
  // Returning zero would drag the mean down by however many such spots there
  // are, silently.
  const Experiment e = simple_experiment(200.0);
  const Vec3 beam = ray_through(e, 500.5, 500.5);
  double variance = 999.0;
  check::is_true(
      !spot_angular_variance(e, one_pixel(500, 500, 1.0f), beam, &variance),
      "one count is refused");
  check::is_true(
      !spot_angular_variance(e, one_pixel(500, 500, 0.0f), beam, &variance),
      "and so is none");

  // And such spots do not enter the mean at all.
  std::vector<Shoebox> boxes = {one_pixel(500, 500, 1.0f)};
  std::vector<Vec3> s1 = {ray_through(e, 500.5, 500.5)};
  std::size_t used = 99;
  const double sigma = beam_divergence(e, boxes, s1, &used);
  check::equal(static_cast<long long>(used), 0, "none used");
  check::close(sigma, 0.0, 0.0, "and no answer invented");
}

TEST(sigma_d_is_the_root_mean_of_the_variances_in_degrees) {
  const Experiment e = simple_experiment(200.0);
  std::vector<Shoebox> boxes;
  std::vector<Vec3> s1;
  // Three identical spots, so the root mean is just the one value.
  for (int i = 0; i < 3; ++i) {
    Shoebox box;
    box.panel = 0;
    box.bbox[0] = 500;
    box.bbox[1] = 502;
    box.bbox[2] = 500;
    box.bbox[3] = 501;
    box.bbox[4] = 0;
    box.bbox[5] = 1;
    box.data = {50.0f, 50.0f};
    const std::uint8_t on = shoebox_mask::kValid | shoebox_mask::kForeground;
    box.mask = {on, on};
    box.background = {0.0f, 0.0f};
    boxes.push_back(box);
    s1.push_back(ray_through(e, 501.0, 500.5));
  }
  double one = 0.0;
  spot_angular_variance(e, boxes[0], s1[0], &one);
  std::size_t used = 0;
  const double sigma = beam_divergence(e, boxes, s1, &used);
  check::equal(static_cast<long long>(used), 3, "all three used");
  check::close(sigma, Scan::degrees(std::sqrt(one)), 1e-15, "root mean, in degrees");
}
