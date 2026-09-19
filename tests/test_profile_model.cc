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

TEST(the_background_is_not_subtracted_which_is_a_departure_from_kabsch) {
  // Kabsch section 3.1 step (v) says to subtract the background before
  // measuring the spread. DIALS does not, and carries a note in its own source
  // saying so. This follows DIALS, because standing in for dials.integrate is
  // the point, and the difference is recorded rather than quietly improved.
  //
  // It costs nothing on a table out of dials.find_spots, where the background
  // is zero. It would matter on one where it is not, and a pedestal does widen
  // the measured spread -- which is what this test now pins.
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
  check::is_true(pedestal > clean,
                 "with the background left in, a pedestal widens the spot");
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

// --------------------------------------------------------------------------
// the reflecting range
// --------------------------------------------------------------------------

TEST(the_recorded_fraction_is_a_density_in_delta) {
  // It has to integrate to one over delta, or the product of these is not a
  // likelihood and maximising it means nothing.
  const double oscillation = Scan::radians(0.1);
  const double sigma = Scan::radians(0.08);
  double total = 0.0;
  const double step = Scan::radians(0.002);
  for (double delta = -Scan::radians(3.0); delta < Scan::radians(3.0); delta += step) {
    total += recorded_fraction(delta, 0.7, sigma, oscillation) * step;
  }
  check::close(total, 1.0, 1e-6, "integrates to one");
}

TEST(a_planted_reflecting_range_comes_back) {
  // The estimator checked against a known answer rather than against DIALS.
  // Samples are drawn from the model itself by inverting its cumulative
  // distribution, so the only question is whether the fit recovers what was
  // put in.
  const double oscillation = Scan::radians(0.1);
  for (double planted_degrees : {0.05, 0.1, 0.3}) {
    const double planted = Scan::radians(planted_degrees);
    const double zeta = 0.8;
    std::vector<RangeSample> samples;
    // The density is a box of width `oscillation` convolved with a Gaussian of
    // width sigma/zeta, so a draw is the sum of a uniform and a normal.
    std::uint64_t state = 12345;
    const auto uniform = [&state]() {
      state = state * 6364136223846793005ULL + 1442695040888963407ULL;
      return static_cast<double>((state >> 11) & ((1ULL << 53) - 1)) /
             static_cast<double>(1ULL << 53);
    };
    for (int i = 0; i < 20000; ++i) {
      const double u1 = std::fmax(uniform(), 1e-12);
      const double u2 = uniform();
      const double normal = std::sqrt(-2.0 * std::log(u1)) *
                            std::cos(2.0 * 3.14159265358979323846 * u2);
      const double box = (uniform() - 0.5) * oscillation;
      samples.push_back({box + normal * planted / zeta, zeta});
    }
    const double found = reflecting_range(samples, oscillation, 0.05);
    check::close(found, planted_degrees, 0.05 * planted_degrees,
                 "the planted range comes back within five per cent");
  }
}

TEST(samples_near_the_rotation_axis_are_dropped) {
  // Their reflecting range is sigma_M / |zeta|, which diverges as zeta goes to
  // zero: they carry no information about sigma_M and would dominate the
  // likelihood if kept.
  const double oscillation = Scan::radians(0.1);
  std::vector<RangeSample> samples;
  for (int i = 0; i < 500; ++i) {
    samples.push_back({Scan::radians(0.02 * ((i % 7) - 3)), 0.8});
  }
  const double clean = reflecting_range(samples, oscillation, 0.05);
  // Add a hundred samples sitting on the axis, with wild offsets.
  for (int i = 0; i < 100; ++i) {
    samples.push_back({Scan::radians(5.0 * ((i % 3) - 1)), 0.001});
  }
  const double polluted = reflecting_range(samples, oscillation, 0.05);
  check::close(polluted, clean, 1e-9, "the cut must remove them entirely");

  // Without the cut they pull the answer DOWN, not up, which is the opposite
  // of what was assumed when this test was written: a sample whose modelled
  // range is sigma/0.001 makes every sigma look equally bad, the likelihood
  // flattens, and the fit settles lower. Measured at 0.0109 against 0.0168.
  // Either way it is wrong, and either way the point is that the cut decides
  // the answer rather than trimming it.
  const double without_cut = reflecting_range(samples, oscillation, 0.0);
  check::is_true(without_cut < 0.8 * clean, "and without it they take over");
}

TEST(one_sample_per_image_not_one_per_reflection) {
  // With one per reflection the gap is bounded by half an oscillation width by
  // construction, and the likelihood is maximised by driving sigma to zero.
  // This checks the gathering really does produce one per image with counts.
  Experiment e = simple_experiment(200.0);
  e.scan.osc_width = 0.1;
  Shoebox box;
  box.panel = 0;
  box.bbox[0] = 500;
  box.bbox[1] = 501;
  box.bbox[2] = 500;
  box.bbox[3] = 501;
  box.bbox[4] = 10;
  box.bbox[5] = 14;  // four images
  const std::uint8_t on = shoebox_mask::kValid | shoebox_mask::kForeground;
  // The criterion is the MASK, not the counts: an image counts if the spot
  // finder marked any pixel on it as valid foreground, whatever its value.
  // The last image here is marked but empty, and DIALS counts it.
  box.data = {5.0f, 20.0f, 20.0f, 0.0f};
  box.mask = {on, on, on, on};
  box.background = {0.0f, 0.0f, 0.0f, 0.0f};

  const double phi = Scan::radians(1.2);
  const std::vector<RangeSample> samples = range_samples(e, box, phi, 0.7);
  check::equal(static_cast<long long>(samples.size()), 4,
               "one per marked image, counts or no counts");

  // An image the spot finder did not mark at all is not a sample.
  box.mask = {on, on, on, 0};
  check::equal(static_cast<long long>(range_samples(e, box, phi, 0.7).size()), 3,
               "an unmarked image is not one");
  // The images are consecutive, so the gaps step by exactly one oscillation.
  check::close(samples[0].delta - samples[1].delta, Scan::radians(0.1), 1e-12,
               "consecutive images are one oscillation apart");
  check::close(samples[1].delta - samples[2].delta, Scan::radians(0.1), 1e-12,
               "and so are the next two");
}

// --------------------------------------------------------------------------
// the captured fraction
// --------------------------------------------------------------------------

TEST(the_kabsch_frame_is_orthonormal_and_e1_is_perpendicular_to_both_beams) {
  // e1 = s1 x s0 is the axis about which the point would cross the Ewald
  // sphere by the shortest route, so it is perpendicular to both beams. Using
  // the difference s1 - s0 instead, which is what this did, gives a vector
  // that is not, and put sigma_M out by fifty per cent.
  const Experiment e = simple_experiment(200.0);
  const Vec3 s1 = ray_through(e, 620.0, 540.0) * (1.0 / e.beam.wavelength);
  const KabschFrame f = kabsch_frame(e, s1);
  check::is_true(f.valid, "a frame exists");
  check::close(f.e1.norm(), 1.0, 1e-12, "e1 is a unit vector");
  check::close(f.e2.norm(), 1.0, 1e-12, "e2 is a unit vector");
  check::close(f.e3.norm(), 1.0, 1e-12, "e3 is a unit vector");
  check::close(f.e1.dot(s1), 0.0, 1e-12, "e1 is perpendicular to s1");
  check::close(f.e1.dot(e.beam.s0()), 0.0, 1e-12, "and to s0");
  check::close(f.e1.dot(f.e2), 0.0, 1e-12, "e1 and e2 are orthogonal");
}

TEST(a_pixel_on_the_beam_sits_at_the_origin_of_the_frame) {
  const Experiment e = simple_experiment(200.0);
  const Vec3 s1 = ray_through(e, 620.5, 540.5) * (1.0 / e.beam.wavelength);
  const KabschFrame f = kabsch_frame(e, s1);
  const Epsilon eps =
      epsilon_of(e, f, e.detector[0], 620.5, 540.5, Scan::radians(1.0),
                 Scan::radians(1.0));
  check::close(eps.e1, 0.0, 1e-12, "eps1 is zero at the reflection");
  check::close(eps.e2, 0.0, 1e-12, "eps2 too");
  check::close(eps.e3, 0.0, 1e-12, "and eps3 when the image is at the angle");
}

TEST(eps3_is_the_rotation_offset_scaled_by_zeta) {
  // Not the offset itself. The scaling is what makes a reflection near the
  // rotation axis, which sweeps through the sphere slowly, comparable with one
  // far from it.
  const Experiment e = simple_experiment(200.0);
  const Vec3 s1 = ray_through(e, 620.0, 540.0) * (1.0 / e.beam.wavelength);
  const KabschFrame f = kabsch_frame(e, s1);
  const double phi = Scan::radians(1.0);
  const double offset = Scan::radians(0.03);
  const Epsilon eps = epsilon_of(e, f, e.detector[0], 620.0, 540.0, phi + offset, phi);
  check::close(eps.e3, Scan::degrees(f.zeta * offset), 1e-12, "zeta times the offset");
  check::is_true(std::abs(f.zeta) < 1.0, "and zeta is less than one here");
}

TEST(the_captured_fraction_counts_what_is_inside) {
  // Two images, one inside a sigma and one outside, so the answer is known
  // without any distributional argument: a half at one sigma and all of it at
  // two.
  Experiment e = simple_experiment(200.0);
  e.scan.osc_width = 0.05;
  const Vec3 s1 = ray_through(e, 620.5, 540.5) * (1.0 / e.beam.wavelength);
  const KabschFrame f = kabsch_frame(e, s1);

  Shoebox box;
  box.panel = 0;
  box.bbox[0] = 620;
  box.bbox[1] = 621;
  box.bbox[2] = 540;
  box.bbox[3] = 541;
  box.bbox[4] = 20;
  box.bbox[5] = 22;
  const std::uint8_t on = shoebox_mask::kValid | shoebox_mask::kForeground;
  box.data = {10.0f, 10.0f};
  box.mask = {on, on};
  box.background = {0.0f, 0.0f};

  // Put the Bragg angle at the centre of the first image, so the second is one
  // oscillation away.
  const double width = Scan::radians(e.scan.osc_width);
  const double phi = Scan::radians(e.scan.osc_start) + 20.5 * width;
  const double step = Scan::degrees(std::abs(f.zeta) * width);

  // A sigma_M just over the step puts both inside; just under puts one in.
  const double generous = 1e6;  // sigma_D large enough that the detector never cuts
  const Capture wide = capture_fractions(e, {box}, {s1}, {phi}, generous, step * 1.01);
  check::close(wide.rotation[0], 1.0, 1e-12, "both images inside one sigma");

  const Capture narrow = capture_fractions(e, {box}, {s1}, {phi}, generous, step * 0.99);
  check::close(narrow.rotation[0], 0.5, 1e-12, "only the nearer one");
  check::close(narrow.rotation[1], 1.0, 1e-12, "and both by two sigma");
}

TEST(the_detector_and_rotation_directions_are_reported_apart) {
  // Which sigma is wrong cannot be read off the combined number: a spot too
  // concentrated on the detector and too spread in rotation would look
  // correct. This checks a cut in one direction does not move the other.
  Experiment e = simple_experiment(200.0);
  e.scan.osc_width = 0.05;
  const Vec3 s1 = ray_through(e, 620.5, 540.5) * (1.0 / e.beam.wavelength);
  Shoebox box;
  box.panel = 0;
  box.bbox[0] = 619;
  box.bbox[1] = 622;
  box.bbox[2] = 540;
  box.bbox[3] = 541;
  box.bbox[4] = 20;
  box.bbox[5] = 21;
  const std::uint8_t on = shoebox_mask::kValid | shoebox_mask::kForeground;
  box.data = {10.0f, 10.0f, 10.0f};
  box.mask = {on, on, on};
  box.background = {0.0f, 0.0f, 0.0f};
  const double phi = Scan::radians(e.scan.osc_start) + 20.5 * Scan::radians(0.05);

  // Rotation wide open, detector tight: the rotation fraction must stay one.
  const Capture c = capture_fractions(e, {box}, {s1}, {phi}, 1e-6, 1e6);
  check::close(c.rotation[0], 1.0, 1e-12, "the rotation direction is untouched");
  check::is_true(c.detector[0] < 1.0, "while the detector direction cuts");
  check::close(c.fraction[0], c.detector[0], 1e-12,
               "and the combined figure follows the one that cuts");
}
