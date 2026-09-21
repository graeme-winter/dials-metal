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

TEST(the_lp_factor_collapses_to_the_unpolarized_form) {
  // P = (1-p) + (2p-1)(u.n)^2 + p(u.s0hat)^2, which at p = 0.5 must become the
  // textbook (1 + cos^2 2theta)/2 with no dependence on the polarization
  // normal at all. That limit is the check that the form is right rather than
  // merely fitted to a column: the three coefficients were solved for by least
  // squares and came back as (1-p), (2p-1) and p, and this is what says those
  // are the physics and not three free numbers.
  Beam beam;
  beam.direction = {0.0, 0.0, 1.0};
  beam.wavelength = 1.0;
  beam.polarization_fraction = 0.5;
  Goniometer goniometer;
  goniometer.axis = {1.0, 0.0, 0.0};

  const Vec3 s0 = beam.s0();
  for (double angle : {0.1, 0.4, 0.9}) {
    // A diffracted beam at 2theta = angle, in the plane containing the axis.
    const Vec3 s1 =
        (s0 * std::cos(angle) + Vec3{0.0, 1.0, 0.0} * std::sin(angle) * s0.norm())
            .normalized() *
        s0.norm();
    const double cos_two_theta = (s1 / s1.norm()).dot(s0 / s0.norm());
    const double lorentz =
        std::fabs(s1.dot(goniometer.lab_axis().cross(s0))) /
        (s1.norm() * s0.norm());
    const double expected =
        lorentz / (0.5 * (1.0 + cos_two_theta * cos_two_theta));
    check::close(lorentz_polarization(beam, goniometer, s1), expected,
                 1e-12 * std::fabs(expected), "unpolarized at p = 0.5");

    // And with p = 0.5 the polarization normal cannot matter.
    Beam turned = beam;
    turned.polarization_normal = {0.0, 0.0, 1.0};
    check::close(lorentz_polarization(turned, goniometer, s1),
                 lorentz_polarization(beam, goniometer, s1), 1e-12,
                 "and the normal makes no difference");
  }
}

TEST(a_polarized_beam_makes_the_normal_matter) {
  // The opposite limit: at p near one the factor must depend on where the
  // diffracted beam sits relative to the polarization plane, or the column is
  // constant and the correction does nothing.
  Beam beam;
  beam.direction = {0.0, 0.0, 1.0};
  beam.wavelength = 1.0;
  beam.polarization_fraction = 0.999;
  Goniometer goniometer;
  goniometer.axis = {1.0, 0.0, 0.0};
  const Vec3 s0 = beam.s0();

  const Vec3 along = (s0 * std::cos(0.6) +
                      Vec3{0.0, 1.0, 0.0} * std::sin(0.6) * s0.norm())
                         .normalized() * s0.norm();
  Beam turned = beam;
  turned.polarization_normal = {0.0, 0.0, 1.0};
  check::is_true(std::fabs(lorentz_polarization(beam, goniometer, along) -
                           lorentz_polarization(turned, goniometer, along)) >
                     1e-6,
                 "a polarized beam notices where the normal points");
}

TEST(the_observed_centroid_is_where_the_signal_is) {
  // dials.scale wants xyzobs.px.value, and it is the only thing in the output
  // that says where the spot actually was rather than where the model put it.
  //
  // A single lit pixel puts the answer at that pixel's centre, which is
  // arithmetic rather than statistics; the half is the pixel centre convention
  // and getting it wrong would bias every centroid in the dataset by half a
  // pixel in the same direction.
  Shoebox box = planted(2.0, 0.0, 25);
  for (std::size_t i = 0; i < box.size(); ++i) {
    if (box.mask[i] & shoebox_mask::kForeground) box.data[i] = 2.0f;
  }
  const std::size_t lit = box.at(7, 8, 1);
  box.data[lit] = 500.0f;
  IntegrateOptions options;
  options.background.tuning = 1e6;
  const IntegratedReflection r = integrate_shoebox(&box, options);
  check::is_true(r.valid && r.centroid_valid, "there is a centroid");
  check::close(r.centroid_fast, box.bbox[0] + 7 + 0.5, 1e-9, "fast");
  check::close(r.centroid_slow, box.bbox[2] + 8 + 0.5, 1e-9, "slow");
  check::close(r.centroid_z, box.bbox[4] + 1 + 0.5, 1e-9, "frame");
  // One pixel has no spread, so the variance of its mean is zero.
  check::close(r.centroid_variance_fast, 0.0, 1e-12, "and no variance");
}

TEST(the_centroid_ignores_the_background_and_the_rim) {
  // Two lit pixels either side of centre weight it to the middle; a bright
  // pixel in the BACKGROUND rim must not pull it, because the rim is not the
  // spot. Taking the centroid over the whole box instead of the foreground put
  // it 0.15 to 0.19 pixels away from DIALS' on real data.
  Shoebox box = planted(1.0, 0.0, 25);
  for (std::size_t i = 0; i < box.size(); ++i) {
    if (box.mask[i] & shoebox_mask::kForeground) box.data[i] = 1.0f;
  }
  box.data[box.at(6, 7, 1)] = 101.0f;
  box.data[box.at(8, 7, 1)] = 101.0f;
  IntegrateOptions options;
  options.background.tuning = 1e6;
  const IntegratedReflection clean = integrate_shoebox(&box, options);
  check::is_true(clean.centroid_valid, "a centroid");
  check::close(clean.centroid_fast, box.bbox[0] + 7.5, 1e-9,
               "midway between the two");

  // Now a hot pixel in the rim, far from the spot.
  Shoebox dirty = box;
  const std::size_t rim = dirty.at(1, 1, 1);
  check::is_true((dirty.mask[rim] & shoebox_mask::kBackground) != 0,
                 "that voxel really is background");
  dirty.data[rim] = 5000.0f;
  const IntegratedReflection spoiled = integrate_shoebox(&dirty, options);
  check::close(spoiled.centroid_fast, clean.centroid_fast, 1e-9,
               "the rim does not move the centroid");
}

TEST(the_resolution_column_comes_from_the_static_cell) {
  // d = 1 / |A h|, from the static cell and the Miller index. Three nearly
  // equal numbers were candidates and only one is the column: against a DIALS
  // integrated.refl the static cell agrees to 1.2e-15 for every reflection,
  // while 1/|s1 - s0| and the scan-varying A are both 2.6e-4 out. A tolerance
  // of a part in a thousand would have accepted any of them.
  const double a = 78.0;
  const Crystal cubic = Crystal::from_real_space({a, 0.0, 0.0}, {0.0, a, 0.0},
                                                 {0.0, 0.0, a});
  check::close(resolution(cubic, 1, 0, 0), a, 1e-9, "the 100 is the cell edge");
  check::close(resolution(cubic, 2, 0, 0), a / 2.0, 1e-9, "the 200 is half it");
  check::close(resolution(cubic, 1, 1, 0), a / std::sqrt(2.0), 1e-9, "the 110");
  check::close(resolution(cubic, 1, 1, 1), a / std::sqrt(3.0), 1e-9, "the 111");
  check::close(resolution(cubic, 0, 0, 0), 0.0, 0.0, "and 000 is refused");
}

TEST(partiality_is_the_part_of_the_rocking_curve_inside_the_box) {
  // The box, not the scan: what was summed is what is in the box. A box
  // spanning plus and minus three sigma of the rocking curve holds 0.9973 of
  // it, not all of it, and dials.scale divides by this.
  Scan scan;
  scan.first_image = 1;
  scan.last_image = 1000;
  scan.osc_start = 0.0;
  scan.osc_width = 0.1;
  const double sigma_m = 0.1;
  const double zeta = 1.0;  // so the rocking curve is sigma_m wide in phi

  // A box of plus and minus three sigma about a reflection in the middle.
  const double phi = scan.phi_from_z(500.0);
  const std::int32_t half = 3;  // 3 images = 0.3 degrees = 3 sigma
  const double full =
      partiality(scan, phi, zeta, sigma_m, 500 - half, 500 + half);
  check::close(full, 0.9973, 0.001, "three sigma holds 0.9973 of the curve");

  // One sigma holds much less, and the difference is not a rounding error.
  const double narrow = partiality(scan, phi, zeta, sigma_m, 499, 501);
  check::close(narrow, 0.6827, 0.01, "one sigma holds 0.6827");

  // A reflection whose box is cut off at the start of the scan is partial.
  const double edge = partiality(scan, scan.phi_from_z(0.0), zeta, sigma_m, 0, 3);
  check::is_true(edge > 0.45 && edge < 0.55,
                 "half the curve when the box starts at the reflection");

  // And it is never outside [0, 1], whatever it is handed.
  check::is_true(partiality(scan, phi, 1e-9, sigma_m, 0, 1000) <= 1.0,
                 "a diverging rocking curve does not exceed one");
}
