// Reference profiles and profile fitting.

#include <cmath>
#include <vector>

#include "../src/reference.h"
#include "check.h"

using namespace mxi;

namespace {

GridSpec small_spec() {
  GridSpec spec;
  spec.n = 3;  // side 7
  spec.sigma_d = 0.02;
  spec.sigma_m = 0.07;
  spec.half_width = 3.0;
  spec.subdivisions = 2;
  return spec;
}

//: A grid holding a Gaussian bump of known total, with a flat background and
//: full coverage, so a fit has an arithmetic answer.
Transformed planted_grid(const GridSpec &spec, double intensity,
                         double background) {
  Transformed t;
  t.valid = true;
  t.data.assign(spec.size(), 0.0);
  t.background.assign(spec.size(), background);
  t.coverage.assign(spec.size(), 1.0);
  const int side = spec.side();
  const double middle = 0.5 * static_cast<double>(side - 1);
  double sum = 0.0;
  std::vector<double> shape(spec.size(), 0.0);
  for (int i3 = 0; i3 < side; ++i3) {
    for (int i2 = 0; i2 < side; ++i2) {
      for (int i1 = 0; i1 < side; ++i1) {
        const double d1 = (i1 - middle) / 1.5;
        const double d2 = (i2 - middle) / 1.5;
        const double d3 = (i3 - middle) / 1.5;
        const double v = std::exp(-0.5 * (d1 * d1 + d2 * d2 + d3 * d3));
        shape[spec.at(i1, i2, i3)] = v;
        sum += v;
      }
    }
  }
  for (std::size_t i = 0; i < shape.size(); ++i) {
    t.data[i] = background + intensity * shape[i] / sum;
  }
  return t;
}

std::vector<double> unit_profile(const GridSpec &spec) {
  const Transformed t = planted_grid(spec, 1.0, 0.0);
  std::vector<double> p(spec.size(), 0.0);
  double sum = 0.0;
  for (std::size_t i = 0; i < p.size(); ++i) {
    p[i] = t.data[i];
    sum += p[i];
  }
  for (double &v : p) v /= sum;
  return p;
}

}  // namespace

TEST(fitting_a_profile_to_itself_returns_the_intensity) {
  // The profile is exactly the shape of the data, so the scale is the answer
  // and there is nothing statistical about it.
  const GridSpec spec = small_spec();
  const std::vector<double> reference = unit_profile(spec);
  for (double intensity : {10.0, 1000.0, 100000.0}) {
    const Transformed t = planted_grid(spec, intensity, 2.0);
    const ProfileFit fit = fit_profile(reference, t);
    check::is_true(fit.valid, "fitted");
    check::close(fit.intensity, intensity, 1e-6 * intensity,
                 "the planted intensity comes back");
    check::close(fit.correlation, 1.0, 1e-9,
                 "and the profile describes it perfectly");
  }
}

TEST(a_fitted_intensity_beats_a_summed_one_on_a_weak_reflection) {
  // Leslie section 6.6: weighting by the profile reduces the variance by
  // sum(P^2) m / (sum P)^2, about two for a typical spot. This is the reason
  // profile fitting exists, so it is worth asserting rather than assuming.
  const GridSpec spec = small_spec();
  const std::vector<double> reference = unit_profile(spec);
  const double background = 3.0;
  const Transformed t = planted_grid(spec, 20.0, background);
  const ProfileFit fit = fit_profile(reference, t);
  check::is_true(fit.valid, "fitted");

  // The summed variance over the same points: every point counted equally.
  double summed_variance = 0.0;
  for (std::size_t i = 0; i < t.data.size(); ++i) summed_variance += t.data[i];
  check::is_true(fit.variance < summed_variance,
                 "the fitted variance is smaller than the summed one");

  // And the predicted ratio, which depends only on the profile's shape.
  double sp = 0.0, spp = 0.0, n = 0.0;
  for (double v : reference) {
    sp += v;
    spp += v * v;
    n += 1.0;
  }
  const double expected = spp * n / (sp * sp);
  check::is_true(expected > 1.5 && expected < 20.0,
                 "the shape predicts a real gain");
}

TEST(the_fit_reduces_to_the_sum_when_the_background_is_negligible) {
  // Leslie section 6.4: for a strong reflection with correct weights the
  // profile-fitted intensity becomes the summed one. A check on the weighting
  // rather than a remark: if the weights were wrong this would not hold.
  const GridSpec spec = small_spec();
  const std::vector<double> reference = unit_profile(spec);
  const Transformed t = planted_grid(spec, 1e6, 1e-6);
  const ProfileFit fit = fit_profile(reference, t);
  double summed = 0.0;
  for (std::size_t i = 0; i < t.data.size(); ++i) {
    summed += t.data[i] - t.background[i];
  }
  check::close(fit.intensity, summed, 1e-4 * summed,
               "fitted equals summed for a strong reflection");
}

TEST(a_grid_point_no_pixel_reached_is_not_a_zero) {
  // Coverage of zero means no measurement. Treating it as an observed zero
  // pulls every intensity down, and the points outside a spot's own footprint
  // are exactly where a profile has its tails.
  const GridSpec spec = small_spec();
  const std::vector<double> reference = unit_profile(spec);
  Transformed full = planted_grid(spec, 500.0, 1.0);
  const ProfileFit whole = fit_profile(reference, full);

  // Blank the outermost plane, as a module gap or a box edge would.
  Transformed gapped = full;
  const int side = spec.side();
  for (int i2 = 0; i2 < side; ++i2) {
    for (int i1 = 0; i1 < side; ++i1) {
      const std::size_t at = spec.at(i1, i2, 0);
      gapped.coverage[at] = 0.0;
      gapped.data[at] = 0.0;
    }
  }
  const ProfileFit partial = fit_profile(reference, gapped);
  check::is_true(partial.valid, "still fitted");
  // The profile knows what belongs in the missing plane, so the intensity
  // should barely move -- which is the other reason to fit rather than sum.
  check::close(partial.intensity, whole.intensity, 0.05 * whole.intensity,
               "a missing plane is filled in by the profile, not lost");
}

TEST(the_profile_is_the_average_shape_not_the_average_spot) {
  // Contributions are normalised before they are added, so one strong
  // reflection does not outvote a hundred weak ones.
  const GridSpec spec = small_spec();
  ReferenceProfiles reference = make_reference(spec, 1, 1);
  // Ninety-nine weak spots of one shape.
  for (int i = 0; i < 99; ++i) {
    check::is_true(add_reference(&reference, 0, planted_grid(spec, 10.0, 0.0)),
                   "added a weak one");
  }
  // And one enormous one, shifted so its shape is different.
  Transformed loud = planted_grid(spec, 1e7, 0.0);
  std::rotate(loud.data.begin(), loud.data.begin() + 1, loud.data.end());
  check::is_true(add_reference(&reference, 0, loud), "added the loud one");
  finalise_reference(&reference, 1);

  // The average must look like the weak ones, which are ninety-nine per cent
  // of the votes, not like the loud one.
  const std::vector<double> weak = unit_profile(spec);
  double worst = 0.0;
  for (std::size_t i = 0; i < weak.size(); ++i) {
    worst = std::max(worst, std::fabs(reference.profile[0][i] - weak[i]));
  }
  check::is_true(worst < 0.02,
                 "the loud spot did not take the profile over");

  double sum = 0.0;
  for (double v : reference.profile[0]) sum += v;
  check::close(sum, 1.0, 1e-9, "and the profile is normalised");
}

TEST(a_region_with_too_few_spots_borrows_the_whole_detector_average) {
  const GridSpec spec = small_spec();
  ReferenceProfiles reference = make_reference(spec, 3, 1);
  for (int i = 0; i < 40; ++i) {
    add_reference(&reference, 4, planted_grid(spec, 100.0, 0.0));
  }
  finalise_reference(&reference, 10);
  double empty_sum = 0.0;
  for (double v : reference.profile[0]) empty_sum += v;
  check::close(empty_sum, 1.0, 1e-9,
               "an empty region still has a usable profile");
  check::close(reference.profile[0][spec.at(3, 3, 3)],
               reference.profile[4][spec.at(3, 3, 3)], 1e-9,
               "and it is the average of what was seen");
}
