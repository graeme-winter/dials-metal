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
  ReferenceProfiles reference = make_reference(spec, 1, 1, 1, 100);
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
  ReferenceProfiles reference = make_reference(spec, 3, 1, 1, 100);
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

TEST(the_profile_varies_smoothly_across_a_cell_boundary) {
  // Taking the nearest cell's profile makes the model jump at a boundary, so
  // two reflections a pixel apart either side of one are fitted with different
  // profiles and their intensities differ by more than their positions
  // warrant. The weights fall linearly with distance, as Leslie section 6.1
  // describes, so the profile has to change continuously instead.
  const GridSpec spec = small_spec();
  Panel panel;
  panel.image_size[0] = 900;
  panel.image_size[1] = 900;
  ReferenceProfiles reference = make_reference(spec, 3, 1, 1, 100);

  // Two neighbouring cells with visibly different profiles.
  for (std::size_t region = 0; region < reference.profile.size(); ++region) {
    Transformed t = planted_grid(spec, 100.0, 0.0);
    if (region % 3 == 1) {
      std::rotate(t.data.begin(), t.data.begin() + 1, t.data.end());
    }
    for (int i = 0; i < 20; ++i) add_reference(&reference, region, t);
  }
  finalise_reference(&reference, 5);

  // Walk across the boundary between the first and second cells in fast.
  const double boundary = 300.0;
  std::vector<double> before =
      profile_at(reference, panel, 0, boundary - 1.0, 450.0, 50.0);
  std::vector<double> after =
      profile_at(reference, panel, 0, boundary + 1.0, 450.0, 50.0);
  // Against the profile's OWN scale, not an absolute number. These profiles
  // are normalised over 343 grid points so their peak is 0.0198, and an
  // absolute threshold of 0.01 asks whether the difference is half the peak --
  // which was this test's first version, and it failed on correct code.
  double peak = 0.0;
  for (double v : before) peak = std::max(peak, v);
  check::is_true(peak > 0.0, "there is a profile at all");

  double jump = 0.0;
  for (std::size_t i = 0; i < before.size(); ++i) {
    jump = std::max(jump, std::fabs(before[i] - after[i]));
  }
  check::is_true(jump < 0.05 * peak,
                 "the profile does not jump at a cell boundary");

  // And it really does change between the cell centres, or the interpolation
  // has simply flattened everything into one profile.
  const std::vector<double> at_first = profile_at(reference, panel, 0, 150.0, 450.0, 50.0);
  const std::vector<double> at_second = profile_at(reference, panel, 0, 450.0, 450.0, 50.0);
  double change = 0.0;
  for (std::size_t i = 0; i < at_first.size(); ++i) {
    change = std::max(change, std::fabs(at_first[i] - at_second[i]));
  }
  check::is_true(change > 0.2 * peak, "but it does change between cell centres");
}

TEST(the_weights_of_the_nearby_cells_sum_to_one) {
  // Including at the corners, where fewer cells exist than the interpolation
  // reaches for: the weights must fall back onto the cells that are there
  // rather than quietly summing to less and scaling every profile down.
  const GridSpec spec = small_spec();
  Panel panel;
  panel.image_size[0] = 900;
  panel.image_size[1] = 900;
  const ReferenceProfiles reference = make_reference(spec, 3, 4, 1, 400);
  for (double x : {0.0, 1.0, 150.0, 449.0, 450.0, 899.0}) {
    for (double y : {0.0, 450.0, 899.0}) {
      for (double z : {0.0, 50.0, 200.0, 399.0}) {
        const std::vector<Neighbour> near =
            neighbours_of(reference, panel, 0, x, y, z);
        double total = 0.0;
        for (const Neighbour &n : near) {
          total += n.weight;
          check::is_true(n.region < reference.region_count(), "a real cell");
        }
        check::close(total, 1.0, 1e-12, "the weights sum to one");
      }
    }
  }
}

TEST(the_scan_is_divided_as_well_as_the_detector) {
  // The profile changes along the scan because the crystal does, which is why
  // it is refined scan-varying. With one profile for a whole scan the fitted
  // intensities drifted from 0.986 of DIALS' at the start to 0.951 at the end.
  const GridSpec spec = small_spec();
  Panel panel;
  panel.image_size[0] = 900;
  panel.image_size[1] = 900;
  ReferenceProfiles reference = make_reference(spec, 1, 4, 1, 400);
  check::equal(static_cast<long long>(reference.region_count()), 4,
               "one detector region, four scan blocks");

  // A different profile in the first block and the last.
  for (int i = 0; i < 20; ++i) {
    add_reference(&reference, reference.region_of(panel, 0, 450.0, 450.0, 10.0),
                  planted_grid(spec, 100.0, 0.0));
    Transformed late = planted_grid(spec, 100.0, 0.0);
    std::rotate(late.data.begin(), late.data.begin() + 2, late.data.end());
    add_reference(&reference,
                  reference.region_of(panel, 0, 450.0, 450.0, 390.0), late);
  }
  finalise_reference(&reference, 5);

  const std::vector<double> early = profile_at(reference, panel, 0, 450.0, 450.0, 10.0);
  const std::vector<double> late = profile_at(reference, panel, 0, 450.0, 450.0, 390.0);
  double peak = 0.0;
  for (double v : early) peak = std::max(peak, v);
  double change = 0.0;
  for (std::size_t i = 0; i < early.size(); ++i) {
    change = std::max(change, std::fabs(early[i] - late[i]));
  }
  check::is_true(change > 0.2 * peak,
                 "the profile at the start of the scan is not the one at the end");
}
