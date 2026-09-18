// mxi_profile: estimate the Gaussian profile model from indexed strong spots.
//
//   mxi_profile refined.expt refined.refl
//
// The reflection table must carry its shoeboxes, which is what dials.find_spots
// writes and what this package now preserves. A table whose shoeboxes have been
// stripped cannot be used: the estimate is made from the pixels.

#include <cmath>
#include <cstdio>
#include <set>
#include <string>
#include <vector>

#include "args.h"
#include "expt.h"
#include "profile_model.h"
#include "refl.h"
#include "shoebox.h"

using namespace mxi;

namespace {

void usage(const char *program) {
  std::printf(
      "usage: %s [options] EXPT REFL\n"
      "\n"
      "Estimate the Gaussian profile model from the strong spots' shoeboxes.\n"
      "\n"
      "  --all             use every reflection, not only the strong ones\n"
      "  --n-sigma N       width of the integration region, in sigmas (3)\n"
      "  --max-centroid-shift N  reject a spot whose observed centroid is more\n"
      "                    than N images from its predicted angle; 0 keeps all (5)\n",
      program);
}

}  // namespace

int main(int argc, char **argv) {
  const std::set<std::string> known = {"--all", "--n-sigma", "--max-centroid-shift"};
  const std::set<std::string> takes_value = {"--n-sigma", "--max-centroid-shift"};
  const Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (args.help) {
    usage(argv[0]);
    return 0;
  }
  if (!args.ok) {
    std::fprintf(stderr, "mxi_profile: %s\n", args.error.c_str());
    return 2;
  }
  if (args.positional.size() != 2) {
    std::fprintf(stderr, "mxi_profile: expected an .expt and a .refl\n");
    usage(argv[0]);
    return 2;
  }

  try {
    const ExperimentList experiments = read_experiments(args.positional[0]);
    const Table reflections = read_reflections(args.positional[1]);
    if (experiments.size() == 0) {
      std::fprintf(stderr, "mxi_profile: no experiments\n");
      return 1;
    }

    const std::vector<Shoebox> boxes = decode_shoeboxes(reflections);
    if (boxes.empty()) {
      std::fprintf(stderr,
                   "mxi_profile: %s has no shoeboxes, and the profile model is "
                   "estimated from the pixels. Use a table that still has "
                   "them; dials.find_spots writes them and this package keeps "
                   "them.\n",
                   args.positional[1].c_str());
      return 1;
    }
    if (!reflections.has("s1")) {
      std::fprintf(stderr,
                   "mxi_profile: %s has no s1 column, so there is nothing to "
                   "measure the spread of each spot about\n",
                   args.positional[1].c_str());
      return 1;
    }

    // Strong reflections only by default, which is what the estimate is
    // defined over. Saying how many were dropped matters: a mean over a
    // different set of spots is a different number.
    const bool everything = args.has("--all");
    const double max_shift = args.number("--max-centroid-shift", 5.0);
    const Column &s = reflections.at("s1");
    const bool has_flags = reflections.has("flags");
    std::vector<Shoebox> selected;
    std::vector<Vec3> s1;
    for (std::size_t i = 0; i < reflections.nrows && i < boxes.size(); ++i) {
      if (!everything && has_flags &&
          (reflections.at("flags").integer(i) & flag::kStrong) == 0) {
        continue;
      }
      selected.push_back(boxes[i]);
      s1.push_back({s.real(i, 0), s.real(i, 1), s.real(i, 2)});
    }

    std::printf("%zu reflections, %zu with shoeboxes, %zu %s\n", reflections.nrows,
                boxes.size(), selected.size(),
                everything ? "used" : "strong and used");

    std::size_t used = 0;
    ProfileModel model;
    model.n_sigma = args.number("--n-sigma", 3.0);
    model.sigma_d = beam_divergence(experiments[0], selected, s1, &used);
    model.n_used = used;

    std::printf("sigma_D (beam divergence) %.9f degrees, from %zu spots\n",
                model.sigma_d, model.n_used);
    // sigma_M, from one sample per image each spot was seen on.
    std::vector<RangeSample> samples;
    if (reflections.has("xyzcal.mm")) {
      const Column &cal = reflections.at("xyzcal.mm");
      const Vec3 s0 = experiments[0].beam.s0();
      const Vec3 axis = experiments[0].goniometer.lab_axis();
      for (std::size_t i = 0; i < reflections.nrows && i < boxes.size(); ++i) {
        if (!everything && has_flags &&
            (reflections.at("flags").integer(i) & flag::kStrong) == 0) {
          continue;
        }
        // A row with no prediction carries uninitialised xyzcal -- denormals,
        // not zeros -- and feeding those to the likelihood put the estimate
        // out by a factor of four.
        if (!has_prediction(reflections, i)) continue;
        // Kabsch step (vii): reject a spot whose observed centroid is far from
        // where the model puts it. Without this the estimate is 0.506 rather
        // than 0.293 degrees, because a handful of spots whose shoebox sits
        // hundreds of images from their predicted angle dominate a likelihood
        // that falls off quadratically.
        //
        // Kabsch says "deviates too much" without a number. Measured here, the
        // answer is flat at 0.293 for any cut between one and ten images and
        // moves only outside that, so within the plateau this is not a knob.
        if (max_shift > 0.0 && reflections.has("xyzobs.px.value") &&
            reflections.has("xyzcal.px")) {
          const double observed = reflections.at("xyzobs.px.value").real(i, 2);
          const double predicted = reflections.at("xyzcal.px").real(i, 2);
          if (std::fabs(observed - predicted) > max_shift) continue;
        }
        const Vec3 beam{s.real(i, 0), s.real(i, 1), s.real(i, 2)};
        const Vec3 difference = beam - s0;
        const double length = difference.norm();
        if (!(length > 0.0)) continue;
        const double zeta = axis.dot(difference / length);
        for (const RangeSample &sample :
             range_samples(experiments[0], boxes[i], cal.real(i, 2), zeta)) {
          samples.push_back(sample);
        }
      }
      model.sigma_m = reflecting_range(
          samples, Scan::radians(experiments[0].scan.osc_width));
      std::printf("sigma_M (reflecting range) %.9f degrees, from %zu images\n",
                  model.sigma_m, samples.size());
    } else {
      std::printf("sigma_M needs xyzcal.mm, which this table does not have\n");
    }
    std::printf("n_sigma %.1f\n", model.n_sigma);
    // Said here rather than only in the documentation, because a number that
    // disagrees with DIALS and does not say so is worse than no number.
    std::printf(
        "\nNOTE: neither number agrees with dials.integrate yet. On 1800\n"
        "images of insulin it records sigma_b = 0.031698 and sigma_m =\n"
        "0.097667, against 0.030356 here (-4.2 per cent) and 0.293 (a factor\n"
        "of three). Both are unexplained. See docs/integration.md.\n");
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_profile: %s\n", e.what());
    return 1;
  }
}
