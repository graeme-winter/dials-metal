// mxi_profile: estimate the Gaussian profile model from indexed strong spots.
//
//   mxi_profile refined.expt refined.refl
//
// The reflection table must carry its shoeboxes, which is what dials.find_spots
// writes and what this package now preserves. A table whose shoeboxes have been
// stripped cannot be used: the estimate is made from the pixels.

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
      "  --n-sigma N       width of the integration region, in sigmas (3)\n",
      program);
}

}  // namespace

int main(int argc, char **argv) {
  const std::set<std::string> known = {"--all", "--n-sigma"};
  const std::set<std::string> takes_value = {"--n-sigma"};
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
    std::printf("sigma_M (reflecting range) not implemented\n");
    std::printf("n_sigma %.1f\n", model.n_sigma);
    // Said here rather than only in the documentation, because a number that
    // disagrees with DIALS and does not say so is worse than no number.
    std::printf(
        "\nNOTE: on 1800 images of insulin this gives 0.030356 where\n"
        "dials.integrate records sigma_b = 0.031698, a difference of 4.2 per\n"
        "cent that is not yet understood. See docs/integration.md.\n");
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_profile: %s\n", e.what());
    return 1;
  }
}
