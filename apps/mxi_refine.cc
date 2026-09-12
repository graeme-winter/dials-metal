// mxi_refine: refine the experimental model against indexed spot centroids.
//
//   mxi_refine indexed.expt indexed.refl

#include <cstdio>
#include <cstdlib>
#include <set>
#include <string>

#include "args.h"
#include "expt.h"
#include "refine.h"
#include "refl.h"

using namespace mxi;

namespace {
void usage() {
  std::printf(
      "usage: mxi_refine INDEXED.expt INDEXED.refl [options]\n"
      "  --no-crystal      hold the crystal fixed\n"
      "  --no-detector     hold the detector fixed\n"
      "  --beam            refine the beam direction too (off: correlated\n"
      "                    with the detector on a single sweep)\n"
      "  --separate        one crystal per experiment instead of one shared\n"
      "  --conditional-depth  mean depth given absorption, not eqn (6)\n"
      "  --scan-varying N  control points in A across each scan (1 = static)\n"
      "  --unit-weights    ignore the centroid variances\n"
      "  --strong-only     build the model from the stronger half only\n"
      "  --z-weight W      scale the weight on the rotation-angle residual\n"
      "  --macrocycles N   (3)\n"
      "  --outlier-sigma S (4; 0 disables rejection)\n"
      "  --output-expt P   (refined.expt)\n"
      "  --output-refl P   (refined.refl)\n");
}
}  // namespace

int main(int argc, char **argv) {
  const std::set<std::string> known = {
      "--no-crystal",   "--no-detector",  "--beam",         "--separate",
      "--macrocycles",  "--outlier-sigma", "--output-expt", "--output-refl",
      "--conditional-depth", "--scan-varying", "--unit-weights",
      "--strong-only",  "--z-weight"};
  const std::set<std::string> takes_value = {
      "--macrocycles", "--outlier-sigma", "--output-expt", "--output-refl",
      "--scan-varying", "--z-weight"};
  const Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (args.help) {
    usage();
    return 0;
  }
  if (!args.ok) {
    std::fprintf(stderr, "mxi_refine: %s\n", args.error.c_str());
    return 2;
  }
  if (args.positional.size() != 2) {
    std::fprintf(stderr,
                 "mxi_refine: expected an .expt and a .refl, got %zu file "
                 "arguments\n",
                 args.positional.size());
    usage();
    return 2;
  }

  RefineOptions options;
  options.verbose = true;
  options.crystal = !args.has("--no-crystal");
  options.detector = !args.has("--no-detector");
  options.beam = args.has("--beam");
  options.shared_crystal = !args.has("--separate");
  options.unit_weights = args.has("--unit-weights");
  options.strong_only = args.has("--strong-only");
  options.macrocycles = static_cast<int>(args.number("--macrocycles", 3));
  options.outlier_sigma = args.number("--outlier-sigma", 4.0);
  options.z_weight = args.number("--z-weight", 1.0);
  const int scan_points = static_cast<int>(args.number("--scan-varying", 1));
  const bool conditional_depth = args.has("--conditional-depth");
  const std::string out_expt = args.value("--output-expt", "refined.expt");
  const std::string out_refl = args.value("--output-refl", "refined.refl");

  try {
    ExperimentList experiments = read_experiments(args.positional[0]);
    Table reflections = read_reflections(args.positional[1]);
    if (conditional_depth) {
      for (Experiment &e : experiments) {
        for (Panel &p : e.detector.panels) p.parallax_conditional = true;
      }
      std::printf("using the conditional absorption depth\n");
    }
    std::printf("%zu experiments, %zu reflections\n", experiments.size(),
                reflections.nrows);

    RefineResult result = refine(experiments, reflections, options);
    if (scan_points > 1 && result.n_used > 0) {
      std::printf("static: rmsd %.4f %.4f %.4f -> now %d control points\n",
                  result.rmsd_x, result.rmsd_y, result.rmsd_z, scan_points);
      options.scan_points = static_cast<std::size_t>(scan_points);
      result = refine(experiments, reflections, options);
    }
    if (result.n_used == 0) {
      std::fprintf(stderr, "mxi_refine: nothing to refine against\n");
      return 1;
    }
    std::printf("refined on %zu reflections (%zu rejected), %d steps\n",
                result.n_used, result.n_rejected, result.iterations);
    std::printf("rmsd %.4f px  %.4f px  %.4f images\n", result.rmsd_x,
                result.rmsd_y, result.rmsd_z);
    for (std::size_t i = 0; i < experiments.size(); ++i) {
      if (!experiments[i].crystal) continue;
      const UnitCell c = experiments[i].crystal->cell();
      std::printf("  [%zu] cell %.4f %.4f %.4f  %.3f %.3f %.3f  V %.1f\n", i, c.a,
                  c.b, c.c, c.alpha, c.beta, c.gamma, c.volume());
    }

    update_predictions(experiments, reflections);
    write_experiments(out_expt, experiments);
    write_reflections(out_refl, reflections);
    std::printf("wrote %s and %s\n", out_expt.c_str(), out_refl.c_str());
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_refine: %s\n", e.what());
    return 1;
  }
}
