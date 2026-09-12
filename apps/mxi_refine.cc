// mxi_refine: refine the experimental model against indexed spot centroids.
//
//   mxi_refine indexed.expt indexed.refl

#include <cstdio>
#include <cstdlib>
#include <string>

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
      "  --macrocycles N   (3)\n"
      "  --outlier-sigma S (4; 0 disables rejection)\n"
      "  --output-expt P   (refined.expt)\n"
      "  --output-refl P   (refined.refl)\n");
}
}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    usage();
    return 2;
  }
  std::string out_expt = "refined.expt";
  std::string out_refl = "refined.refl";
  RefineOptions options;
  options.verbose = true;
  bool conditional_depth = false;
  int scan_points = 1;

  for (int i = 3; i < argc; ++i) {
    const std::string arg = argv[i];
    const auto next = [&]() -> const char * {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "mxi_refine: %s needs a value\n", arg.c_str());
        std::exit(2);
      }
      return argv[++i];
    };
    if (arg == "--no-crystal") options.crystal = false;
    else if (arg == "--no-detector") options.detector = false;
    else if (arg == "--beam") options.beam = true;
    else if (arg == "--separate") options.shared_crystal = false;
    else if (arg == "--conditional-depth") conditional_depth = true;
    else if (arg == "--scan-varying") scan_points = std::atoi(next());
    else if (arg == "--macrocycles") options.macrocycles = std::atoi(next());
    else if (arg == "--outlier-sigma") options.outlier_sigma = std::atof(next());
    else if (arg == "--output-expt") out_expt = next();
    else if (arg == "--output-refl") out_refl = next();
    else if (arg == "-h" || arg == "--help") { usage(); return 0; }
    else {
      std::fprintf(stderr, "mxi_refine: unknown option '%s'\n", arg.c_str());
      return 2;
    }
  }

  try {
    ExperimentList experiments = read_experiments(argv[1]);
    Table reflections = read_reflections(argv[2]);
    if (conditional_depth) {
      for (Experiment &e : experiments) {
        for (Panel &p : e.detector.panels) p.parallax_conditional = true;
      }
      std::printf("using the conditional absorption depth\n");
    }
    std::printf("%zu experiments, %zu reflections\n", experiments.size(),
                reflections.nrows);

    // Static first, always. Scan-varying control points started from an
    // unrefined model absorb errors that belong to the detector, and the
    // result fits well and means nothing.
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
