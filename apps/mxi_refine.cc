// mxi_refine: refine the experimental model against indexed spot centroids.
//
//   mxi_refine indexed.expt indexed.refl

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <string>

#include "args.hh"
#include "../src/expt.hh"
#include "refine.hh"
#include "../src/refl.hh"

namespace mxi {

namespace {
void usage() {
  std::printf(
      "usage: mxi_refine INDEXED.expt INDEXED.refl [options]\n"
      "  --no-crystal      hold the crystal fixed\n"
      "  --jacobian-threads N  threads for the Jacobian; 0 is one per core (0)\n"
      "  --no-detector     hold the detector fixed\n"
      "  --beam            refine the beam direction too (off: correlated\n"
      "                    with the detector on a single sweep)\n"
      "  --separate        one crystal per experiment instead of one shared\n"
      "  --conditional-depth  mean depth given absorption, not eqn (6)\n"
      "  --scan-varying N  control points in A across each scan (1 = static)\n"
      "  --unit-weights    ignore the centroid variances\n"
      "  --strong-only     build the model from the stronger half only\n"
      "  --analytic        analytical derivatives, not finite differences\n"
      "  --timing          where the time went, by phase\n"
      "  --normal-threads N  threads for the normal equations; 0 is one per\n"
      "                   core (0). A reduction, so a threaded run differs from\n"
      "                   a serial one in the last bits; 1 to avoid that\n"
      "  --detector-in-scan-varying  keep refining the detector during the\n"
      "                    scan-varying pass; it is degenerate with the cell\n"
      "  --min-volume V    drop reflections whose rotation angle is not\n"
      "                    determined by the data; 0 keeps them all (0.05)\n"
      "  --z-weight W      scale the weight on the rotation-angle residual\n"
      "  --macrocycles N   (3)\n"
      "  --outlier-sigma S (4; 0 disables rejection)\n"
      "  --output-expt P   (refined.expt)\n"
      "  --output-refl P   (refined.refl)\n");
}
}  // namespace

namespace {
double now_wall() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}
}  // namespace

int run_program(int argc, char **argv) {
  const double t_start = now_wall();
  double t_read = 0.0;
  double t_write = 0.0;
  const std::set<std::string> known = {
      "--no-crystal",   "--no-detector",  "--beam",         "--separate",
      "--macrocycles",  "--outlier-sigma", "--output-expt", "--output-refl",
      "--conditional-depth", "--scan-varying", "--unit-weights",
      "--strong-only",  "--z-weight",  "--analytic", "--min-volume", "--detector-in-scan-varying", "--jacobian-threads", "--timing", "--normal-threads"};
  const std::set<std::string> takes_value = {
      "--macrocycles", "--outlier-sigma", "--output-expt", "--output-refl",
      "--scan-varying", "--z-weight", "--min-volume", "--jacobian-threads", "--normal-threads"};
  const Arguments args = parse_arguments(argc, argv, known, takes_value);
  // 0 means one per core, 1 means none. Exposed because a threading change
  // that cannot be switched off cannot be measured against its absence.
  g_normal_threads =
      static_cast<std::size_t>(args.number("--normal-threads", 0.0));
  g_jacobian_threads =
      static_cast<std::size_t>(args.number("--jacobian-threads", 0.0));
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
  options.analytic = args.has("--analytic");
  options.macrocycles = static_cast<int>(args.number("--macrocycles", 3));
  options.outlier_sigma = args.number("--outlier-sigma", 4.0);
  options.z_weight = args.number("--z-weight", 1.0);
  options.min_volume = args.number("--min-volume", 0.05);
  const int scan_points = static_cast<int>(args.number("--scan-varying", 1));
  const bool conditional_depth = args.has("--conditional-depth");
  const std::string out_expt = args.value("--output-expt", "refined.expt");
  const std::string out_refl = args.value("--output-refl", "refined.refl");

  try {
    const double t_read_start = now_wall();
    ExperimentList experiments = read_experiments(args.positional[0]);
    Table reflections = read_reflections(args.positional[1]);
    t_read = now_wall() - t_read_start;
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
      // The detector is held where the static pass put it.
      //
      // A scan-varying crystal and a refinable detector distance are
      // degenerate: the cell scales with the distance, and a crystal free at
      // every control point can pay for a smaller residual by moving the
      // detector. Measured on 1800 images of insulin, eighteen control points:
      // the distance drifts 0.27 mm and the cell volume falls 0.59 per cent,
      // for five thousandths of a pixel. Holding the detector keeps the volume
      // within 0.04 per cent of what DIALS reports and costs 0.005 px.
      //
      // There is nothing scan-varying about a detector in any case. It does
      // not move during a sweep, so letting it move while the crystal is free
      // only gives crystal drift somewhere else to go.
      if (!args.has("--detector-in-scan-varying")) {
        options.detector = false;
        std::printf(
            "holding the detector for the scan-varying pass; "
            "--detector-in-scan-varying to refine it too\n");
      }
      result = refine(experiments, reflections, options);
    }
    if (result.n_used == 0) {
      std::fprintf(stderr, "mxi_refine: nothing to refine against\n");
      return 1;
    }
    std::printf("refined on %zu reflections (%zu rejected as outliers, %zu "
                "with an undetermined rotation angle), %d steps\n",
                result.n_used, result.n_rejected, result.n_ill_conditioned,
                result.iterations);
    // Said plainly, because a residual averaged over a different set of
    // reflections is not comparable with anything -- including dials.refine,
    // which applies the same cutoff and reports over what is left.
    std::printf("rmsd is over those %zu reflections, not all %zu\n",
                result.n_used, reflections.nrows);
    std::printf("rmsd %.4f px  %.4f px  %.4f images\n", result.rmsd_x,
                result.rmsd_y, result.rmsd_z);
    for (std::size_t i = 0; i < experiments.size(); ++i) {
      if (!experiments[i].crystal) continue;
      const UnitCell c = experiments[i].crystal->cell();
      std::printf("  [%zu] cell %.4f %.4f %.4f  %.3f %.3f %.3f  V %.1f\n", i, c.a,
                  c.b, c.c, c.alpha, c.beta, c.gamma, c.volume());
    }

    // s1, rlp and entering follow the refined model, as they do in DIALS.
    // xyzobs.mm deliberately does not: it is what the spot finder measured
    // through the model as imported, and recomputing it here would silently
    // change the observations refinement was just fitted to.
    set_refinement_flags(result, reflections);
    add_reciprocal_columns(experiments, reflections);
    update_predictions(experiments, reflections);
    const double t_write_start = now_wall();
    write_experiments(out_expt, experiments);
    write_reflections(out_refl, reflections);
    t_write = now_wall() - t_write_start;
    std::printf("wrote %s and %s\n", out_expt.c_str(), out_refl.c_str());

    if (args.has("--timing")) {
      const double total = now_wall() - t_start;
      const auto line = [&](const char *name, double seconds) {
        std::printf("  %-24s %7.3f s  %5.1f%%\n", name, seconds,
                    total > 0.0 ? 100.0 * seconds / total : 0.0);
      };
      std::printf("\ntiming\n");
      line("read", t_read);
      line("build the target rows", g_observations_seconds);
      line("the jacobian", g_jacobian_seconds);
      line("the normal equations", g_normal_seconds);
      line("the solve", g_solve_seconds);
      line("the trial residuals", g_residual_seconds);
      line("outlier rejection", g_outlier_seconds);
      line("write", t_write);
      std::printf("  %-24s %7.3f s\n", "total", total);
    }
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_refine: %s\n", e.what());
    return 1;
  }
}

}  // namespace mxi

int main(int argc, char **argv) { return mxi::run_program(argc, argv); }
