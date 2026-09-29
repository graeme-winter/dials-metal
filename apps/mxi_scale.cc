// mxi_scale: put every observation of one sweep on one scale, refine an error
// model, and report merging statistics. The physical model of Beilsten-Edmands
// et al. (2020), with B-splines for the scale and decay; see
// docs/scaling.md.

#include <cmath>
#include <cstdio>
#include <set>
#include <string>

#include "args.hh"
#include "expt.hh"
#include "log_mirror.hh"
#include "parallel.hh"
#include "refl.hh"
#include "resolution.hh"
#include "scale.hh"
#include "symmetry.hh"
#include "timing.hh"

namespace mxi {

namespace {

void usage() {
  std::printf(
      "usage: mxi_scale [options] INTEGRATED_EXPT INTEGRATED_REFL\n"
      "  --space-group NAME  the space group, \"I 2 3\" or \"197\"; default "
      "the\n"
      "                    crystal's own\n"
      "  --change-of-basis CB   to that group's setting first, in DIALS'\n"
      "                    notation: \"b+c,a+c,a+b\" from the primitive cell "
      "of\n"
      "                    a body-centred cubic lattice\n"
      "  --d-min D         leave out reflections beyond D A\n"
      "  --no-absorption   no absorption surface, whatever the sweep\n"
      "  --profile-only    profile-fitted intensities alone, not a mix with\n"
      "                    summation chosen by Rmeas\n"
      "  --shells N        resolution shells in the table (20)\n"
      "  -o PATH           scaled reflections (scaled.refl)\n"
      "  --output-expt PATH   the models, reindexed (scaled.expt)\n"
      "  --threads N       threads to use (0: the machine's); the answer is "
      "the\n"
      "                    same for any number\n"
      "  --timing          where the time goes\n");
}

} // namespace

int run_program(int argc, char **argv) {
  const std::set<std::string> known = {
      "--threads", "--timing",        "--space-group",  "--change-of-basis",
      "--d-min",   "--no-absorption", "--profile-only", "--shells",
      "-o",        "--output-expt"};
  const std::set<std::string> takes_value = {
      "--threads", "--space-group", "--change-of-basis", "--d-min", "--shells",
      "-o",        "--output-expt"};
  const Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (args.help) {
    usage();
    return 0;
  }
  if (!args.ok) {
    std::fprintf(stderr, "mxi_scale: %s\n", args.error.c_str());
    return 2;
  }
  if (args.positional.size() != 2) {
    std::fprintf(stderr,
                 "mxi_scale: expected an .expt and a .refl, got %zu files\n",
                 args.positional.size());
    usage();
    return 2;
  }
  Timing timing(args.has("--timing"));
  set_parallel_threads(static_cast<std::size_t>(args.number("--threads", 0.0)));
  try {
    double mark = Timing::now();
    const auto phase = [&](const char *name) {
      const double t = Timing::now();
      timing.add(name, t - mark);
      mark = t;
    };
    ExperimentList experiments = read_experiments(args.positional[0]);
    Table reflections = read_reflections(args.positional[1]);
    phase("reading");
    if (experiments.size() != 1 || !experiments[0].crystal)
      throw std::runtime_error("one sweep with a crystal, for now");
    const SpaceGroup group =
        args.has("--space-group")
            ? SpaceGroup::from_name(args.value("--space-group", ""))
            : SpaceGroup::from_hall(experiments[0].crystal->space_group_hall);
    if (args.has("--change-of-basis")) {
      const std::string cb = args.value("--change-of-basis", "");
      reindex(experiments, reflections, ChangeOfBasis::parse(cb), group);
      std::printf("Reindexed by %s\n", cb.c_str());
    } else {
      experiments[0].crystal->space_group_hall = group.hall();
    }
    const UnitCell cell = experiments[0].crystal->cell();
    std::printf("Space group %s, Laue class %s; cell %.3f %.3f %.3f A, %.3f "
                "%.3f %.3f deg\n",
                group.name().c_str(), group.laue().c_str(), cell.a, cell.b,
                cell.c, cell.alpha, cell.beta, cell.gamma);

    ScaleRunOptions options;
    options.combine = !args.has("--profile-only");
    options.absorption = !args.has("--no-absorption");
    options.d_min = args.number("--d-min", 0.0);
    phase("reindexing");
    const ScaleRun run = scale_sweep(experiments, reflections, group, options);
    phase("scaling");
    for (const auto &[name, seconds] : run.timing)
      timing.add(name, seconds, 1);
    const ScaleData &data = run.data;
    if (data.size() == 0)
      throw std::runtime_error("no observations fit to scale");

    const ScaleModelShape &shape = run.model.shape();
    std::printf("\n%zu observations of %zu reflections, of %zu rows\n",
                data.size(), data.unique.size(), reflections.nrows);
    std::printf("Scaling model: %zu scale, %zu decay and %zu absorption "
                "parameters, %zu in "
                "all; fitted on %zu observations\n",
                shape.scale_points, shape.decay_points,
                harmonic_count(shape.lmax), run.model.size(), run.fitted_on);
    for (std::size_t k = 0; k < run.fits.size(); ++k)
      std::printf("  fit %zu: %d steps, target %.6g to %.6g%s\n", k + 1,
                  run.fits[k].iterations, run.fits[k].target_start,
                  run.fits[k].target_end,
                  run.fits[k].converged ? "" : ", not converged");
    double cmin = HUGE_VAL, cmax = 0.0, bmin = HUGE_VAL, bmax = -HUGE_VAL;
    for (std::size_t i = 0; i < shape.scale_points; ++i) {
      cmin = std::fmin(cmin, run.model.parameters[i]);
      cmax = std::fmax(cmax, run.model.parameters[i]);
    }
    for (std::size_t i = 0; i < shape.decay_points; ++i) {
      bmin = std::fmin(bmin, run.model.parameters[run.model.first_decay() + i]);
      bmax = std::fmax(bmax, run.model.parameters[run.model.first_decay() + i]);
    }
    std::printf("  scale %.4f to %.4f", cmin, cmax);
    if (shape.decay_points > 0)
      std::printf("; relative B %.3f to %.3f A^2", bmin, bmax);
    std::printf("\n");
    if (run.i_mid == 0.0)
      std::printf("Intensities: profile fitted\n");
    else if (std::isinf(run.i_mid))
      std::printf("Intensities: summation\n");
    else
      std::printf(
          "Intensities: profile fitted, crossing to summation at I = %.0f\n",
          run.i_mid);
    if (run.covariance.ok) {
      double big = 0.0, total = 0.0;
      for (std::size_t k = 0; k < run.g_variance.size(); ++k) {
        const double f = std::sqrt(run.g_variance[k]) / run.g[k];
        big = std::fmax(big, f);
        total += f;
      }
      std::printf("Scale uncertainties: sigma(g)/g %.4f on average, %.4f at "
                  "most; goodness of "
                  "fit %.3f on %zu degrees of freedom; carried into each "
                  "variance before the "
                  "error model\n",
                  total / static_cast<double>(run.g_variance.size()), big,
                  run.covariance.goodness_of_fit,
                  run.covariance.degrees_of_freedom);
    } else {
      std::printf("Scale uncertainties: the normal matrix could not be "
                  "inverted; none carried\n");
    }
    std::printf("Error model: a = %.4f, b = %.4f, from %zu observations\n",
                run.error_model.a, run.error_model.b, run.error_model.used);
    std::printf("%zu outliers\n", run.outliers);

    MergingShell all;
    const int shells = static_cast<int>(args.number("--shells", 20.0));
    const std::vector<MergingShell> table = merging_statistics(
        data, run.g, group, *experiments[0].crystal, shells, &all);
    phase("merging statistics");
    std::printf("\nMerging statistics, Friedel mates merged\n");
    std::printf("  %6s %6s %7s %6s %6s %6s %8s %7s %6s %7s %6s %6s\n", "d_max",
                "d_min", "#obs", "#uniq", "mult.", "%comp", "<I>", "<I/sI>",
                "r_mrg", "r_meas", "r_pim", "cc1/2");
    const auto row = [](const MergingShell &m) {
      std::printf("  %6.2f %6.2f %7zu %6zu %6.2f %6.2f %8.1f %7.1f %6.3f %7.3f "
                  "%6.3f %6.3f\n",
                  m.d_max, m.d_min, m.observations, m.unique, m.multiplicity,
                  100.0 * m.completeness, m.mean_i, m.i_over_sigma, m.rmerge,
                  m.rmeas, m.rpim, m.cc_half);
    };
    for (const MergingShell &m : table)
      row(m);
    row(all);

    // The resolution limit, as dials.estimate_resolution finds it: a tanh in
    // d*^2 through CC half in 50 bins of equal count, where it crosses 0.3;
    // and where CC half stops being significant at the 0.1 level.
    std::size_t wilson = 0;
    const std::vector<ResolutionBin> bins = cc_half_bins(
        data, run.g, *experiments[0].crystal, 50, 10, 0.1, &wilson);
    const ResolutionEstimate res = estimate_resolution(bins, 0.3);
    phase("the resolution limit");
    std::printf("\nRemoving %zu Wilson outliers with E^2 >= 16.0\n", wilson);
    if (res.d_min_cc_half > 0.0)
      std::printf("Resolution cc_half:       %.2f\n", res.d_min_cc_half);
    else
      std::printf("Resolution cc_half:       none found\n");
    if (res.d_min_significance > 0.0)
      std::printf("Resolution cc_half_significance_level:    %.2f\n",
                  res.d_min_significance);
    else
      std::printf("Resolution cc_half_significance_level:    none found\n");

    // "Suggested": the same summary, cut at the CC half limit.
    MergingShell suggested;
    const bool cut = res.d_min_cc_half > 0.0;
    if (cut) {
      ScaleData to_limit = data;
      for (std::size_t i = 0; i < to_limit.size(); ++i) {
        const Miller &m = to_limit.unique[to_limit.group[i]];
        if (experiments[0].crystal->d_spacing(m[0], m[1], m[2]) <
            res.d_min_cc_half)
          to_limit.outlier[i] = true;
      }
      merging_statistics(to_limit, run.g, group, *experiments[0].crystal,
                         shells, &suggested);
    }
    phase("the suggested cut");

    // The summary dials.scale ends with: overall, the lowest shell and the
    // highest, and the data cut at the suggested limit.
    {
      const MergingShell &low = table.front(), &high = table.back();
      std::printf("\n            -------------Summary of merging "
                  "statistics--------------\n\n");
      std::printf("%-44s %8s %7s %7s %9s\n", "", "Overall", "Low", "High",
                  "Suggested");
      const auto line = [&](const char *what, double MergingShell::*f,
                            const char *fmt, double scale) {
        char x[4][32];
        const MergingShell *m[4] = {&all, &low, &high, &suggested};
        for (int k = 0; k < 4; ++k)
          std::snprintf(x[k], sizeof x[k], fmt, scale * (m[k]->*f));
        std::printf("%-44s %8s %7s %7s %9s\n", what, x[0], x[1], x[2],
                    cut ? x[3] : "");
      };
      const auto count = [&](const char *what, std::size_t MergingShell::*f) {
        if (cut)
          std::printf("%-44s %8zu %7zu %7zu %9zu\n", what, all.*f, low.*f,
                      high.*f, suggested.*f);
        else
          std::printf("%-44s %8zu %7zu %7zu\n", what, all.*f, low.*f, high.*f);
      };
      const auto only = [&](const char *what, double v, const char *fmt) {
        char x[32];
        std::snprintf(x, sizeof x, fmt, v);
        std::printf("%-44s %8s\n", what, x);
      };
      line("High resolution limit", &MergingShell::d_min, "%.2f", 1.0);
      line("Low resolution limit", &MergingShell::d_max, "%.2f", 1.0);
      line("Completeness", &MergingShell::completeness, "%.1f", 100.0);
      line("Multiplicity", &MergingShell::multiplicity, "%.1f", 1.0);
      line("I/sigma", &MergingShell::i_over_sigma, "%.1f", 1.0);
      line("Rmerge(I)", &MergingShell::rmerge, "%.3f", 1.0);
      line("Rmerge(I+/-)", &MergingShell::rmerge_anom, "%.3f", 1.0);
      line("Rmeas(I)", &MergingShell::rmeas, "%.3f", 1.0);
      line("Rmeas(I+/-)", &MergingShell::rmeas_anom, "%.3f", 1.0);
      line("Rpim(I)", &MergingShell::rpim, "%.3f", 1.0);
      line("Rpim(I+/-)", &MergingShell::rpim_anom, "%.3f", 1.0);
      line("CC half", &MergingShell::cc_half, "%.3f", 1.0);
      line("Anomalous completeness", &MergingShell::anom_completeness, "%.1f",
           100.0);
      line("Anomalous multiplicity", &MergingShell::anom_multiplicity, "%.1f",
           1.0);
      line("Anomalous correlation", &MergingShell::cc_anom, "%.3f", 1.0);
      only("Anomalous slope", all.anom_slope, "%.3f");
      only("dF/F", all.df_over_f, "%.3f");
      only("dI/s(dI)", all.di_over_sig_di, "%.3f");
      count("Total observations", &MergingShell::observations);
      count("Total unique", &MergingShell::unique);
    }

    // dials.scale's columns and flags, so that dials.merge and dials.export
    // take the table on.
    write_scaling(reflections, data, run.g, run.g_variance);
    const std::string out_refl = args.value("-o", "scaled.refl");
    const std::string out_expt = args.value("--output-expt", "scaled.expt");
    write_reflections(out_refl, reflections);
    write_experiments(out_expt, experiments);
    std::printf("\nWrote %s and %s\n", out_refl.c_str(), out_expt.c_str());
    phase("writing");
    timing.report(stdout);
  } catch (const std::exception &error) {
    std::fprintf(stderr, "mxi_scale: %s\n", error.what());
    return 1;
  }
  return 0;
}

} // namespace mxi

int main(int argc, char **argv) {
  // Mirrored to mxi_scale.log in the working directory, as DIALS writes
  // dials.scale.log; not for a run that only asks for help.
  if (!mxi::only_asks_for_help(argc, argv))
    mxi::mirror_to_log("mxi_scale.log");
  return mxi::run_program(argc, argv);
}
