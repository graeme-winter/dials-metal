// mxi_index: index a strong spot list against an imported experiment.
//
//   mxi_index imported.expt strong.refl [options]
//
// Writes indexed.expt and indexed.refl, so the result can be compared against
// DIALS' own with `mxeq check indexed`.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include <set>

#include "args.h"
#include "expt.h"
#include "index.h"
#include "refine.h"
#include "refl.h"

using namespace mxi;

namespace {

void usage() {
  std::printf(
      "usage: mxi_index IMPORTED.expt STRONG.refl [options]\n"
      "  --d-min D        resolution limit (default: from the data)\n"
      "  --max-cell A     longest cell edge (default: from spot spacing)\n"
      "  --grid N         FFT grid size, power of two (default: from d_min)\n"      "  --timing         where the time went, by phase\n"
      "  --tolerance T    how far an index may fall from an integer (0.3)\n"
      "  --candidates N   basis vectors taken from the peak list (30)\n"
      "  --output-expt P  (default indexed.expt)\n"
      "  --output-refl P  (default indexed.refl)\n"
      "  --macrocycles N  assign/refine/re-assign cycles (3)\n"
      "  --all-reflections  refine on everything, not the stronger half\n"
      "  --quiet\n");
}

}  // namespace

int main(int argc, char **argv) {
  const std::set<std::string> known = {
      "--d-min",       "--max-cell",    "--grid",          "--tolerance",
      "--candidates",  "--output-expt", "--output-refl",   "--quiet",
      "--macrocycles", "--all-reflections", "--timing"};
  const std::set<std::string> takes_value = {
      "--d-min",      "--max-cell",    "--grid",        "--tolerance",
      "--candidates", "--output-expt", "--output-refl", "--macrocycles"};
  const Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (args.help) {
    usage();
    return 0;
  }
  if (!args.ok) {
    std::fprintf(stderr, "mxi_index: %s\n", args.error.c_str());
    return 2;
  }
  if (args.positional.size() != 2) {
    std::fprintf(stderr,
                 "mxi_index: expected an .expt and a .refl, got %zu file "
                 "arguments\n",
                 args.positional.size());
    usage();
    return 2;
  }

  IndexOptions options;
  options.verbose = !args.has("--quiet");
  options.d_min = args.number("--d-min", 0.0);
  options.max_cell = args.number("--max-cell", 0.0);
  options.grid = static_cast<std::size_t>(args.number("--grid", 0));
  options.tolerance = args.number("--tolerance", 0.3);
  options.n_candidates = static_cast<std::size_t>(args.number("--candidates", 30));
  options.macrocycles = static_cast<int>(args.number("--macrocycles", 3));
  options.refine_on_strong = !args.has("--all-reflections");
  const std::string out_expt = args.value("--output-expt", "indexed.expt");
  const std::string out_refl = args.value("--output-refl", "indexed.refl");

  try {
    ExperimentList experiments = read_experiments(args.positional[0]);
    Table reflections = read_reflections(args.positional[1]);
    if (options.verbose) {
      std::printf("%zu experiments, %zu reflections\n", experiments.size(),
                  reflections.nrows);
      for (const std::string &d : reflections.dropped()) {
        std::printf("  not read: %s\n", d.c_str());
      }
    }

    // Before indexing, because these are observations expressed through the
    // detector model and dials.find_spots writes them from the model as
    // imported. Computing them from the model indexing has just refined would
    // be a difference from DIALS dressed up as a correction -- and DIALS never
    // recomputes them either, which is exactly why xyzobs.mm.value must not be
    // used as a join key between two runs.
    add_observed_columns(experiments, reflections);

    const IndexResult result = index(experiments, reflections, options);
    if (result.n_indexed == 0) {
      std::fprintf(stderr,
                   "mxi_index: no lattice found. Try --max-cell, or --d-min to "
                   "restrict the resolution range used.\n");
      return 1;
    }

    const UnitCell cell = result.crystal.cell();
    std::printf("indexed %zu of %zu (%.1f%%)  rmsd %.4f\n", result.n_indexed,
                result.n_total, 100.0 * result.fraction_indexed(),
                result.rmsd_index);
    std::printf("cell %.4f %.4f %.4f  %.3f %.3f %.3f   volume %.1f\n", cell.a,
                cell.b, cell.c, cell.alpha, cell.beta, cell.gamma, cell.volume());

    set_indexed_flags(reflections);
    add_reciprocal_columns(experiments, reflections);
    update_predictions(experiments, reflections);
    write_experiments(out_expt, experiments);
    write_reflections(out_refl, reflections);
    std::printf("wrote %s and %s\n", out_expt.c_str(), out_refl.c_str());

    if (args.has("--timing")) {
      const IndexTiming &t = result.timing;
      const auto line = [&](const char *name, double seconds) {
        std::printf("  %-22s %7.3f s  %5.1f%%\n", name, seconds,
                    t.total > 0.0 ? 100.0 * seconds / t.total : 0.0);
      };
      std::printf("\ntiming\n");
      line("reciprocal points", t.reciprocal_points);
      line("max cell", t.max_cell);
      line("candidate vectors", t.candidate_vectors);
      line("  the transform", t.fft);
      line("  the peak search", t.peak_search);
      line("  the rest of it",
           t.candidate_vectors - t.fft - t.peak_search);
      line("choose basis", t.choose_basis);
      std::printf("    %zu triples scored, %zu skipped as degenerate\n",
                  t.triples_scored, t.triples_skipped);
      line("fit and reduce", t.fit_and_reduce);
      line("macrocycles", t.macrocycles);
      line("  copy and select", t.subset_copy);
      line("  refinement", t.refine);
      line("  reassignment", t.reassign);
      line("    the jacobian", g_jacobian_seconds);
      line("    the normal equations", g_normal_seconds);
      std::printf("  %-22s %7.3f s\n", "indexing total", t.total);
    }
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_index: %s\n", e.what());
    return 1;
  }
}
