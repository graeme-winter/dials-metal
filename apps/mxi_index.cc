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

#include "expt.h"
#include "index.h"
#include "refl.h"

using namespace mxi;

namespace {

void usage() {
  std::printf(
      "usage: mxi_index IMPORTED.expt STRONG.refl [options]\n"
      "  --d-min D        resolution limit (default: from the data)\n"
      "  --max-cell A     longest cell edge (default: from spot spacing)\n"
      "  --grid N         FFT grid size, power of two (default: from d_min)\n"
      "  --tolerance T    how far an index may fall from an integer (0.3)\n"
      "  --candidates N   basis vectors taken from the peak list (30)\n"
      "  --output-expt P  (default indexed.expt)\n"
      "  --output-refl P  (default indexed.refl)\n"
      "  --quiet\n");
}

double number(const char *text, const char *what) {
  char *end = nullptr;
  const double v = std::strtod(text, &end);
  if (end == text) {
    std::fprintf(stderr, "mxi_index: %s needs a number, got '%s'\n", what, text);
    std::exit(2);
  }
  return v;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    usage();
    return 2;
  }
  const std::string expt_path = argv[1];
  const std::string refl_path = argv[2];
  std::string out_expt = "indexed.expt";
  std::string out_refl = "indexed.refl";
  IndexOptions options;
  options.verbose = true;

  for (int i = 3; i < argc; ++i) {
    const std::string arg = argv[i];
    const auto next = [&](const char *what) -> const char * {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "mxi_index: %s needs a value\n", what);
        std::exit(2);
      }
      return argv[++i];
    };
    if (arg == "--d-min") options.d_min = number(next("--d-min"), "--d-min");
    else if (arg == "--max-cell") options.max_cell = number(next("--max-cell"), "--max-cell");
    else if (arg == "--grid") options.grid = static_cast<std::size_t>(number(next("--grid"), "--grid"));
    else if (arg == "--tolerance") options.tolerance = number(next("--tolerance"), "--tolerance");
    else if (arg == "--candidates") options.n_candidates = static_cast<std::size_t>(number(next("--candidates"), "--candidates"));
    else if (arg == "--output-expt") out_expt = next("--output-expt");
    else if (arg == "--output-refl") out_refl = next("--output-refl");
    else if (arg == "--quiet") options.verbose = false;
    else if (arg == "-h" || arg == "--help") { usage(); return 0; }
    else {
      std::fprintf(stderr, "mxi_index: unknown option '%s'\n", arg.c_str());
      return 2;
    }
  }

  try {
    ExperimentList experiments = read_experiments(expt_path);
    Table reflections = read_reflections(refl_path);
    if (options.verbose) {
      std::printf("%zu experiments, %zu reflections\n", experiments.size(),
                  reflections.nrows);
      for (const std::string &d : reflections.dropped()) {
        std::printf("  not read: %s\n", d.c_str());
      }
    }

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

    write_experiments(out_expt, experiments);
    write_reflections(out_refl, reflections);
    std::printf("wrote %s and %s\n", out_expt.c_str(), out_refl.c_str());
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_index: %s\n", e.what());
    return 1;
  }
}
