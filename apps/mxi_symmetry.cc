// mxi_symmetry: the Laue group and space group of one sweep's integrated
// data, and the data reindexed into it -- as dials.symmetry does. The lattice's
// metric symmetry by Le Page's method; each element and each subgroup scored
// as Evans (2011); screw axes from the absences.

#include <cmath>
#include <cstdio>
#include <set>
#include <string>

#include "args.hh"
#include "expt.hh"
#include "laue.hh"
#include "log_mirror.hh"
#include "refl.hh"
#include "symmetry.hh"

namespace mxi {

namespace {

void usage() {
  std::printf(
      "usage: mxi_symmetry [options] INTEGRATED_EXPT INTEGRATED_REFL\n"
      "  --max-delta D     the lattice's symmetry to D degrees of obliquity "
      "(2)\n"
      "  --output-expt PATH   the models, reindexed (symmetrized.expt)\n"
      "  --output-refl PATH   the reflections, reindexed (symmetrized.refl)\n");
}

} // namespace

int run_program(int argc, char **argv) {
  const std::set<std::string> known = {"--max-delta", "--output-expt",
                                       "--output-refl"};
  const Arguments args = parse_arguments(argc, argv, known, known);
  if (args.help) {
    usage();
    return 0;
  }
  if (!args.ok) {
    std::fprintf(stderr, "mxi_symmetry: %s\n", args.error.c_str());
    return 2;
  }
  if (args.positional.size() != 2) {
    std::fprintf(stderr,
                 "mxi_symmetry: expected an .expt and a .refl, got %zu files\n",
                 args.positional.size());
    usage();
    return 2;
  }
  try {
    ExperimentList experiments = read_experiments(args.positional[0]);
    Table reflections = read_reflections(args.positional[1]);
    if (experiments.size() != 1 || !experiments[0].crystal)
      throw std::runtime_error("one sweep with a crystal, for now");
    const UnitCell cell = experiments[0].crystal->cell();
    std::printf("Cell %.3f %.3f %.3f A, %.3f %.3f %.3f deg\n", cell.a, cell.b,
                cell.c, cell.alpha, cell.beta, cell.gamma);

    const P1Intensities merged = merge_in_p1(experiments, reflections);
    P1Intensities normalised = merged;
    normalise(normalised);
    const std::vector<Rotation> lattice =
        lattice_symmetry(cell, args.number("--max-delta", 2.0));
    const LaueScores scores = score_laue_groups(normalised, lattice);
    std::printf("%zu reflections merged in P1; the lattice has %zu rotations, "
                "%zu symmetry "
                "elements, %zu subgroups\n",
                merged.size(), lattice.size(), scores.elements.size(),
                scores.groups.size());
    std::printf("E(CC; S) %.3f, sigma(CC) factor %.3f\n", scores.cc_true,
                scores.cc_sig_fac);

    std::printf("\nScoring each symmetry element\n");
    std::printf("  %-4s %5s %7s %6s %7s %10s\n", "", "order", "pairs", "CC",
                "Z-CC", "likelihood");
    for (std::size_t k = 0; k < scores.elements.size(); ++k) {
      const ElementScore &e = scores.elements[k];
      std::printf("  %-4c %5d %7zu %6.3f %7.2f %10.3f\n",
                  static_cast<char>('a' + k), e.order, e.pairs, e.cc, e.z,
                  e.likelihood);
    }

    std::printf("\nScoring all possible subgroups\n");
    std::printf("  %-14s %3s %10s %7s %6s %6s %5s %5s  %-20s %s\n",
                "Patterson group", "", "likelihood", "NetZcc", "Zcc+", "Zcc-",
                "CC", "CC-", "reindex", "elements");
    std::optional<Setting> best;
    for (std::size_t k = 0; k < scores.groups.size(); ++k) {
      const GroupScore &g = scores.groups[k];
      const std::optional<Setting> s = reference_setting(g.rotations, cell);
      if (k == 0)
        best = s;
      std::string letters;
      for (std::size_t j = 0; j < g.contains.size(); ++j)
        letters += g.contains[j] ? static_cast<char>('a' + j) : '.';
      std::printf(
          "  %-14s %3s %10.3f %7.2f %6.2f %6.2f %5.2f %5.2f  %-20s %s\n",
          s ? s->group.name().c_str() : "?", k == 0 ? "***" : "", g.likelihood,
          g.z_net, g.z_for, g.z_against, g.cc_for, g.cc_against,
          s ? s->cb_text.c_str() : "?", letters.c_str());
    }
    if (!best)
      throw std::runtime_error(
          "the most likely group could not be put in a reference setting");

    // Screw axes from the absences, in the chosen group's own setting.
    std::vector<Miller> hkl;
    for (const Miller &h : merged.hkl)
      hkl.push_back(best->cb.apply(h));
    const SpaceGroupChoice choice =
        choose_space_group(hkl, merged.i, merged.sigma, best->group);
    std::printf(
        "\nSpace groups with Patterson group %s, judged by their absences\n",
        best->group.name().c_str());
    std::printf("  %-14s %8s %10s  %s\n", "space group", "tested", "<I/sigma>",
                "");
    for (const AbsenceTest &t : choice.candidates)
      std::printf("  %-14s %8zu %10.2f  %s\n", t.group.name().c_str(), t.tested,
                  t.mean_i_over_sigma,
                  t.consistent ? "consistent" : "inconsistent");
    std::printf("\nRecommended space group: %s\n",
                choice.chosen.name().c_str());
    for (const std::string &n : choice.indistinguishable)
      std::printf("  indistinguishable by its absences from %s\n", n.c_str());
    std::printf("Reindex operator: %s\n", best->cb_text.c_str());

    reindex(experiments, reflections, best->cb, choice.chosen);
    const UnitCell after = experiments[0].crystal->cell();
    std::printf("Cell %.3f %.3f %.3f A, %.3f %.3f %.3f deg\n", after.a, after.b,
                after.c, after.alpha, after.beta, after.gamma);
    const std::string out_expt =
        args.value("--output-expt", "symmetrized.expt");
    const std::string out_refl =
        args.value("--output-refl", "symmetrized.refl");
    write_experiments(out_expt, experiments);
    write_reflections(out_refl, reflections);
    std::printf("\nWrote %s and %s\n", out_expt.c_str(), out_refl.c_str());
  } catch (const std::exception &error) {
    std::fprintf(stderr, "mxi_symmetry: %s\n", error.what());
    return 1;
  }
  return 0;
}

} // namespace mxi

int main(int argc, char **argv) {
  // Mirrored to mxi_symmetry.log in the working directory, as DIALS writes
  // dials.symmetry.log; not for a run that only asks for help.
  if (!mxi::only_asks_for_help(argc, argv))
    mxi::mirror_to_log("mxi_symmetry.log");
  return mxi::run_program(argc, argv);
}
