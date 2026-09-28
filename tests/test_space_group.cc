#include <cmath>
#include <random>

#include "../src/laue.hh"
#include "check.hh"

namespace mxi {

TEST(the_candidates_for_a_laue_group_are_its_chiral_space_groups) {
  // For m-3 with I centring: I 2 3 and I 21 3, not the centrosymmetric I m -3.
  const std::vector<SpaceGroup> i =
      space_groups_with_patterson(SpaceGroup::from_name("I m -3"));
  check::equal(static_cast<long long>(i.size()), 2, "two");
  check::is_true(i[0].name() == "I 2 3" && i[1].name() == "I 21 3",
                 i[0].name() + ", " + i[1].name());
  // For mmm, primitive: the four Sohncke groups of 222.
  const std::vector<SpaceGroup> p =
      space_groups_with_patterson(SpaceGroup::from_name("P m m m"));
  check::equal(static_cast<long long>(p.size()), 4,
               "P 2 2 2, P 2 2 21, P 21 21 2, P 21 21 21");
}

namespace {

//: Merged intensities in an orthorhombic cell whose true group is `truth`:
//: strong where it allows a reflection, noise of unit sigma where it forbids.
void planted(const char *truth, std::vector<Miller> *hkl,
             std::vector<double> *i, std::vector<double> *s, unsigned seed) {
  const SpaceGroup g = SpaceGroup::from_name(truth);
  std::mt19937 rng(seed);
  std::exponential_distribution<double> wilson(1.0 / 400.0);
  std::normal_distribution<double> noise(0.0, 1.0);
  for (int h = 0; h <= 10; ++h)
    for (int k = 0; k <= 10; ++k)
      for (int l = 0; l <= 10; ++l) {
        if (h == 0 && k == 0 && l == 0)
          continue;
        const Miller m{h, k, l};
        hkl->push_back(m);
        const double truth_i = g.absent(m) ? 0.0 : 20.0 + wilson(rng);
        const double sigma = std::sqrt(truth_i + 25.0);
        i->push_back(truth_i + sigma * noise(rng));
        s->push_back(sigma);
      }
}

} // namespace

TEST(the_space_group_planted_is_the_space_group_found) {
  const SpaceGroup mmm = SpaceGroup::from_name("P m m m");
  for (const char *truth : {"P 2 2 2", "P 2 2 21", "P 21 21 2", "P 21 21 21"}) {
    std::vector<Miller> hkl;
    std::vector<double> i, s;
    planted(truth, &hkl, &i, &s, 5);
    const SpaceGroupChoice c = choose_space_group(hkl, i, s, mmm);
    check::is_true(c.chosen.name() == truth, std::string("planted ") + truth +
                                                 ", found " + c.chosen.name());
  }
}

TEST(what_centring_hides_is_reported_as_indistinguishable) {
  // Under I centring the 21 axes of I 21 3 forbid nothing more: both
  // candidates are consistent with anything, and the choice is the lower.
  std::vector<Miller> hkl;
  std::vector<double> i, s;
  const SpaceGroup i23 = SpaceGroup::from_name("I 2 3");
  for (int h = 0; h <= 8; ++h)
    for (int k = 0; k <= 8; ++k)
      for (int l = 0; l <= 8; ++l) {
        const Miller m{h, k, l};
        if ((h || k || l) && !i23.absent(m)) {
          hkl.push_back(m);
          i.push_back(100.0);
          s.push_back(10.0);
        }
      }
  const SpaceGroupChoice c =
      choose_space_group(hkl, i, s, SpaceGroup::from_name("I m -3"));
  check::is_true(c.chosen.name() == "I 2 3",
                 "I 2 3, the lower number: " + c.chosen.name());
  check::is_true(c.indistinguishable.size() == 1 &&
                     c.indistinguishable[0] == "I 21 3",
                 "with I 21 3 reported as indistinguishable");
}

} // namespace mxi
