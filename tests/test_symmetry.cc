#include <cmath>
#include <set>

#include "../src/symmetry.hh"
#include "check.hh"

namespace mxi {

TEST(every_equivalent_of_a_reflection_in_m3_is_one_reflection) {
  // I 2 3, Laue class m-3: 24 operations counting Friedel mates -- the three
  // cyclic permutations of (h, k, l), times sign changes in pairs, times
  // inversion. All 24 images of (1, 2, 3) are one reflection.
  const SpaceGroup g = SpaceGroup::from_name("I 2 3");
  check::is_true(g.laue() == "m-3", "the Laue class is m-3: " + g.laue());
  check::equal(static_cast<long long>(g.order()), 12,
               "12 operations without centring");
  const Miller h{1, 2, 3};
  const Miller u = g.unique(h);
  std::set<Miller> images;
  const int signs[4][3] = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
  for (int c = 0; c < 3; ++c) {
    const Miller p{h[c % 3], h[(c + 1) % 3], h[(c + 2) % 3]};
    for (const auto &s : signs) {
      for (int inv : {1, -1}) {
        const Miller e{inv * s[0] * p[0], inv * s[1] * p[1], inv * s[2] * p[2]};
        images.insert(e);
        check::is_true(g.unique(e) == u,
                       "every image of (1,2,3) is one reflection");
      }
    }
  }
  check::equal(static_cast<long long>(images.size()), 24,
               "and there are 24 of them");
  // A non-cyclic permutation is m-3m's symmetry, not m-3's: a different
  // reflection. This is what would catch the wrong Laue class.
  check::is_true(g.unique({1, 3, 2}) != u,
                 "(1,3,2) is another reflection in m-3");
  const SpaceGroup m3m = SpaceGroup::from_name("I 4 3 2");
  check::is_true(m3m.unique({1, 3, 2}) == m3m.unique(h),
                 "but the same one in m-3m");
}

TEST(centring_forbids_reflections_and_p1_pairs_only_friedel_mates) {
  const SpaceGroup g = SpaceGroup::from_name("I 2 3");
  check::is_true(g.absent({1, 0, 0}),
                 "h + k + l odd is absent under I centring");
  check::is_true(!g.absent({1, 1, 0}), "and even is not");
  const SpaceGroup p1 = SpaceGroup::from_hall(" P 1");
  check::equal(static_cast<long long>(p1.number()), 1,
               "an .expt's ' P 1' is P 1");
  check::is_true(p1.unique({1, 2, 3}) == p1.unique({-1, -2, -3}),
                 "Friedel mates merge");
  check::is_true(p1.unique({1, 2, 3}) != p1.unique({1, 2, -3}),
                 "and nothing else does");
  check::is_true(SpaceGroup::from_hall(g.hall()).name() == g.name(),
                 "a Hall symbol written out reads back as the same group");
}

TEST(a_change_of_basis_takes_the_primitive_cell_to_the_conventional) {
  // DIALS' b+c,a+c,a+b: the primitive cell of a body-centred cubic lattice to
  // the conventional one. hkl' = M^T hkl, so (1,2,3) -> (2+3, 1+3, 1+2).
  const ChangeOfBasis cb = ChangeOfBasis::parse("b+c,a+c,a+b");
  check::is_true(cb.apply({1, 2, 3}) == Miller{5, 4, 3}, "(1,2,3) -> (5,4,3)");
  const Miller n = cb.apply({1, 0, 0});
  check::equal(static_cast<long long>((n[0] + n[1] + n[2]) % 2), 0,
               "every reindexed index obeys I centring");
  // Its inverse, with fractions, undoes it; and a fractional index is refused.
  const ChangeOfBasis back =
      ChangeOfBasis::parse("-1/2a+1/2b+1/2c,1/2a-1/2b+1/2c,1/2a+1/2b-1/2c");
  check::is_true(back.apply({5, 4, 3}) == Miller{1, 2, 3},
                 "the inverse undoes it");
  bool refused = false;
  try {
    back.apply({1, 0, 0});
  } catch (const std::exception &) {
    refused = true;
  }
  check::is_true(refused,
                 "an index with no integral image is refused, not rounded");

  // A crystal: a cubic cell of 78 A, taken to its primitive cell and back.
  Crystal c;
  c.A = Mat3::identity() * (1.0 / 78.0);
  back.apply(c);
  const UnitCell primitive = c.cell();
  check::close(primitive.a, 78.0 * std::sqrt(3.0) / 2.0, 1e-9,
               "primitive a = a sqrt(3)/2");
  check::close(primitive.alpha, std::acos(-1.0 / 3.0) * 180.0 / std::acos(-1.0),
               1e-9, "and 109.47 degrees");
  cb.apply(c);
  const UnitCell cubic = c.cell();
  check::close(cubic.a, 78.0, 1e-9, "and back to 78 A");
  check::close(cubic.beta, 90.0, 1e-9, "cubic");
}

TEST(a_change_of_basis_is_refused_when_it_cannot_be_one) {
  for (const char *bad : {"a,b", "a+b,a+b,c", "a,b,d", "a b,b,c"}) {
    bool refused = false;
    try {
      ChangeOfBasis::parse(bad);
    } catch (const std::exception &) {
      refused = true;
    }
    check::is_true(refused, std::string("refused: ") + bad);
  }
}

} // namespace mxi
