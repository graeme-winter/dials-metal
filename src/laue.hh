#pragma once

// Choosing the Laue group: each symmetry element of the lattice scored by the
// correlation of the reflections it relates, and each subgroup by the product
// of its elements' likelihoods -- Evans (2011), Acta Cryst. D67, 282-292,
// appendices A1 and A2, as dials.symmetry implements them.

#include <cstddef>
#include <string>
#include <vector>

#include "symmetry.hh"

namespace mxi {

//: Intensities merged in P1 with Friedel mates merged, one per reflection.
struct P1Intensities {
  std::vector<Miller> hkl; //: the P -1 unique index
  std::vector<double> i, sigma, d;
  std::size_t size() const { return hkl.size(); }
};

//: Divided by <I> in shells of equal count by resolution -- the
//: quasi-normalisation E^2, where dials.symmetry fits an anisotropic
//: maximum-likelihood model -- so that correlations compare like with like.
void normalise(P1Intensities &data, std::size_t per_shell = 200);

struct ElementScore {
  Rotation rotation;
  int order = 1;
  std::size_t pairs = 0;
  double cc = 0.0, sigma_cc = 0.0, z = 0.0;
  double p_given_present = 0.0, p_given_absent = 0.0, likelihood = 0.0;
  //: For pooling into a group's CC: sums of x, y, xx, yy, xy over the pairs.
  double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
};

struct GroupScore {
  std::vector<Rotation> rotations;
  std::vector<bool> contains; //: per element
  double likelihood = 0.0, z_for = 0.0, z_against = 0.0, z_net = 0.0;
  double cc_for = 0.0, cc_against = 0.0;
};

struct LaueScores {
  double cc_true = 0.0, cc_sig_fac = 0.0, e_cc_true = 0.0, cc_identity = 0.0;
  std::vector<ElementScore> elements;
  std::vector<GroupScore> groups; //: most likely first
};

//: The distinct symmetry elements of a lattice's rotations, as dials.symmetry
//: counts them: each rotation axis once, a rotation and its inverse together,
//: and the square of a four-fold as a two-fold of its own. 17 for m-3m.
std::vector<Rotation> symmetry_elements(const std::vector<Rotation> &lattice);
int rotation_order(const Rotation &r);

//: Score every element and every subgroup of the lattice's rotations.
LaueScores score_laue_groups(const P1Intensities &normalised,
                             const std::vector<Rotation> &lattice,
                             unsigned seed = 1);

//: Evans (2011) A1: p(CC; S), a Cauchy centred on E(CC; S), and p(CC; !S),
//: the same averaged over a true CC with density (1 - x^2)^(1/2) on [0, 1],
//: both truncated to [-1, 1].
double p_cc_given_present(double cc, double sigma_cc, double expected);
double p_cc_given_absent(double cc, double sigma_cc);

} // namespace mxi
