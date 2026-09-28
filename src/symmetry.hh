#pragma once

// Space groups and changes of basis, for scaling and what follows it.
//
// gemmi does the group theory -- operators, the reciprocal asymmetric unit,
// absences -- behind this header, so that it is included in one file here
// rather than everywhere a Miller index is grouped.

#include <array>
#include <memory>
#include <string>

#include "expt.hh"
#include "refl.hh"

namespace mxi {

using Miller = std::array<int, 3>;

class SpaceGroup {
public:
  //: From a name, "I 2 3" or "P 1 21 1", or a number, "197".
  static SpaceGroup from_name(const std::string &name);
  //: From a Hall symbol, as an .expt holds it: " P 1", " I 2 2 3".
  static SpaceGroup from_hall(const std::string &hall);

  std::string name() const; //: Hermann-Mauguin, "I 2 3"
  std::string hall() const; //: as an .expt holds it
  std::string laue() const; //: the Laue class, "m-3"
  int number() const;
  //: Operations of the point group, without centring: 12 for I 2 3.
  std::size_t order() const;

  //: The symmetry-unique index of hkl, Friedel mates merged: the index that
  //: groups every observation of one reflection for scaling.
  Miller unique(const Miller &hkl) const;
  //: Whether hkl is the I+ of its Friedel pair -- an odd MTZ ISYM -- or
  //: centric, where the two are one reflection.
  bool friedel_plus(const Miller &hkl) const;
  //: Whether hkl is forbidden by centring or a screw axis or glide.
  bool absent(const Miller &hkl) const;
  bool centric(const Miller &hkl) const;

private:
  struct Impl;
  std::shared_ptr<const Impl> impl_;
};

//: A change of basis in DIALS' notation, the new axes in terms of the old:
//: "b+c,a+c,a+b" takes the primitive cell of a body-centred cubic lattice to
//: its conventional one. Coefficients may be fractions, "-1/2a+1/2b+1/2c".
class ChangeOfBasis {
public:
  static ChangeOfBasis parse(const std::string &text);
  //: The column of M for each new axis: new axis j = sum_i M(i, j) old axis i.
  Mat3 M;
  //: hkl in the new basis, M^T hkl. Refused, not rounded, if it is not
  //: integral: that means the new cell is not a supercell of the old.
  Miller apply(const Miller &hkl) const;
  //: The reciprocal matrix in the new basis, A M^-T; every scan point too.
  void apply(Crystal &crystal) const;
};

//: Every row's miller_index in the new basis, and the crystals of the
//: experiments; with the space group recorded on each crystal.
void reindex(ExperimentList &experiments, Table &reflections,
             const ChangeOfBasis &cb, const SpaceGroup &group);

} // namespace mxi
