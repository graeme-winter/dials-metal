# Symmetry

STATUS: `mxi_symmetry` works on one sweep. Not yet: more than one sweep, and so
the indexing ambiguity between them. `docs/outstanding.md` has the whole list.

`mxi_symmetry` finds the Laue group and space group of one sweep's integrated
data and writes the data reindexed into them, as dials.symmetry does.

## Running it

    mxi_symmetry integrated.expt integrated.refl     # symmetrized.expt, symmetrized.refl

`--max-delta D` sets the lattice's tolerance in degrees of obliquity (2). It
prints what it merged and to what resolution, each symmetry element's score,
every subgroup's, the space groups judged by their absences, and the choice.

## How it works, and against dials.symmetry

1. **The lattice's symmetry**, by gemmi's implementation of Le Page's method,
   within 2 degrees; its subgroups, from the closures of every pair of
   rotations; each named in its reference setting by a search over short
   lattice vectors along the group's axes -- and for a group of one axis those
   it reverses or that lie perpendicular to it -- since gemmi names a group
   only in a tabulated setting.
2. **Each element and each subgroup scored** as Evans (2011) A1 and A2, read
   from dials.symmetry's source: elements counted as it counts them, 17 for
   m-3m; E(CC; S) and sigma(CC) estimated as it does. Intensities merged in P1
   and quasi-normalised by resolution, where dials.symmetry fits an anisotropic
   maximum-likelihood model.
3. **The space group** from the absences, among the chiral groups with that
   Patterson group: consistent if what it forbids beyond its centring has mean
   I/sigma of 3 or less, and of those the one explaining most; a tie to the
   lowest number, the others reported as indistinguishable.

On the 300 image sweep, integrated here:

                          here              dials.symmetry
    Patterson group       I m -3            I m -3
      likelihood          1.000             0.879
      NetZcc              9.42              6.79
      CC, CC-             0.99, 0.04        0.89, 0.18
    reindex operator      b+c,a+c,a+b       b+c,a+c,a+b
    space group           I 2 3             I 2 3
                          (I 21 3 indistinguishable by absences under I)

The candidates are the same Patterson groups with the same multiplicities as
dials.symmetry lists; monoclinic ones are named C 1 2/m 1, the reference
setting, where dials.symmetry gives I 1 2/m 1. mxi_scale on symmetrized.expt
with no options gives merging statistics identical, and scaled.refl
byte-identical, to mxi_scale given the space group and change of basis by hand.
The cell is written reindexed but not constrained: 77.93 77.90 77.89 A where
dials.symmetry reports 77.898.

## A dataset measured to the corner: two simplifications undone

On a sweep dials.symmetry called I m -3, mxi_symmetry chose P -1: the true
three-folds scored CC 0.27 where dials.symmetry has 0.97. The lattice's rotations
were right -- they keep the cell's metric to 0.2 per cent, transposed they are
167 per cent off -- and the difference was in the data scored:

* **No resolution limit.** dials.symmetry scores to a limit from the data, the
  finer of CC half above 0.6 and <I>/<sigma> above 4 (1.68 A there); this scored
  to the detector's corner, 176220 reflections with mates merged where
  dials.symmetry kept 173252 with mates apart, most of them beyond the
  diffraction. Quasi-normalised, noise weighs as much as signal in a correlation.
  Now the same limit, with dials.symmetry's other filters: observations with
  I/sigma below -5, and Wilson outliers of E^2 16 or more.
* **Friedel mates merged.** The identity then compared each reflection with
  itself and gave CC 1 exactly, and E(CC; S), the CC a true element is expected
  to reach, was set from that. With mates apart the identity compares I(h) with
  I(-h) and measures it: 0.987 on the 300 image sweep, where dials.symmetry has
  0.979.

A test plants m-3 in strong reflections to 3.5 A and pure noise beyond: without
the limit it finds P -1, as that sweep did, and with it I m -3. The 300 image
sweep still gives I m -3, b+c,a+c,a+b and I 2 3. The limits differ from
dials.symmetry's -- 1.85 A from <I>/<sigma> there, 2.46 in dials.symmetry, which
fits a curve where this takes the last shell above the threshold.

After the fix that sweep gives I m -3 at likelihood 0.998, from 1.71 A where
dials.symmetry used 1.68, E(CC; S) 0.871 against 0.879, and I 2 3.

## The indexing ambiguity

On that sweep dials.symmetry reindexed by a+b,a+c,-b-c and this by b+c,a+c,a+b.
In the conventional cell dials.symmetry's axes are these as (c, b, -a): a
four-fold about b, which is in the lattice's m-3m and not in the crystal's m-3.
The two tables are on opposite sides of the indexing ambiguity m-3 has in a
cubic lattice. For one sweep either is right and scales the same; symmetrized.refl
is not index for index the same as dials.symmetry's, and sweeps -- or a sweep
and a reference -- must be put on the same side before they are merged. That
wants a reference, or dials.cosym's method.
