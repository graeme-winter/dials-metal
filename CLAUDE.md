# dials-metal-index -- working notes

## What this is

The indexing, refinement and prediction stages of the standalone Metal
pipeline. Reads `strong.refl`, writes `indexed.expt` and `indexed.refl`.
Separate from `mxeq`, which is the referee and must never depend on a pipeline
it judges.

## Hard constraints

**No third-party dependencies.** Same rule as the spot finder's hand-rolled
msgpack writer: this builds on a laptop and on a beamline machine with nothing
installed. Vec3, Mat3, the Cholesky solver and the test harness are all here
for that reason, and each is under a hundred lines.

**Double, not float, throughout the geometry.** Positions are wanted to a
millipixel over a 4000-pixel detector, one part in 4e6, and float has seven
digits. Apple GPUs have no doubles at all, so any Metal port of prediction has
to demonstrate the precision it achieves against this reference rather than
assume it — the summation-order differences already seen in the spot finder's
centroids are the floor, not the ceiling.

**Prediction is the oracle.** Build order is prediction, then indexing, then
refinement, because prediction generates ground-truth reflection lists to test
the other two against.

## The convention trap

Every geometry convention is unvalidated. `docs/conventions.md` is the list.

The important part: **the closed-loop test cannot detect a wrong convention.**
Predicting and then mapping back applies the error twice, once in each
direction, and it cancels. A green suite here proves self-consistency only.

This is the same failure that hit the `.refl` reader in `mxeq`: seventy tests
passed against a format that was wrong in two ways, because the only check was
a round trip against its own writer. Only a file written by something else
proves a convention.

Required to close it: one real `indexed.expt` + `indexed.refl` pair. The
highest-value single check is mapping real `xyzobs.px.value` into reciprocal
space and comparing against `A * miller_index` — millipixel agreement validates
five of the seven conventions at once.

## Things got right for a reason

**The setting matrix is the plain inverse of the real-space matrix.** Real
space vectors are the rows, so their product with A is the identity. The
inverse transpose passes on any cubic cell because both are diagonal; the test
uses a triclinic cell for exactly this reason. Same bug, same fix, as in
`mxeq`.

**Panel intersection has a nanopixel edge tolerance.** A ray aimed at the exact
corner of a panel returns -4e-13 pixels, because the intersection is solved in
millimetres and divided back. A bare `< 0` test rejects it. That loses
reflections predicted on an edge, and on a tiled detector loses rays striking
the seam between two panels, which then belong to neither. A nanopixel is nine
orders below anything physical and three above the round-off.

**The goniometer inverse negates the angle rather than inverting the matrix.**
Exact for the rotation part, and avoids a determinant on a matrix that is
orthogonal by construction. It does assume `fixed` and `setting` are
orthogonal, which is true for goniometers and is not checked.

**`Mat3::inverse` sets a flag instead of returning infinities.** A caller that
ignores it gets an answer that looks wrong, rather than NaNs that propagate
silently through a refinement and come out as a plausible-looking cell.

## Style

Row-major Mat3, because that is how DIALS serialises a matrix into a flat nine,
and a transpose hidden in the I/O layer is the easiest possible day to lose.

Tests register themselves by static constructor, so there is no list to forget
to add to. Comments explain why, not what.
