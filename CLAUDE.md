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

**Thresholds come from measurement, not from comfort.** The real-data tests
assert against numbers that were printed first: parallax 0.96 px over the
sample, residual median 1.9e-4 and max 7.1e-4, agreement with DIALS' `rlp` at
4.3e-5. Median and maximum are asserted separately, because the residual
distribution has a tail and bounding only the maximum bounds only the tail.

**A fixture that cannot diffract is a silent test, not a failing one.** The
synthetic panel was at `+z`, consistent only with the wrong sign of s0.
Correcting the sign made it predict nothing, and the prediction tests passed
over an empty list. They now assert counts.

**Double, not float, throughout the geometry.** Positions are wanted to a
millipixel over a 4000-pixel detector, one part in 4e6, and float has seven
digits. Apple GPUs have no doubles at all, so any Metal port of prediction has
to demonstrate the precision it achieves against this reference rather than
assume it — the summation-order differences already seen in the spot finder's
centroids are the floor, not the ceiling.

**Prediction is the oracle.** Build order is prediction, then indexing, then
refinement, because prediction generates ground-truth reflection lists to test
the other two against.

## The convention trap, and what it caught

The closed-loop test **cannot** detect a wrong convention: predicting and
mapping back applies the error twice, once in each direction, and it cancels.
A green suite proves self-consistency only.

It hid two real errors until a real `.expt` arrived. `docs/conventions.md` has
the detail; the short version:

- **`s0 = -direction / wavelength`.** The sign was backwards. dxtbx points the
  beam direction back towards the source.
- **The parallax correction was missing entirely.** Worth 1.58 pixels on a
  0.45 mm sensor. Not a constant offset — a smooth function of scattering
  angle, which refinement partly absorbs and leaves as a radial residual.

Plus two structural surprises: the goniometer is multi-axis
(`axes`/`angles`/`scan_axis`, not `rotation_axis`/`fixed_rotation`) and
`scan.properties.oscillation` is a per-image array, not `[start, width]`.

`tests/test_real_geometry.cc` embeds forty real reflections and the geometry
that produced them, regenerable by `tests/make_real_data.py`. Both parallax
directions reproduce DIALS to the last bit; `entering` agrees on 100% of 13072
reflections; the full chain lands on `A h` with median residual 1.9e-4.

**Still open, and a green suite says nothing about any of them:** the
goniometer decomposition (every angle is zero in this dataset, so both
rotations collapse to the identity), `first_image != 1`, multi-panel
detectors, and the `hierarchy` block, which is ignored.

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
