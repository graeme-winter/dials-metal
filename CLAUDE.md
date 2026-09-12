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

The goniometer decomposition is now closed too, against the l-cysteine
four-sweep data in `tests/real_cysteine.h`. Insulin could never have done it:
every setting angle there is zero, so both rotations are the identity and any
arrangement passes.

**Test the wrong versions, not just the right one.** `test_multi_axis.cc`
asserts that swapping fixed and setting is caught by exactly two sweeps, that
dropping the fixed rotation is caught by exactly two, and that folding the
scan-axis angle into the fixed rotation is caught by three -- including sweep
0, which is blind to the other two failures and is the only thing that pins
that rule. Without those, a green suite would not distinguish a correct
composition from a lucky one.

**Goniometer axes run from the sample outwards to the laboratory.** `axes[0]`
carries the sample; each subsequent axis carries the one before it; the last is
bolted to the floor. So an axis further out applies later and multiplies on the
**left**. The first implementation accumulated the other way, composing the
stack inside out -- and all 48 tests passed, because with one axis below the
scan axis the two orders are the same matrix. Pinned now against the physical
arrangement on a synthetic three-axis goniometer, including a check that the
two orders really do differ so the test is not vacuous.

**"Three circles" is not the condition for testing the composition order.** A
four-sweep phi/chi/omega dataset does put two axes below the scan axis, but phi
is zero in every sweep, so its rotation is the identity and both orders give
the same matrix to the last digit. What is needed is two axes below the scan
axis at *simultaneously non-zero* angles. Asserted, so that substituting better
data makes the test fail rather than silently start passing for a new reason.

**Still open, and a green suite says nothing about any of them:** the
composition order *against data*, `first_image != 1`, multi-panel detectors,
and the `hierarchy` block, which is ignored.

## Indexing: three things that produced plausible wrong answers

**Peak centroids are not a basis.** The FFT grid step is around 0.8 Angstrom,
so vectors straight off the transform are good to a few parts in a thousand,
which is a tenth of an index at the detector edge. The least squares fit of A
to the indexed reflections is what makes the cell numerically right rather than
merely recognisable. It tightens its tolerance over four rounds; fitting
throughout at the acceptance tolerance lets reflections a quarter of an index
out pull the basis as hard as ones that are exact.

**The fraction indexed is basis-dependent and cannot compare two bases of the
same lattice.** Acceptance asks whether every fractional index is within
tolerance, and a unimodular change of basis mixes the components: (0.2, 0.2,
0.2) in one basis is (0.4, 0, 0.2) in another describing the same lattice. A
guard that only accepted the reduction if it indexed at least as many
reflections duly rejected the correct cell. The invariant to test is the
volume.

**Coincident reciprocal lattice points are input, not noise.** Pooling several
sweeps puts the same reflection, measured at different goniometer settings, at
the same place in the crystal frame — that is the point of measuring more than
one sweep. Counting those pairs in the nearest-neighbour cell estimate made it
diverge to 9e15 Angstrom and find no candidates at all.

## Resolved: the 0.45 per cent cell edge

Refinement fixed it, and the explanation is the detector. Before refinement two
edges matched DIALS to 0.02 per cent and one was 0.45 per cent short; after
refining crystal and detector together all three agree with each other and the
asymmetry is gone. Refining with the detector *held fixed* leaves the asymmetry
in place, which is the evidence: the imported detector model was wrong, and the
cell was absorbing the error anisotropically because indexing had no other
parameter to put it in.

## The degeneracy that replaces it

Detector distance and cell scale are very nearly degenerate. A detector further
away with a proportionally larger cell predicts the same spots in the same
places, so the two refinements settle at different points on a flat valley:

    distance  170.34 mm (here)  vs  170.07 (DIALS)   0.16 per cent
    cell      +0.18 per cent per edge, +0.55 in volume

which is arithmetically the same statement twice. Both fit the observations
about equally well, and through the reindexing operator the two models predict
the same positions to 0.13 px median. So this is not a bug to find; it is a
direction the data barely constrains. Anything that claims to have fixed it
should be checked against whether the residual actually improved.

Refining the beam as well moves it by almost nothing (170.340 against 170.344),
so the beam is not what breaks the degeneracy.

## Refinement notes

**Derivatives are numerical, deliberately.** Fifteen parameters over thirteen
thousand reflections is sixteen predictions per iteration and takes a second. A
wrong analytical derivative does not crash — it converges smoothly to the wrong
answer and reports a small residual doing it. If profiling ever demands the
analytical version, this one stays as the thing it is checked against.

**Choose the Ewald root by proximity to the observation, not by the `entering`
flag.** The flag assumes a model already close enough to trust, which at the
start of refinement it is not. Choosing wrongly puts a reflection tens of
images away, so there is a test for it.

**Rotate the detector about the panel centre, not the laboratory origin.** A
rotation about a point 200 mm away is mostly a translation, and the two
parameter groups would then be so correlated that the normal matrix is nearly
singular.

**Reject outliers between macrocycles, never within one.** Rejecting while the
model is still moving throws away reflections for being far from a prediction
that was wrong.

## Planned: how indexing should handle multiple sweeps

Index across **all** sweeps at once, then immediately split and refine against
each sweep individually, keeping the bulk matrix common.

The reason is visible in the l-cysteine residuals. A joint index must assume
one UB and perfect goniometry, and the goniometer does not return to precisely
the same place between sweeps, so no single matrix fits all four -- residuals
come out a hundredfold worse than insulin. That is not a failure of indexing;
it is the constraint doing what it must. Joint indexing is still the right
first step, because it is what guarantees a consistent basis across sweeps and
avoids four independently-chosen and mutually reindexed lattices. The
constraint is then broken at the first opportunity rather than carried forward.

Consequence for the code: the indexer works on the pooled reciprocal lattice
points from every sweep, but each sweep keeps its own goniometer, scan and
detector throughout, and the output is a list of experiments sharing a crystal
rather than one experiment.

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
