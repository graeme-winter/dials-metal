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

## Not the parallax correction, and how that was settled

A natural suspicion about the cell and distance discrepancy is that the
millimetre-to-pixel mapping is missing from the prediction path: omitting it
would make predictions fall short radially, and refinement would push the
detector further out to compensate, enlarging the cell. The magnitude fits --
the correction is 0.048 mm median, which at a typical 50 mm radius is about a
tenth of a per cent in distance, against the 0.16 observed.

It is not that. Predicting from DIALS' own refined model reproduces DIALS' own
`xyzcal.px` **exactly**, 0.0000 px median and maximum. Turning the correction
off displaces predictions by 0.68 px median, 1.0 px worst.

The test that should have settled this immediately had a tolerance of 0.5 px,
which is less than the effect it was guarding against -- so it would have
passed for more than half the reflections with the correction missing
altogether. It now asserts exactness, and a companion test measures what
switching the correction off does, so the bound is known to be sensitive to the
thing it guards. **A tolerance looser than the error it exists to catch is not
a weak test, it is not a test.**

## Residual structure: what is measured, and two claims retracted

There is real structure in the post-refinement residuals at the 0.2 px level.
Binned eight ways with 1600 reflections each, the medians are far above the
noise on the mean. `mxi_residuals` prints it.

**Retracted: "x and y disagree about the distance by 0.6 mm."** That came from
scanning the detector distance while holding the *whole panel rigid*, so the
crystal had to absorb every lateral and angular error too. Measured properly --
fitting the radial residual against radius separately along each axis, with the
panel refined -- the implied distance differs between fast and slow by 0.04 mm,
not 0.6. The 0.6 was an artefact of the scan's own constraint.

**Partly retracted: "a clean radial gradient" on this refinement's own
residuals.** There the tangential component mirrors the radial almost exactly,
which is a directional pattern pushed through a radial decomposition. On DIALS'
own output, with the affine part removed, the radial gradient is clean and the
claim stands -- but it stands on DIALS' residuals, not on this code's, and the
two should not have been conflated.

**Refuted: the absorption depth.** Eqn (6) is the unconditional first moment
and counts photons that pass through the sensor as contributing depth zero. The
conditional mean, dividing by 1 - exp(-mu t), is 24 per cent larger at normal
incidence and varies with angle, so it cannot be absorbed into the distance.
Implemented behind `parallax_conditional` and tested: rmsd goes from
0.3514/0.3388 to 0.3516/0.3409 and the cell moves by 0.015 per cent. It is not
the cause.

## The radial residual, measured properly

Measured on DIALS' own refined insulin model, 12907 indexed reflections, so it
is DIALS' model inadequacy with none of this code's mixed in. `mxi_residuals`
prints all of it.

**It is real and it is radial.** Median radial residual runs monotonically from
-0.155 px at the beam centre to +0.067 px at the edge, with a fitted slope of
+3.0e-4 px per px, about 0.48 px across the radius range.

**No change to the detector can remove it.** The affine part of the residual
field -- two scales, a rotation about the beam, a shear, a translation, which
is everything a flat panel can express -- comes to almost nothing: scales of
0.01 per cent, 0.07 mrad of rotation, 0.035 mrad of shear, worth 0.054 px at
the panel corner. Removing it leaves the radial gradient essentially untouched,
-0.155 to +0.067. Refinement was free to move the distance and did; the
gradient is not a distance error.

**It is geometric, not a centroid artefact.** Fitting the slope within each
intensity quartile separately gives 3.4, 2.9, 3.3 and 4.7e-4 across a factor of
fifty in intensity. A centroid-estimation bias -- thresholding, background,
shoebox truncation -- would scale strongly with signal to noise. It does not.
The mild rise in the strongest quartile tracks bounding-box depth (3.2e-4 at
one image, 5.0e-4 at seven) and is plausibly second order.

**Not the absorption depth.** See below; tested and refuted.

**Not fluorescence escape either.** Si K-alpha is 1.74 keV. A 13 keV photon
that loses it still deposits 11.26 keV, far above a threshold set at half the
photon energy, so it is counted at its own pixel -- and the escaped 1.74 keV
photon falls below threshold and is never counted at all. A photon-counting
detector with a half-energy threshold is immune to this by construction.

**Partly the module tiling.** This is an Eiger2 16M: 4 x 1028 + 3 x 12 = 4148
fast, 8 x 512 + 7 x 38 = 4362 slow, confirmed against the data by zero spots in
every predicted gap out of 13766. DIALS imports it as ONE flat panel, so
nothing about the individual modules' positions or tilts can be represented,
and whatever is wrong with the tiling has nowhere to go but the residuals.

`mxi_residuals --modules 1028,12,512,38` shows it as a step at each boundary:

    fast at 1028   step -0.127 px
    fast at 2068   step -0.203
    fast at 3108   step -0.121

Three boundaries out of three, same sign, 0.12 to 0.20 px, each many times the
standard error on medians of hundreds. The slow boundaries are four of five
negative and noisier.

How much it accounts for, with controls, on median |residual| of 0.2795 px:

    per-module constant            0.2412  (13.7 per cent)
    per-module affine              0.1856  (33.6)
    control, diagonal bands        0.2508  (10.3)
    control, random groups         0.2800  (-0.2)

The control matters. A per-module constant barely beats carving the detector
into the same number of arbitrary diagonal bands, so rigid module offsets are
not the right description -- it is mostly absorbing smooth spatial structure.
The per-module *affine* is much better, which points at each module having its
own tilt rather than its own position.

Removing per-module constants also drops the radial slope from 3.0e-4 to
1.9e-4, so tiling is about a third of the radial gradient and something else is
the rest.

**The actionable part:** this detector is modelled as a single panel. A
multi-panel model, one panel per module, would give refinement somewhere to put
this. That is a change to how the NXmx file is imported, not to any algorithm.

So: a radial, geometric, intensity-independent displacement of about a fifth of
a pixel peak to peak, in a direction the detector model cannot represent.
Candidates not yet tested: the parallax offset is evaluated at the front-face
intersection rather than solved self-consistently (matching DIALS exactly is no
defence -- both could be wrong the same way, though the one-shot to converged
difference measured 0.0005 px, which is too small); a non-planar sensor; and
the interaction depth distribution being something other than a single
exponential, which for silicon near an absorption edge it is not.

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
