# dials-metal -- working notes

## What this is

The standalone Metal pipeline downstream of spot finding: indexing, refinement
and prediction. Reads `strong.refl`, writes `indexed.expt` and `indexed.refl`.

Three independent pieces share this repository and are not one program:
`spotfinder/`, the pipeline in `src/` and `apps/`, and `mxeq/`. Each has its own
notes; `spotfinder/CLAUDE.md` and `mxeq/CLAUDE.md` are theirs and this file does
not restate them.

`mxeq` must stay independent of everything else here. It is the referee, and a
referee that depends on the thing it judges is not one. Nothing in `src/` or
`apps/` may import it, and it must never import them; it reads files, which is
the whole point.

**The spot finder has not been built or tested in the environment these notes
were written in**, which has no HDF5 and no way to install it. Its CMake is
included unmodified and guarded at the top level, and the only path verified
here is the one where it is skipped. Anything said about it below comes from
reading it, not from running it.

Two overlaps, neither resolved:

* `spotfinder/src/refl.cc` writes reflection tables and so does `src/refl.cc`.
  The spot finder's is specialised to its Spot type, mine is a general reader
  and writer validated against real DIALS files. Consolidating them is worth
  doing and is exactly the kind of change that quietly breaks a format that
  currently works, so: not in the same commit as the merge, and with the spot
  finder's own output compared before and after.
* `spotfinder/src/expt.cc` reads experiment lists too, but only for the scan
  range, the panel size and the identifier, and its header says plainly that it
  is not trying to be dxtbx. So the geometry conventions still live in exactly
  one place, `src/geometry.h`, which was the thing worth checking before
  merging.

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

## Where the cell and distance difference from DIALS actually came from

Not the refinement weights, which was the obvious suspect. On this detector the
centroid variances are nearly constant -- median 0.086, 0.087, 0.085, with a
floor at 1/12, the variance of a uniform distribution over one pixel -- so
statistical weighting is very nearly unit weighting. Measured: unit weights
give a distance of 170.080 against 170.088 for statistical. Eight microns.

**Given DIALS' reflection set, this refinement reproduces DIALS.**

    refinement input          distance     cell
    DIALS' indexed.refl       170.0714     67.424 67.468 67.475
    DIALS' own answer         170.0730     67.425 67.474 67.468
    our indexed.refl          170.3452     67.526 67.560 67.605

Two microns in distance and six thousandths of an Angstrom in cell. The entire
discrepancy is in the indexing, not the refinement.

**What the indexer admits that DIALS does not:** 380 reflections, all of them
weak. Median intensity 24 against 241 for the rest, median n_signal 4 pixels
against 27. The 13064 both index agree on every single Miller index.

**Outlier rejection cannot remove them.** At 2, 2.5, 3 and 4 sigma the refined
distance is 170.348, 170.347, 170.345, 170.344 -- it does not move. That is the
part worth remembering: *a misindexed population that is internally consistent
does not produce outliers, it moves the model.* Refinement adjusts until those
reflections fit, and then nothing about them looks anomalous. Rejection only
catches reflections that disagree with their neighbours, and these agree with
each other.

Tightening the assignment tolerance does work -- 0.30 and 0.20 both give
170.35, 0.15 gives 170.076, 0.10 gives 170.038 -- but lowering the default is
the wrong conclusion, because DIALS' own `hkl_tolerance` is also 0.3. The
difference is that DIALS *iterates*: it assigns, refines, re-assigns against
the improved model, and discards what stops fitting. This indexer assigns once.

**The fix is a macrocycle**, not a tighter tolerance, and it is now in:
assign, refine on the stronger half, re-assign under the improved model,
repeat. Strength is measured against the dataset's own median so the criterion
travels between detectors and spot finders.

    stage                                distance    volume vs DIALS
    assign once, refine on everything    170.345        +0.55 %
    indexing macrocycles                 170.101        +0.09 %
    ... and mxi_refine --strong-only     170.092        +0.08 %
    DIALS                                170.073          --

It converges in one cycle on this data and `rmsd_index` falls from 0.0757 to
0.0318. The same option exists on `mxi_refine`, because refining on everything
afterwards puts the bias straight back -- 170.101 becomes 170.243. It is off by
default there: throwing away half the data should have to be asked for.

## Is the radial residual the same in DIALS and here? Partly

Compared directly, on the same 12896 reflections, both models' residuals binned
by radius:

     radius      DIALS     here        radius      DIALS     here
     50- 150    -0.270   -0.250       800-1000    -0.007   +0.011
    150- 250    -0.227   -0.214      1000-1200    +0.051   -0.003
    250- 350    -0.184   -0.156      1200-1500    +0.098   -0.082
    350- 450    -0.140   -0.096
    450- 600    -0.070   -0.031
    600- 800    -0.045   +0.006

**Shared and robust:** an inward radial displacement at low angle, reaching
-0.27 px nearest the beam and decaying monotonically to zero by about 800 px.
The two models agree to 0.02 px in the innermost bins. It is independent of
signal to noise -- -0.125, -0.108, -0.109, -0.117 px across a factor of 4.7 in
I/sigma -- so it is geometric, not a centroiding or background bias.

**Not shared:** everything beyond about 800 px. Excluding the innermost bin the
two profiles correlate at 0.12, and at 0.03 after the affine part is removed;
over the outer eight bins alone the correlation is -0.16, and at 1200-1500 px
the two models disagree in sign. Peak to peak, DIALS 0.235 px against 0.156
here.

**So the single "radial slope" figure describes neither model well.** It mixes
one real shared feature at low angle with refinement-specific structure at high
angle, and its value -- +3.0e-4 for DIALS, +1.09e-4 here -- depends mostly on
where each refinement put the detector. Earlier entries in this file quote that
slope as though it were one phenomenon. It is not.

What is genuinely unexplained is therefore narrower than previously stated: an
intensity-independent inward radial displacement of a quarter of a pixel within
about twenty degrees of the beam, decaying to nothing beyond it. Module tiling
is a separate matter and lives at all radii.


### It is a fixed hardware fingerprint

The steps repeat. Across four independently refined sweeps of a separate
four-sweep insulin experiment on the same detector, the slow-axis pattern comes
out the same shape every time:

    boundary   sweep 0   sweep 1   sweep 2   sweep 3
      1062      +0.045    +0.070    +0.022    +0.090
      1612      -0.187    -0.097    -0.128    -0.142
      2162      -0.321    -0.362    -0.366    -0.326
      2712      +0.354    +0.239    +0.330    +0.361
      3262      -0.319    -0.132    -0.317    -0.348

and comparing that experiment against the single-sweep one -- different data,
different crystal orientation, separately refined detector models -- gives a
**correlation of 0.950 over seven boundaries, with an rms difference of 0.063
px**. The steps themselves are 0.229 px rms, which at 75 micron pixels is
**17 microns**.

Seventeen microns is a module placement tolerance, not a modelling error. This
is a property of the hardware, measurable from diffraction data, stable between
experiments, and currently absorbed into the residuals because the detector is
described as one flat panel.

**The actionable part:** a multi-panel model, one panel per module, would give
refinement somewhere to put this. That is a change to how the NXmx file is
imported, not to any algorithm. The per-module affine result suggests the
modules want tilts and not just positions.

### The falsification test, passed

l-cysteine on a PILATUS 2M at I19: 1475 x 1679, 172 micron pixels,
3 x 487 + 2 x 7 fast and 8 x 195 + 7 x 17 slow, confirmed from the data by zero
spots in every predicted gap out of 15477 where arbitrary bands of the same
width hold 46 and 141.

**It shows no module steps.** Median absolute step at its nine real boundaries
is 0.030 px, against 0.038 at sixty arbitrary slow rows and 0.024 at forty
arbitrary fast columns. The per-sweep values flip sign between sweeps. Nothing
at the Eiger2 spacings either, as there should not be.

The test has the power to have found it. With around four hundred reflections
each side and a 0.56 px spread the standard error on a median step is about
0.035 px, so an Eiger2-sized step of 0.12 to 0.35 px would have shown at four
to ten sigma.

One thing had to be fixed first. As supplied, this dataset is indexed with one
UB across all four sweeps and perfect goniometry, which is false, and the
residuals are dominated by it: median 1.05 px, with only 5218 of 11686
reflections inside a 1.5 px clip. Refining with a crystal per sweep --
`mxi_refine --separate` -- takes it to 0.56 px with 10791 reflections, and only
then is the test sensitive enough to mean anything. **A null result from an
underpowered test is not evidence of absence**, and this one would have been
underpowered by a factor of two.

So the fingerprint is detector-specific: present and reproducible on the Eiger2
at 0.229 px rms (17 microns), absent on the PILATUS 2M below about 0.06 px
(10 microns). An artefact of the analysis would appear on both.

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

## Scan-varying refinement

Control points in A, evenly spaced over the scan, interpolated by a uniform
cubic B-spline clamped at both ends. Not DIALS' Gaussian smoother, which
weights three sample points per observation; both are smooth and local, and
they are not the same model.

It was linear at first. Cubic because a crystal does not change direction
abruptly at arbitrary points in a scan, and linear interpolation says it does --
continuous, but with a derivative that jumps at every control point. A B-spline
is C2.

B-spline rather than an interpolating spline because the support is local: four
control points per evaluation, so the Jacobian is banded rather than dense.
That is worth more than it sounds -- see `docs/gpu.md`.

The cost is that in the interior a control point is a coefficient, not the
value of the model at its own position. Only the two ends interpolate, by the
clamping. Anything reading an interior control point as "the setting matrix
there" is wrong.

What changed and what did not: on real l-cysteine the fit is identical, 0.2101
px at nine control points against 0.2100 for linear and DIALS' 0.2096. On
synthetic data with a planted drift the spline recovers it exactly at five
control points where linear needed three and could only approximate -- because
the clamped ends are pinned by the data near them, where linear's outermost
points wandered.

It matches DIALS. On l-cysteine, median |d| against number of control points:

    1 (static)  0.590      5  0.211
    3           0.252      9  0.210
                          15  0.217     DIALS (Gaussian smoother)  0.2096

Five to nine points reproduces DIALS to the fourth decimal, and fifteen is
worse -- the drift is captured and what is left is noise being fitted.

**Static first, always.** Control points started from an unrefined model
absorb errors that belong to the detector, fit well, and mean nothing.

**Outlier rejection is required, not optional.** About four per cent of
reflections have forward and reverse maps that disagree under a scan-varying
model. This was attributed to near-tangential geometry; **that attribution is
withdrawn**. The paper's own criterion for near-tangential -- the volume
(e x r_phi) . s0 below 0.05 -- removes 5.4 per cent of reflections and only 12
per cent of the disagreements. Nor is it convergence of the forward iteration,
nor reflections whose two Ewald roots are close. The population is
unexplained; rejection removes it and refinement then works, which is a
workaround. Once the model varies over the scan, the forward
and reverse maps select different roots for those, and their residuals run to
tens of images. With them in, a planted 0.5 degree drift is not recovered at
all; with them rejected it comes back as 0.5000. The fraction is 3.5 per cent
at a drift of 0.02 degrees and 4.3 at 0.5, so it is a property of those
reflections and not of how much the crystal moved.

**The outermost control points were weakly constrained under linear
interpolation**, by whatever lay at the very ends of the scan. Clamping the
B-spline fixed that: the end control points are on the curve, so the data near
the ends pins them directly.

## A bug the round-trip test could not see

The forward map (`predict`) and the reverse map (`centroid_residual`) must be
one model seen from two directions. For a scan-varying crystal they were not:
the iteration that finds the setting matrix was seeded from the start of the
scan for *both* Ewald roots, so a reflection late in the sweep converged on the
wrong one. Forward and reverse disagreed by 0.62 px rms, and since refinement
minimises the reverse map against data made by the forward one, no amount of
refining could reach the truth.

`prediction_round_trips_through_reciprocal_space` could not catch it: it goes
through `reciprocal_lattice_point`, not through the target function refinement
actually minimises. **Test the function being minimised, not a cousin of it.**
`prediction_and_the_refinement_target_agree_exactly` now does.

It was found by asking the dullest possible question -- does the truth model
reproduce its own predictions -- which should have been the first test written
and was not.

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

## Command lines take options in any position

All three programs originally read their file names from `argv[1]` and
`argv[2]` and scanned for options from `argv[3]` onwards, so

    mxi_refine --beam indexed.expt indexed.refl

took `--beam` as the experiment list, `indexed.expt` as the reflections, and
reported "unknown option 'indexed.refl'". The error named an argument that was
not the problem, which is worse than no message.

`src/args.h` splits a command line into positionals and options wherever they
appear, and is in `src/` rather than `apps/` so that it can be tested. Argument
parsing is exactly the kind of code that never gets tested because it looks too
simple to get wrong.

An unknown option is refused rather than ignored. Silently accepting a
misspelling is how a run comes to use settings nobody chose -- `--strong-ponly`
would otherwise have been dropped and the refinement would have used every
reflection while the operator believed otherwise.

## dials.* asks the flags, not the Miller indices

The `flags` column is a bitmask and every DIALS tool filters on it. Indexing
here set correct Miller indices and left the flags alone, so the output
processed perfectly and was then invisible to every selection downstream --
which is a failure mode worth naming, because nothing errored.

`mxi_index` now sets `indexed` (1 << 2) and `mxi_refine` sets
`used_in_refinement` (1 << 3), both cleared as well as set so that the flag and
the Miller index cannot disagree. An indexed insulin row now carries 36, strong
| indexed, which is what a real DIALS indexed.refl carries; after refinement
the fitted ones carry 44.

**`Table::int_column` REPLACES a column with a zeroed one.** That is right for
a derived column recomputed in full -- leaving stale predictions in the rows a
pass skips would be worse than clearing them -- and wrong for any column that
must be read before it is written. Setting the indexed bit that way threw away
the strong bit dials.find_spots had set, and nothing failed until something
filtered on it. `modify_int_column` returns the existing column;
`flags` is the one place that needs it.

`RefineResult::rows_used` reports which rows the fit actually used, after
outlier rejection and after the ill-conditioned ones were dropped, because
which reflections a residual was averaged over is not a detail and this is
where DIALS records it.

`centroid_outlier` is bit 17, settled from a real DIALS indexed.refl rather
than guessed: 8230 rows carry it, every one indexed, none also marked
used_in_refinement, and their median |xyzcal - xyzobs| is 0.843 px against
0.310 for the rest. The four values DIALS writes are 32, 36, 44 and 131108, and
this now writes the same four.

DIALS sets `used_in_refinement` on exactly 18000 of 75406 indexed reflections
on a 180 degree scan -- a hundred per degree, its sampling default. This sets
it on everything it fitted, which is a real difference in what the flag means
between the two and is not a defect in either.

## Hold the detector during the scan-varying pass

A scan-varying crystal and a refinable detector distance are degenerate in
position. Measured on a truth that fits exactly, both scaled by two tenths of a
per cent:

    distance and cell together     0.062  0.072  0.318
    distance alone                 0.778  1.100  0.000
    cell alone                     0.836  1.152  0.318

Twelve times smaller in position when they move together. The rotation angle is
*not* degenerate -- it sees the cell and not the distance -- which is why a
static refinement pins the pair, and why a scan-varying one does not: a hundred
and sixty crystal parameters can absorb the angular residual by drifting the
orientation, leaving the cell and the distance free to slide together.

On 1800 images of insulin, eighteen control points:

    --beam                       distance 169.977  V 236088   rmsd 0.309 0.295 0.320
    --beam --scan-varying 18     distance 169.705  V 234704   rmsd 0.233 0.218 0.193
    ... detector held            distance 170.010  V 236014   rmsd 0.231 0.234 0.192
    dials.refine                                   V 236092   rmsd 0.206 0.208 0.203

The distance drifted 0.27 mm and the volume fell 0.59 per cent, for five
thousandths of a pixel. Holding the detector brings the cell to within 0.013 per
cent of DIALS' edges and 0.03 per cent of its volume.

**The static pass has already placed the detector**, which is what makes
holding it safe rather than a constraint on a quantity nobody has determined. A
detector does not move during a sweep, so there is nothing scan-varying about
it; letting it move while the crystal is free only gives crystal drift
somewhere else to go. `--detector-in-scan-varying` restores the old behaviour.

A first attempt to test this refined a static truth with the detector free and
with it held, and asserted that holding it kept the cell closer. It does not,
and should not: with exact data and a detector starting in the wrong place,
refining it recovers the truth exactly and holding it forces the error into the
cell. The degeneracy is invisible there because there is a unique exact answer.
The test now measures the degenerate direction itself.

## The rotation angle has to be determined before it is worth fitting

Waterman eqn (40) divides by the volume of the parallelepiped formed by the
rotation axis, the reciprocal lattice vector and the beam. Near the rotation
axis it goes to zero: the angle at which such a reflection diffracts is
arbitrarily sensitive to the model, and the Lorentz factor has the same
asymptote, so its observed angular centroid is poorly determined too. Fitting
them puts noise into the rotation-angle residual that no model can remove.

DIALS discards below 0.05 by default. This did not, and the cost shows in the
rotation-angle residual first:

    insulin, 30 degrees, static     0.310 0.290 0.245  ->  0.274 0.254 0.220
    l-cysteine, 170 deg, sv 9       0.168 0.212 0.133  ->  0.161 0.194 0.124

The fraction removed depends on the geometry, not only on the cutoff: 2.4 per
cent of indexed insulin over thirty degrees, 0.1 per cent of l-cysteine, whose
small cell puts every reflection far from the axis.

**A residual averaged over a different set of reflections compares with
nothing**, so `mxi_refine` now says how many it dropped and how many it
averaged over. dials.refine applies the same cutoff and reports over what is
left, which is part of why its numbers looked better than they were.

## Which model each observed column was computed through

`dials.refine` stops without `xyzobs.mm.value`, and reasonably: DIALS measures
its residual in millimetres and radians -- Waterman eqn (25) is in X, Y and phi
-- so that column is the observation it minimises against. `dials.index`
produces it and this did not.

The columns split by **which model they were computed through**, and the split
was measured against DIALS' own output rather than assumed:

    xyzobs.mm.value       the model as IMPORTED    median difference 3e-17
    xyzobs.mm.variance    the model as imported    essentially exact
    s1, rlp, entering     the CURRENT model        3e-3 from imported,
                                                   9e-5 from DIALS' refined

`dials.find_spots` writes the millimetre centroids through the detector as
imported and nothing recomputes them; `dials.index` calls
map_centroids_to_reciprocal_space after refining, so s1 and rlp follow the
refined model. **That is why xyzobs.mm.value is stale after refinement and must
never be a join key** -- a fact recorded here long before the reason for it was
understood.

Recomputing the millimetre centroids from the refined model would look like a
correction and would be a silent change to the observations refinement had just
been fitted to.

`entering` is `s1 . (m2 x s0) > 0`, the opposite sign to the test as usually
written, because s0 here points source to sample and dxtbx's beam direction
points the other way. Determined against DIALS' flags: 13760 of 13766 agreed,
the rest having a triple product within rounding of zero. The wrong sign gave
six.

## Write back everything that was there, not only what we model

`dials.refine` then failed with

    IndexError: list index out of range        experiment_list.py:603

because the experiment said `"imageset": 0` and the imageset list had been
written out empty. The imageset block is the only link from an experiment list
to the images; `profile`, `scaling_model` and `history` are equally not this
package's to discard.

`ExperimentList` now keeps the document it was read from, and writing replaces
the models it understands and leaves everything else alone. A list built in
memory has no source, and then the absent models are **null** -- never `-1`,
which is a valid index into a Python list and reaches for the last element of
one that may be empty.

Carrying a block through then exposed the other half of the number problem
below: the reader did not record whether a value had been written with a
decimal point, so an `"imageset": 0` came back out as `"imageset": 0.0` and
dxtbx indexed a list with it. Both halves are needed, and each is useless
alone.

## `0` is not `0.0`, and dxtbx knows the difference

`dials.refine` refused an `indexed.expt` written here:

    DXTBX_ASSERT(obj_type == "float") failure      dxtbx scan.cc:80

A scan starting at zero degrees has an oscillation array beginning 0.0, and the
JSON writer collapsed any whole-valued double to an integer, so the file said
`[0, 0.1, 0.2, ...]`. dxtbx reads the type of an array from its first element
and requires float.

The fix is not to write every number with a decimal point -- `image_range` of
`[1.0, 300.0]` would be just as wrong. Whether a number is a count or a
measurement is known at the point it is constructed and nowhere else, so
`json::Value` now carries it: built from `int`, `long`, `long long` or
`std::size_t` it writes without a point, built from `double` it writes with one.

**Nothing here could have caught this, and that is the recurring lesson.** The
document round-tripped through this reader perfectly, because this reader parses
`0` and `0.0` into the same double. A round-trip test cannot detect a wrong
convention; only a file written or read by something else can. That is now five
times: the `.refl` container shape, the geometry conventions, the scan
oscillation form, the number types, and the blocks dropped on write.

`tests/test_expt_format.cc` therefore asserts on characters rather than values,
and was checked by reverting the fix -- three of its tests fail without it.

## Relative precision belongs to the arithmetic, not to the answer

This file once claimed that float32 would leave about four digits in a
numerical derivative of the target. It leaves none. The reasoning was that a
1e-6 relative parameter step changes the residual by 1.5e-3 of its own size and
float epsilon is 1.2e-7, so four digits survive.

The residual is a difference of detector positions of order two thousand
pixels. Its absolute error in float32 is epsilon times the *position*, about
2.4e-4 px, not epsilon times the residual. The change being measured is about
5e-4 px. Signal and noise are the same size.

Measured rather than estimated, in `tests/test_precision.cc`: median relative
error 1.00 at the step the refinement uses, and 1.4e-2 at the best step float
can manage. The analytical derivative in the same precision gives 1.3e-7.
Analytical derivatives are therefore a precondition for a device port and not
an optimisation.

Both `target.h` and `derivatives_t.h` are templated on the scalar type for this
reason, and both are checked against the double-precision originals first --
the derivatives bit for bit. Writing a second implementation to measure the
first is worse than useless: a disagreement could be either thing.

## Style

Row-major Mat3, because that is how DIALS serialises a matrix into a flat nine,
and a transpose hidden in the I/O layer is the easiest possible day to lose.

Tests register themselves by static constructor, so there is no list to forget
to add to. Comments explain why, not what.
