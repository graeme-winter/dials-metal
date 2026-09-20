# Integration

Not started. This is the plan, what it is bound by, and what is missing before
parts of it can be written at all.

## Where this stands

Not started. The profile model is implemented and the tools to judge it exist;
nothing integrates yet.

| | |
| --- | --- |
| profile model | implemented, `mxi_profile`; does not yet match DIALS |
| the region | implemented, `mxi_mask`; box and ellipsoid |
| forward model | implemented, `mxi_forward`; agrees with the data to four per cent |
| background | not started; needs Parkhurst et al. (2016) |
| summation | not started; needs Leslie (1999) |
| profile fitting | not started; needs the reference profiles of Kabsch section 3.3 |

Neither sigma agrees with DIALS. **Run `mxi_profile` for the numbers**: it
prints both, and the disagreement with each, every time. They are not repeated
here or in the README, because a number copied into prose is wrong from the
first time either estimator changes and nothing checks it.

What the forward model says about them is the most useful reading: rendered
onto the pixels, the observed spot is about three per cent narrower than the
model on the detector and five per cent narrower in rotation, so both sigmas
are a little large and by far less than the raw comparisons suggested.

Read `## Rendering the model onto the pixels` before anything below it. Several
earlier sections reach conclusions it overturns, and they are marked.

## What this follows

DIALS, transcribed, as the spot finder transcribes `DispersionExtendedThreshold`
window for window. That decision was taken earlier and it still holds: a
transcription can stand in for `dials.integrate`, while an integrator built on
this project's own ideas would be a fast pipeline that agrees with nothing
anyone runs.

Its cost is real and worth restating. Every entry in the diff checklist from
`dials_algorithm_map.md` becomes a specification rather than an interesting
divergence -- the 2D-disc foreground mask, the inert `filter.threshold`
parameter, LP and QE stored as columns and applied at scaling rather than at
integration, the GLM background. Some of those are things a research pipeline
would deliberately do differently. Divergences go in as run-time options that
default off, each with its reason recorded, so that a disagreement with DIALS is
always a choice and never a drift.

## What integration is, per Winter et al. (2018), §3.4

Three steps, in order:

1. **Profile parameters.** A three-dimensional Gaussian in a local
   reciprocal-space frame, with two parameters: `sigma_D`, the extent on the
   detector face, and `sigma_M`, the extent over images. Estimated from the
   indexed strong spots, per Kabsch (2010a).
2. **Background.** Shoeboxes are read from the images, and the background under
   the peak is modelled from the non-peak pixels around it. The default is a
   robust generalised linear model that assumes the counts are Poisson, which
   matters below one count per pixel where a normal approximation is biased
   (Parkhurst et al., 2016).
3. **Intensity.** Summation of background-subtracted pixels in the peak region,
   with errors from Poisson statistics (Leslie, 1999); and profile fitting,
   where the shoebox is transformed into the local reciprocal-space frame and
   fitted against a reference profile (Kabsch, 2010a).

Two details from the paper that are choices rather than consequences, and so
have to be copied deliberately:

* **DIALS differs from XDS in the transform.** Counts are distributed onto the
  reciprocal-space grid by computing the overlap of each detector pixel with
  the transformed grid point using Sutherland–Hodgman polygon clipping, rather
  than by XDS' assignment. This changes the profile and therefore the fitted
  intensity.
* **Blocks overlap by half.** Images are integrated in blocks whose start is
  aligned with the centre of the preceding block, so that most reflections are
  whole within one block. Reference profiles are built per block, at several
  points across the detector, each strong reflection contributing to its
  nearest profiles with a Gaussian weight by distance.

## What is still missing

Winter et al. (2018) defers to others for everything numerical. Two of those
are still not here, and each blocks one stage:

| | for | blocks |
| --- | --- | --- |
| Parkhurst et al. (2016) | the Poisson GLM background | background |
| Leslie (1999) | summation error estimates | summation |

Kabsch (2010a) and DIALS' own `calculator.py` are here, and between them they
settled the profile model. From the paper:

* §2.3, the `{e1, e2, e3}` frame and the mapping of a pixel to `(eps1, eps2,
  eps3)`, with `zeta = m2 . e1` correcting for the path length through the
  Ewald sphere;
* §3.1, the reflection mask `|eps1| <= delta_D/2`, `|eps2| <= delta_D/2`,
  `|eps3| <= delta_M/2`, and the estimators: `sigma_D^2` as the mean of the
  per-spot variances of the intensity-weighted beam directions, and `sigma_M`
  by maximising the likelihood of the observed offsets under `R(Delta,
  sigma_M/zeta)`;
* §3.3, the profile grid, the `f_3j` fractions that split a frame's counts
  between grid planes, and the 5x5 subdivision of each pixel that DIALS
  replaces with polygon clipping;
* §3.4, the fitted intensity `I = sum (c - b) p / v / sum p^2 / v`, with
  `v = b + I p` iterated from `v = b`, three cycles.

**Pixels: this package had them all along and was skipping them.** A
`strong.refl` carries a `Shoebox<>` column -- 13.6 MB of it for insulin -- and
the reader dropped it as an undecoded type. `src/shoebox.h` now decodes it. The
layout was derived from a real `dials.find_spots` file and checked against all
13766 of its records:

    int32          panel
    int32 x 6      bbox as x0, x1, y0, y1, z0, z1, half open
    uint8          a flag, 2 in every record seen
    float32 x N    data          N = (x1-x0)(y1-y0)(z1-z0)
    uint8  x N     mask          0 and 5, which is Valid | Foreground
    float32 x N    background    all zero out of dials.find_spots

Every bounding box agrees with the table's own `bbox` column and the records
consume the blob to the byte: 1467208 voxels holding 14163216 counts.

So the profile model can be estimated here, on real data, with no images and no
HDF5. What still needs the images is integration proper -- the shoeboxes of the
*predicted* reflections, which are a superset of the strong ones and are not in
any file this package has.

## sigma_D and sigma_M, against DIALS' own source

    mxi_profile refined.expt refined.refl

With `calculator.py` to hand the recipe is no longer guesswork. Three things it
settled, all of which had been wrong here:

**The reflections.** DIALS does not use every strong spot. It selects those
flagged `used_in_refinement`, then cuts on `|zeta| >= 0.05`. On 1800 images of
insulin that is 70425 of 78618, and it moves sigma_b from -4.2 per cent of the
DIALS value to -2.9 and sigma_M from a factor of three to twenty per cent.

**zeta uses the CROSS product.** `e1 = s1 x s0`, normalised -- the axis about
which the point would cross the Ewald sphere by the shortest route, Kabsch
section 2.3 after Schutt & Winkler. This had `s1 - s0`, which is the reciprocal
lattice vector and points somewhere else entirely. Correcting it moved sigma_M
from +50 to +21 per cent.

**Two places DIALS departs from Kabsch**, both followed here because standing
in for `dials.integrate` is the point:

* the background is NOT subtracted before the angular spread is measured, where
  Kabsch section 3.1 step (v) says to. DIALS' source carries a note saying so.
  It costs nothing on a table out of `dials.find_spots`, where the background is
  zero, and it would matter on one where it is not;
* an image contributes to the reflecting range if the spot finder marked any
  pixel on it as valid foreground, whatever the counts. Requiring counts as
  well is a different criterion that happens to select the same images here.

DIALS' `R` is a partiality rather than a density -- it is not divided by the
oscillation width as Kabsch writes it -- but that is a constant in the log and
the argument of the maximum is identical.

Both still disagree; `mxi_profile` prints by how much. The figure quoted here
was 0.030786 for sigma_b until the pixel-to-millimetre conventions were matched
two sections below, which moved it -- which is the argument for not quoting it
at all.

The remaining candidate is in code this has not read: `Shoebox::beam_vectors`,
which decides exactly which lab coordinate a pixel maps to. Parallax applied
the other way is not it -- that moves sigma_b to +24 per cent and inverting it
to -14, so DIALS is doing neither.

The estimators themselves are checked against planted values rather than
against DIALS: a known angular spread comes back exactly, and samples drawn
from the reflecting-range model with a known sigma come back within five per
cent at 0.05, 0.1 and 0.3 degrees.

## Where the centre is, and the half pixel that was hiding in it

Every measurement here maps a pixel into the frame of a reflection:

    lab   = origin + mm_fast * fast + mm_slow * slow
    s'    = lab * |s1| / |lab|                       elastic, so |s'| = |s1|
    eps1  = degrees( e1 . (s' - s1) / |s1| )
    eps2  = degrees( e2 . (s' - s1) / |s1| )
    eps3  = degrees( zeta * (phi_image - phi_calculated) )

The inputs are the panel geometry, the beam, the rotation axis, `s1`, and
`phi`. What was wrong was the first line: `mm` came from a plain multiplication
by the pixel size, while the `s1` it is compared against was built from
`xyzobs.mm`, which is parallax corrected. Two conventions, differing by about
half a pixel radially.

A width cannot see that, which is why it survived every test here. Comparing
each spot's centroid with its own centre in the frame does see it:

    eps1   mean -0.000007 deg   -0.000 sigma      sd 0.014 sigma
    eps2   mean +0.013506 deg   +0.439 sigma      sd 0.110 sigma
    eps3   mean +0.000176 deg   +0.001 sigma      sd 0.127 sigma

A systematic 0.44 sigma in the radial direction and nothing in the other two --
which is the signature of a radial displacement, not of anything physical.
Mapping the pixel through `px_to_mm` instead takes all three to zero within a
millionth of a degree.

**It costs agreement with DIALS.** sigma_D moves from 0.0308 to 0.0274 degrees,
from -2.9 per cent of the DIALS value to -13.6. That is recorded rather than
tuned away: a mapping that is demonstrably inconsistent cannot be the right one
however well its number happens to agree, and the remaining difference is now a
cleaner question than it was.

It also means `s1` in a DIALS reflection table is the OBSERVED scattering
vector, not the predicted one -- computed from `xyzobs.mm`. So these
measurements are centred on the observation, and the centroid sitting at zero
is a check on the arithmetic rather than a result. The mask in `mxi_mask` is a
different matter: it is built on the PREDICTED `s1`, which is what it must be,
since the point of a mask is to say where the model expects the signal.

## Rendering the model onto the pixels

    mxi_forward refined.expt refined.refl --out forward.txt

Every comparison before this one took a number from the data through a pixel
grid that truncates and quantises it, and a number from the model in closed
form. At widths below a pixel that difference is most of what was being
measured, which is why every conclusion carried the same caveat.

This renders instead. The model is integrated over the same pixels, the same
images and inside the same mask as the observation, and both are reduced the
same way, in pixels and images. Whatever the grid does to the data it does to
the model. What is left is the model being wrong.

The answer is that it is barely wrong at all, over 61047 spots:

    sensor model            shift fast   shift slow   obs/model fast  slow    z
    depth distribution        -0.008      -0.002          0.967     0.964   0.948
    mean depth, no smear      -0.008      -0.002          0.979     0.987   0.946
    no sensor                 -0.007      -0.051          0.984     1.004   0.946

**The positions agree to a hundredth of a pixel** once the sensor is in.
Leaving it out leaves a systematic 0.05 pixel shift in the slow direction and
nothing in fast, which is the parallax displacement and is the size it should
be.

**The widths agree to within four per cent**, not the factors that every
earlier measurement suggested. The captured-fraction test said one sigma held
0.74 where a Gaussian holds 0.47; rendered onto the grid, the observed spot is
3.3 per cent narrower than the model on the detector and 5 per cent narrower in
rotation. Nearly all of that apparent disagreement was the comparison, not the
model.

Taken at face value it says `sigma_D` should be about 0.0265 rather than 0.0274
and `sigma_M` about 0.113 rather than 0.119 -- small corrections in the same
direction as everything else, and much smaller than they looked.

### And the anisotropy is predicted

    anisotropy slow / fast    observed 1.074
                              model    1.092   with the depth distribution
                              model    1.076   with the mean depth only
                              model    1.062   with no sensor

The model predicts an anisotropy of the same size as the observed one, and
gets closer as more of the sensor goes in. So the elongation is geometry and
absorption acting on an isotropic Gaussian, not a missing physical effect: the
projection onto a flat detector and the depth at which a photon stops are
enough to produce it.

That is the answer to a question asked several ways in this document and
answered wrongly twice. The cubic crystal does argue for isotropy, and the
spots are consistent with an isotropic model -- once the model is compared with
the data on the data's own terms.

## Superseded: the anisotropy hunt

Everything from here to the end of this section was measured by comparing a
number taken from the data with a number taken from the model in closed form.
At widths of two pixels that comparison is dominated by the grid, and its
conclusions did not survive rendering the model onto the pixels instead: the
spots are consistent with an isotropic Gaussian once model and data are reduced
the same way.

It is kept because the reasoning is sound given what was measured, and because
three separate hypotheses were tested and rejected on evidence that turned out
to be an artefact of the method. That is worth being able to recognise again.

### The anisotropy, remeasured

With the pixel-to-millimetre conventions matched, over 61047 spots with more
than fifty counts:

    in the reflection frame   tangential 0.017363   radial 0.021149   ratio 1.237
    in the laboratory         along axis 0.018159   across   0.019496  ratio 1.051

**It survives the fix.** The earlier figure was about 1.5 from the aggregate
profile and 1.20 from the per-spot moments; it is now 1.24. The half-pixel
convention error was not the cause.

**It is not the sensor.** The ratio is flat with obliquity -- 1.248, 1.249,
1.236, 1.226, 1.216, 1.253 from 2 to 31 degrees -- and an absorption-depth
smear must grow with obliquity. This is the third measurement to say so, and
the cleanest.

**It is not the beam either.** A synchrotron beam is routinely wider in one
direction than the other, which would be fixed in the laboratory; resolved
along and across the rotation axis the ratio is 1.05, against 1.24 in the
reflection's own frame. Whatever it is prefers the radial direction, not a
laboratory one. Worth saying because the cubic crystal argues for isotropic
mosaicity, which is `sigma_M`, and says nothing about `sigma_D`, which is the
beam.

**Nor is it the crystal.** Diffraction from the front and the back of the
crystal starts from points separated along the beam, and at a scattering angle
those rays land `extent * tan(2 theta)` apart -- a hundred micrometres at thirty
degrees is fifty-eight, purely radial. The difficulty is that seen from the
crystal this goes as `sin(2 theta) cos(2 theta) / distance`, which is exactly
the sensor's dependence, so the two cannot be separated by this data; only
their sum can be tested, and only if the excess follows that shape.

It does not. Fitting a single source extent to `w2^2 - w1^2` gives an answer
that is not single:

    two-theta     excess      implied extent (T = sigma sqrt 12)
      2 - 11    8.27e-05        593 um
     11 - 15    1.51e-04        577 um
     15 - 17    1.56e-04        492 um
     17 - 19    1.64e-04        454 um
     19 - 20    1.54e-04        406 um
     20 - 22    1.50e-04        375 um
     22 - 24    1.47e-04        348 um
     24 - 31    1.59e-04        337 um

The excess is nearly flat while `sin(2 theta) cos(2 theta)` rises sixfold, so
the implied extent falls from 593 to 337 micrometres across the range. A real
source extent is a property of the crystal and cannot do that.

The magnitudes are worth having anyway. At the widest angle a 100 micrometre
crystal would contribute 1.4e-05 square degrees, the sensor alone predicts
2.5e-04, and 1.6e-04 is measured. So the crystal term is about a tenth of what
is seen and the sensor term alone already overshoots it.

**And it is not simply radial.** Split by azimuth around the beam centre the
radial ratio runs from 1.01 to 1.42, which a purely radial effect cannot do.
That may be confounded with radius, since the azimuth bins cover different
parts of the face, and it is not resolved.

The caveat from before still holds and still matters: the widths are 0.69 and
0.84 pixels. A second moment below the sampling interval is not a reliable
shape, and none of these numbers should be modelled until they can be
reproduced on data where a spot is several pixels across.

### Is the model any good? Ask the data, not DIALS

The Gaussian is supposed to hold essentially all of a spot's density by three
sigma. So measure it: over the flagged shoeboxes, what fraction of the counts
lies within n sigma, for n from one to four, in the detector directions and the
rotation direction separately.

    mxi_profile --compare --compare-sigma-b 0.031698 --compare-sigma-m 0.097667 \
                refined.expt refined.refl

On 1800 images of insulin, 70425 spots. A one-dimensional Gaussian holds
0.6827, 0.9545, 0.9973; two independent directions hold the square of that,
0.466, 0.911, 0.995.

                          1       2       3       4
    ours     detector  0.7578  0.9768  0.9990  1.0000
    ours     rotation  0.8653  0.9848  0.9988  0.9999
    DIALS    detector  0.7726  0.9801  0.9993  1.0000
    DIALS    rotation  0.8017  0.9668  0.9952  0.9995

**Both models hold by three sigma**, which is the first thing the test was for
and both pass: 0.999 against an expected 0.995 on the detector, 0.999 and 0.995
in rotation.

**The rotation direction says DIALS is right and we are not.** A correct sigma
would hold 0.6827 at one sigma. DIALS holds 0.8017 and this holds 0.8653 -- both
too concentrated, meaning both sigmas are larger than a Gaussian fitted to the
core, and ours is the further out. That is the same +21 per cent seen against
DIALS' number, now confirmed against the data rather than against DIALS.

**The detector direction says ours is marginally better**, 0.7578 against
0.7726 where 0.466 is expected, which matches being 2.9 per cent smaller.

Two things this does not say. The expected fractions assume the spot really is
Gaussian, and 0.76 where 0.47 is expected is far too large to be a small
sigma error: the real profile is much more peaked than a Gaussian, which is
why profile fitting uses learned reference profiles rather than the analytic
form. And the shoeboxes are the spot finder's, cut at its threshold, so the
density outside the box is not in the denominator at all and every fraction
here is an overestimate. Neither affects the comparison between the two sets of
sigmas, which is over identical pixels.

### Looking at the spots

    mxi_grid  refined.expt refined.refl --out grids.txt --n 5 --neighbours 300
    python3   docs/plot_grid.py grids.txt kabsch.png

Each strong spot's density on a grid in `(eps1, eps2, eps3)`, beside the
average of the spots nearest it on the detector, beside the average of all of
them. A single spot is a handful of pixels across and its grid is coarse and
noisy; the reference is what says whether a feature belongs to the spot or to
the sampling.

Two departures from Kabsch section 3.3, both because this measures the data
rather than applies the model. A pixel is subdivided in the detector plane and
its counts shared between the subdivisions, five ways per axis as Kabsch does,
because the data are badly undersampled and assigning a whole pixel to one grid
point turns the result into a staircase -- `--subdivisions 1` shows that. And
along `e3` an image's counts are shared between grid planes by the plain
geometric overlap of its angular range with theirs, NOT by the Gaussian weights
of Kabsch's `f_3j`: placing the counts with the model would beg the question
the picture is asked to answer.

The picture as a number. The average profile falls to a tenth of its peak at

    eps1   1.20 sigma_b        a Gaussian falls to a tenth at 2.146 sigma
    eps2   1.80 sigma_b
    eps3   1.20 sigma_m

So the spots are far narrower than the model that is supposed to describe
them, in every direction, and `eps3` worst -- which is the same conclusion the
captured-fraction test reached and the same one the comparison with DIALS
reached, now visible.

`eps1` and `eps2` differ from each other by half again, which a single
`sigma_D` cannot express: the model is isotropic on the detector face and the
spots are not.

### Is the anisotropy the sensor? Apparently not, and everything here is a pixel wide

    mxi_grid refined.expt refined.refl --out grids.txt --map map.txt
    python3 docs/plot_anisotropy.py map.txt anisotropy.png 0.0253

`e2` lies in the scattering plane, radially on the detector, and `e1` across
it. A photon absorbed at a random depth in the sensor is recorded further out
than where its ray entered, so that smear is radial: it belongs to `eps2`
alone. That makes it testable -- the excess `w2^2 - w1^2` should equal the
predicted smear squared and grow with obliquity as absorption says.

It does not. Over 61047 spots with more than fifty counts, binned by obliquity:

    obliquity     w2^2 - w1^2     predicted^2     ratio
     2 - 11 deg     8.3e-05         4.2e-05        1.98
    14 - 17         1.6e-04         1.1e-04        1.36
    19 - 20         1.5e-04         1.7e-04        0.93
    24 - 31         1.6e-04         2.5e-04        0.64

The measured excess is flat at about 1.5e-4 while the prediction rises eight
fold. At low obliquity there is twice as much anisotropy as the sensor can
account for, and at high obliquity two thirds as much. **Whatever it is, it is
not the depth of absorption**: that effect is real and is in there, but it has
the wrong shape.

The map of `w2/w1` across the detector face says the same thing more plainly:
it is not radial. A sensor effect on a flat detector must be, since it depends
only on the angle at which the ray strikes. This has large smooth patches that
do not centre on the beam.

**And the measurement sits at the sampling limit.** One pixel subtends 0.0253
degrees at 170 mm. The measured second moments are

    width1, tangential   median 0.69 pixels
    width2, radial       median 0.84 pixels
    sigma_b                     1.22 pixels
    a single lit pixel          0.29 pixels

So a spot is about two pixels across and its second moment is within a factor
of three of what a single lit pixel would give on its own. At that scale a
second moment is not a shape: it depends on where the spot centre falls within
a pixel, and the shoebox mask truncates it. An apparent anisotropy correlated
with position could be produced by that alone, and the non-radial map is as
consistent with a sampling artefact as with anything physical.

That caveat applies backwards as well. The earlier statement that `eps1` and
`eps2` differ by half again, and that the profile falls to a tenth of its peak
at 1.2 sigma, are both measurements made below the pixel scale. They are
evidence that the model is too wide, which is corroborated independently by the
captured fractions and by DIALS' own number; they are NOT evidence about the
shape of the underlying spot, which this data cannot resolve.

Two things worth doing before modelling any of it: repeat this on data with
finer sampling, where the spot is several pixels across; and test the sampling
hypothesis directly by measuring the width against the sub-pixel position of
the spot centre, which should show nothing if the widths are real.

## Seeing the model on the images

    mxi_mask refined.expt refined.refl -o masked.refl \
             --d-min 1.6 --first-image 0 --last-image 20
    dials.image_viewer refined.expt masked.refl

Predicts the reflections, works out each one's integration region, and writes a
shoebox whose mask marks it: `Valid | Foreground` inside the n-sigma region,
`Valid | Background` around it. **The pixel values are left at zero.** Nothing
here has read an image, and invented counts would be worse than none -- the
point is to see where the model says the signal is, drawn over the real image
by the viewer.

The bounding box comes from inverting the region rather than guessing a size:
the four corners of the `(eps1, eps2)` square are offset from `s1` in the
frame's own tangent plane and intersected with the panel, and the rotation
half-width is `n sigma_M / |zeta|` turned into images. It is then clipped to
the panel and the scan, so a box hanging off an edge is kept as the part of it
that is real.

An image range is not optional in practice. A whole 1800-image sweep at 1.6
Angstrom is 129000 reflections and about a gigabyte of empty shoebox; the first
attempt at this was killed by the machine. The boxes are also built and encoded
one at a time rather than all held at once.

### The region is a box, and the model is an ellipsoid

`|eps1| <= n sigma_D` and `|eps2| <= n sigma_D` and `|eps3| <= n sigma_M` is
what Kabsch section 3.1 writes for the mask, and what DIALS uses, so it is the
default. But it is a BOX in Kabsch space, and the model is a three-dimensional
Gaussian, whose surface of constant density is the ELLIPSOID

    (eps1/sigma_D)^2 + (eps2/sigma_D)^2 + (eps3/sigma_M)^2 <= n^2

These are not the same set and the difference is not small. At a corner of the
box all three coordinates are at n sigma at once, so the Gaussian is at
exp(-3n^2/2) there: 1.4e-6 of its peak at n = 3. There are eight such corners
and they hold nothing at all.

    volume          ellipsoid is pi/6 of the box, 0.5236
    density held    0.9919 in the box, 0.9707 in the ellipsoid

So the box buys two per cent more of the density with ninety per cent more
volume, all of it in places the model says are empty. `--shape ellipsoid`
draws the other one; on insulin it marks 23.9 per cent of the voxels where the
box marks 43.7, a ratio of 0.546 against the 0.5236 expected, the difference
being that a voxel is in or out as a whole and the boxes are only a dozen
pixels across.

Mapped onto the image the ellipsoid is not an ellipse: the transform from
Kabsch space to pixels and images is not a similarity, and it varies across the
detector with obliquity and with zeta. That distortion is the thing worth
looking at, and it is what the box hides.

What it shows on insulin, twenty images at 1.6 Angstrom: 1445 boxes, typically
12 by 12 pixels by 8 images. Those are big boxes for spots that are two pixels
across, which is the same conclusion as everything else in this section, in the
form most likely to be believed.

## Order of work

1. `sigma_D` and `sigma_M` from the indexed strong spots, using the shoebox
   pixels that are already in the file. Directly checkable: `dials.integrate`
   writes both into the `profile` block of its output experiment list, so a
   DIALS-written `integrated.expt` is the oracle.
2. Bounding boxes from the profile model, and the foreground mask. Checkable
   against the `bbox` column DIALS already writes.
3. The Kabsch coordinate transform and the polygon clipping, against synthetic
   shoeboxes with known overlaps.
4. Background, summation, profile fitting. Needs pixels to validate.

Step 1 has an oracle and no dependencies, which is why it is first.

## The precision contract, decided in advance

Settled per stage before any of it is written, not discovered afterwards:

* **Summation integration sums integer counts.** Done in integers, exactly, as
  the spot finder's window sums already are. There is no precision question and
  there should not be one.
* **The coordinate transform and the clipping** are geometry in double, and
  port to float with the same argument as the refinement target: measure the
  analytical path, do not assume it.
* **Background and profile fitting** are least-squares, and get a tolerance
  rather than a byte-for-byte test.

## How it will be compared

Both modes, from the start:

* **Pinned** -- this stage gets DIALS' upstream output. Isolates it. Without
  this, a hundredth of a degree of orientation difference shifts every shoebox
  and a correct integrator looks broken.
* **Cascade** -- this stage gets our own upstream output. The number that
  matters.

The join key at this boundary is (`miller_index`, `entering`, frame), and the
metric is the pull distribution, `dI / sqrt(sigma_a^2 + sigma_b^2)` -- its mean,
width and tails -- together with the correlation. Not percentage differences on
intensities.
