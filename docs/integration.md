# Integration

**Read `## Order of work, revised` first.** Several sections above it reach
conclusions that later work overturned; they are marked as superseded and kept
because the reasoning is instructive.

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

## Order of work, revised

Steps 1 to 3 of the original plan are done and what they found changed the
rest. What follows is written knowing that.

### What the profile work settled

**The Gaussian model is nearly right, and every measurement that said otherwise
was measuring the pixel grid.** Spots here are two pixels across. Comparing
their moments with a Gaussian in closed form said the model was wrong by
factors and that the spots were anisotropic; integrating the same model over
the same pixels, inside the same mask, put the positions within a hundredth of
a pixel and the widths within four per cent, and predicted the anisotropy from
projection and absorption alone.

That is the thing to carry into integration, because integration is nothing but
comparing a model with pixels. **Anything compared with the data has to be
rendered into the data's own terms first.**

**`sigma_D` and `sigma_M` do not match DIALS** and `mxi_profile` says by how
much every time it runs. The forward model says the truth is a few per cent
below both, so neither is exactly right and the difference between them matters
less than it looked. Integration takes the model as a parameter rather than
inheriting a number.

**The region is a box in Kabsch space and the model is an ellipsoid.** The box
is what Kabsch section 3.1 writes and what DIALS uses, so it is the default,
but its eight corners hold 1.4e-6 of the peak and the ellipsoid is pi/6 of its
volume. Once background is summed from that volume the difference stops being
cosmetic.

**The pixels were in the file all along** -- the `Shoebox<>` column -- and
`src/shoebox.h` reads them.

### Checked against a DIALS integrated.refl

`integrated.refl` has no shoeboxes, DIALS discards them, but it has `bbox`,
`num_pixels.foreground`, `num_pixels.background` and `num_pixels.valid`, which
is enough to check the masking with no pixels at all.

**Predictions agree exactly.** Over 4018 reflections on the first sixty images:
100 per cent within half an image in `z`, 100 per cent on the `entering` flag,
and the box `z` extents agree 99.3 per cent of the time.

That only worked once it was compared against the right file. Predicting from
`refined.expt` and joining to `integrated.refl` matched almost nothing, because
something between them reindexes: `real_space_a` is (20.46, -14.89, 62.45) in
one and (-20.50, 14.81, -62.46) in the other, so the Miller indices are not
comparable across that step. It looked like a prediction bug and was a join
bug. **Pin the geometry to the file being compared against.**

**The box was wrong, and wrong in a way that blocks integration.** DIALS':

    num_pixels: foreground 462, background 3094, valid 3542
    foreground + background = valid   for 100 per cent of rows
    valid / box volume = 1.000

It is a measurement box in Leslie's sense: the foreground region plus a rim,
the rim being 87 per cent of the volume. `mxi_mask` built the foreground alone,
which has no background in it, and the background estimate has to come from
somewhere. `box_scale` widens the box on the detector; at 1.9 the extents match
DIALS. The rotation direction is untouched, since those extents already agreed.

**One thing still disagrees.** Our foreground is 372 voxels where DIALS' is
462. The `z` extents match, so it is on the detector, and 1.11 per axis is the
wrong size to be rounding. Recorded rather than chased, with the oracle that
will settle it.

### Real images, and what the background looks like

A 30 degree sweep of insulin on a 4M detector, with its own `refined.*` and
`integrated.*`: 300 frames of 2162 by 2068, bitshuffle compressed, 297 MB.
`h5py` with `hdf5plugin` reads it here; one frame decompresses in 0.02 s.

**65535 is the bad-pixel marker**, not a count. It is 5.8 per cent of the frame
-- module gaps and dead pixels -- and reading it as data puts 17 billion counts
on a frame that has 1.65 million.

**The pixel addressing is `data[frame, slow, fast]`**, verified rather than
assumed: over the forty strongest reflections, the centroid of the extracted
box against `xyzcal.px` is

    fast   median -0.046  rms 0.146 px
    slow   median -0.145  rms 0.303 px
    frame  median -0.150  rms 1.111 images

Swapping fast and slow would have shown as a large offset in both.

**The background is 0.3 counts per pixel**, which is the regime Parkhurst's
paper is about: most pixels are 0 or 1, and a normal approximation to a Poisson
that small is not one. On 1482 reflections, taking the outer shell of each box
as background so the region is not chosen by value, against DIALS' robust GLM:

    DIALS (robust GLM)   0.3044 counts/pixel
    plain mean           0.2993    -1.3%
    5 per cent truncated 0.2334   -22.5%
    3 sigma clipped      0.2343   -14.2%
    median               0.0000  -100.0%

Which is Parkhurst's Table 1 and Figure 3 reproduced on this data: every
traditional outlier rejection biases the background DOWN, the median collapses
to zero because most pixels are zero, and the GLM sits beside the unrejected
mean. Since the intensity is the foreground minus the background, a background
biased low is an intensity biased high, for every reflection in the dataset.

A first attempt at this measurement used "the lowest n pixels" as the
background region and produced a 48 per cent bias for the plain mean. That was
the sampling, not the estimator: selecting the low tail guarantees a low
answer, and the median coming out at exactly zero should have been the giveaway
before the numbers were read.

### Background: implemented and checked against DIALS

`src/background.h`, Parkhurst's robust GLM with a Poisson link and Huber
weights at c = 1.345. Against DIALS' own `background.mean` on 2000 reflections
of the real 30 degree sweep:

    DIALS       median 0.2997 counts/pixel
    ours        median 0.2963
    difference  median -1.35 per cent, rms 3.8, all within 20
    correlation 0.9972, converged 100 per cent, median 4 iterations

**The remaining difference is the region, not the estimator**, which was
checked rather than assumed. Widening the background region towards the number
of pixels DIALS uses moves it monotonically towards zero:

    pixels used   1120    1596    2016    2380     (DIALS 2936)
    difference   -1.46%  -1.40%  -0.89%  -0.54%
    correlation  0.9946  0.9973  0.9986  0.9991

So the estimator agrees and what is left is that the background region here is
a shell of the box rather than everything outside the foreground mask. That is
`mxi_mask`'s job and is the next thing to connect.

**The expectations are summed, not solved.** Appendix B gives closed forms for
C1 and C2 in regularized gamma functions, to avoid summing over the
distribution. At 0.3 counts a pixel the distribution is over by twenty, so the
sum is a dozen terms and is exact; the closed forms are an optimisation for a
regime this is not in, and they are the part of the paper whose exponents did
not survive the PDF's text layer. Transcribing them from a guess would have put
an error where no test could find it.

The tests check the reason rather than the answer: with the tuning constant
large the estimator must be the plain mean and the corrections must vanish to
their known values, C1 to zero and C2 to sqrt(mu); C1 and C2 are compared with
an independent summation using factorials rather than the production
recurrence; the fit is unbiased on clean Poisson data at four background
levels; one outlier in five hundred moves the mean by a factor of five and the
GLM by less than a tenth; and an all-zero background is answered rather than
iterated towards, since there is no log of zero.

### Summation integration, on real images

    mxi_integrate integrated.expt master.nxs -o mine.refl \
                  --sigma-b B --sigma-m M --d-min 1.8 --save-shoeboxes

Predicts, builds each measurement box, fills it from the images through the
spot finder's NXmx reader, fits the background and sums the foreground. The
output is a DIALS reflection table, so `intensity.sum.value` compares directly.
On 2504 reflections of the real 30 degree sweep:

    background.mean      ours 0.3042   DIALS 0.3019   correlation 0.9991
    intensity.sum.value  ours  155.0   DIALS  160.4   correlation 0.8753
    intensity.sum.var    ours  310.1   DIALS  348.6   correlation 0.8730

    on the 1112 with I/sigma > 10 in DIALS:
      ratio ours/DIALS   median 0.9817, 10th 0.9410, 90th 1.0038

**The background is right and the intensity is two per cent low**, which is
what the foreground mask being small predicts: 372 voxels against DIALS' 462,
measured before any of this touched a pixel. A foreground that misses the
skirts of a spot misses signal, and a two per cent deficit on strong
reflections is the size of that difference. So the mask discrepancy is no
longer a curiosity to be chased when convenient -- it is the outstanding error.

**The profile model is estimated rather than given.** In order of preference:
the command line, then this package's own estimate from a `.refl` of strong
spots with shoeboxes, then the `profile` block of the `.expt` if it carries one
-- which is now read. Carrying two numbers by hand between two programs is a
way to run the second on the wrong ones.

Which exposes something worth not mistaking for progress:

    DIALS' sigmas   sigma_b 0.028676  sigma_m 0.092108   ratio to DIALS 0.9817
    our estimates   sigma_b 0.027315  sigma_m 0.128501   ratio to DIALS 0.9943

Our own estimates give intensities closer to DIALS' than DIALS' own numbers do.
That is compensation, not correctness: our `sigma_m` is forty per cent larger,
which widens the foreground in the rotation direction and offsets the mask
being too small on the detector. Two errors of opposite sign look like
agreement, and the only reason this is visible at all is that both were
measured separately first.

**The images come from the `.expt`.** Its imageset block already says where
they are; being told the path a second time on the command line is how the two
come to disagree. `--images` overrides it for data that have moved since the
file was written, and the program prints which it used.

That found a real bug in the NXmx reader, reproduced before it was fixed. When
the virtual dataset maps back into the master file itself -- HDF5 writes `"."`
for that -- the reader set the path to the master and then resolved it against
the master's own directory a second time. A master named `../ins10_1.nxs`
became `../../ins10_1.nxs` and failed to open.

The first fixture written for it did not reproduce the failure: with the master
at `../sub/self.h5` the doubled resolution gives `../sub/../sub/self.h5`, which
is the same file. It only escapes when the master is one level up, which is
exactly the case reported. A fixture that does not fail on the old code is not
a test of anything.

**Every prediction that gets no shoebox is counted, by reason.** On a ten
rotation sweep 749786 predictions became 367050 boxes and there was no way to
tell whether that was the detector, the rotation axis, or a bug. The reasons
are now printed:

    21521 reflections predicted
    21032 shoeboxes to fill
      471 a corner misses the detector
      18 zeta below the cut

The reasons are: no such panel, no Kabsch frame, zeta below the cut, a corner
of the region missing the detector plane, a box spanning more images than any
reflection should, and clipped away entirely by the panel or the scan. Each is
reachable and each has a test, because a reason that is never produced is a
string nobody has checked.

**A multi-turn sweep was predicted once, not once per turn.** The reason
counters found it on their first outing: 749786 predictions of a ten rotation
sweep, of which 378278 had a `z` off the end of the scan, against DIALS'
7516507. The ratio is 10.02, which is the number of turns.

Two faults, one line apart:

    if (span < 2.0 * kPi) { phi = wrap_from(phi, lo); ... }

A scan of a full turn or more skipped the wrap entirely, on the grounds that
everything is inside it. Everything is -- but only after wrapping, and the
caller needs phi in the scan's own coordinates or `z_from_phi` puts it in the
wrong turn. That is the 378278.

And a reciprocal lattice point that crosses the Ewald sphere at `phi` crosses
it again at `phi + 2 pi`, so ten rotations record every reflection ten times
over -- which is the entire reason for collecting them. Each turn's root is
converged from its own seed rather than taken as `phi + 2 pi k`, because with a
scan-varying crystal the setting matrix at turn nine is not the one at turn
zero.

The patch went into `predict_indices` first and changed nothing, because
`predict` has its own copy of the loop and is the one that runs. There is now
one emitter that both call.

**It walks the frames once, with shoeboxes open across them.** Every shoebox
cannot exist at once: ten rotations of insulin is 7.4 million reflections and,
at about 3500 voxels each, 233 GB of pixels.

Blocks of images were the first attempt and were wrong. A reflection near the
rotation axis spans `n_sigma sigma_m / |zeta|` in phi, which at the zeta cut is
about a hundred and fifty frames; one of those in a block forces the whole
block to read that far, and the next block reads it again. On three hundred
frames it read 2951 of them at a block of ten, 817 at fifty -- which is how
this was noticed, a three hundred frame sweep reporting 817 frames read.

So the frames are walked once in order. A shoebox opens when the frame it
starts on comes round, takes a slice from every frame it spans, and is
integrated and released on the frame it ends. Every frame is read exactly once,
and memory is bounded by what is open at the time -- which is what the block
size was trying to bound anyway, and bounds it by the data rather than by a
guess.

    frames read     300 of 300
    peak memory      86 MB
    open at once    931 shoeboxes

with the Miller indices, intensities, variances and backgrounds identical to
the block version.

The first pass is unchanged: work out every bounding box and keep nothing but
the box, 24 bytes a reflection rather than 31 kB.

**Where the time goes**, from `--timing`, on the thirty degree sweep:

    the profile model             0.171 s    4.4%
    prediction                    0.251 s    6.4%
    bounding boxes                0.010 s    0.2%
    opening shoeboxes             0.876 s   22.5%
    fetching frames               0.281 s    7.2%
    decompressing                 0.835 s   21.5%
    filling shoeboxes             0.312 s    8.0%
    background and summation      1.049 s   26.9%
    writing                       0.038 s    1.0%
    total                         3.894 s

It was 21 seconds when the timing went in, and finding out where meant adding
the phases twice: the first breakdown accounted for 15 per cent of the run and
the second for 35, and each gap was somewhere the clock had not been put. A
breakdown that does not add up is not a breakdown.

Two things it found, neither of which was the images.

**Every chunk was fetched twice.** A map from frame number to key was built by
reading every key to ask what frame it was, and then the loop read them all
again -- 18 of the 21 seconds. The loop is now driven by the frames as they
arrive, and checks they arrive in order rather than assuming it.

**Building the masks was 63 per cent of what was left.** `build_shoebox`
evaluated the Kabsch mapping once per voxel, and `eps1` and `eps2` depend only
on where a pixel is on the detector while `eps3` depends only on which image it
is. So the mapping was being run `nz` times for every pixel, and it has a
parallax correction with an `exp()` in it. Computing the face once and the
images once is `nx*ny + nz` evaluations instead of `nx*ny*nz`: 5.344 s to
0.876, with every intensity, variance, background and foreground count
identical.

    21.0 s -> 8.5 s -> 3.9 s

### Threading: the reads were the wrong thing to thread

On ten rotations, 310 seconds:

    decompressing               158.107 s   51.0%
    background and summation     53.300 s   17.2%
    opening shoeboxes            46.137 s   14.9%
    prediction                   25.739 s    8.3%
    filling shoeboxes            11.755 s    3.8%

Threading only the fetch and decompress settles at about 116 per cent of one
core on a sixteen core machine, and the arithmetic says why: 46 + 12 + 53 is
111 seconds of work on the thread that owns the shoeboxes, so the readers
finish their lookahead and wait. No number of reader threads moves that.

**The ordering that forced the pipeline was not needed.** A frame writes only
its own z plane of a shoebox, so two frames of the same box can be filled by
two threads without touching the same double. Nothing needs ordering; the boxes
need only to exist.

So the scan is cut into windows. Within one: the masks are built in parallel,
the frames are fetched, decompressed and filled in parallel, and the boxes are
integrated in parallel. Only the window boundary is serial.

**The window must be long compared with a shoebox**, which the first two
attempts were not. Bounded at a thousand frames it swallowed a three hundred
frame sweep whole -- 21032 boxes, 650 MB, and the allocation cost more than the
parallelism saved. Bounded at four thousand boxes instead it became six windows
of fifty frames and read 748 frames instead of 321, with decompression going
from 1.13 seconds to 2.66: the overlapping-block problem again, in a new place.
Frames are the bound and boxes are a safety net, at twenty thousand.

**What this costs and what it cannot show.** This container has one core, where
the window version is 5.2 seconds against the streaming version's 4.5 and peaks
at 843 MB against 86. Both are real and neither is offset here; the case for it
is that three phases now parallelise instead of one, and that is an argument
rather than a measurement. The results are identical at every thread count, and
the numbers to watch elsewhere are the total, the CPU share, and whether frames
read stays near the frame count -- if it climbs, the window is too short.

`--threads`, `--window` and `--max-boxes` are the knobs.

### Every turn collapsed into the first, and a diagnostic found it

    27366 frames read (3627 wanted by a shoebox, each read 7.55 times)

3627 frames of a 36000 image sweep is one rotation. Every shoebox of all ten
turns was landing in the first turn's frames, so the same pixels were
integrated ten times over and every intensity of that run was wrong.

Two faults, and the second hid behind the first.

`wrap_from` returns a value in `[lo, lo + 2 pi)`, and `build` called it
unconditionally -- so a phi correctly placed in turn nine was wrapped straight
back to turn zero. It now wraps only a phi that is OUTSIDE the scan; one
already inside it is already in the turn its caller meant.

And `converge_root` answers with a crossing from `ewald_intersections`, which
is in a principal 2 pi interval whatever it was seeded with, so the turn was
lost on the way out too. `emit_turns` now puts it back, taking the turn nearest
the seed.

**The test that should have caught this used a scan-static crystal**, where the
turn was restored by hand with `c.phi = seed` -- so it passed, while every real
experiment, all of which have scan-varying crystals, was broken. There is now a
test on a scan-varying crystal with the same `A` at every scan point, so the
geometry is identical to the static case and only the code path differs; it
fails on the old predictor and asserts that every turn holds its own share of
the predictions rather than only that the count is right. A count can be right
while every one of them is in the wrong place.

The diagnostic that found it exists because `27366 frames read` was not
interpretable on its own, which was itself a complaint about the output rather
than a hypothesis about the code.

### The re-read rate, and where it comes from

With the turns no longer collapsed, ten rotations wants all 36000 frames and
reads 59134 of them, 1.64 times each. That is arithmetic rather than a puzzle:

    7309317 boxes / 20000 per window = 365 windows
      own span       99 frames
      actually read 162 frames
      overlap        63 frames past each window end

and a shoebox at the zeta cut of 0.05 spans 82 images, because its extent in
phi is `n_sigma sigma_m / |zeta|`. So the overlap IS the near-axis reflections,
and the lever is window length against memory:

    --max-boxes  20000:   99 frames per window, re-read 1.65 x, 0.6 GB
    --max-boxes  40000:  197 frames per window, re-read 1.32 x, 1.3 GB
    --max-boxes  60000:  296 frames per window, re-read 1.22 x, 1.9 GB
    --max-boxes 120000:  591 frames per window, re-read 1.11 x, 3.8 GB

Raising `--min-zeta` shortens the long boxes instead and costs reflections.

### The timing report was summing to 124 per cent

Fetching 77.1 per cent and decompressing 47.4 per cent of the same run: both
were the busiest thread's wall clock inside one parallel region, so they
overlap and the percentages are not shares of anything. The region now reports
its own wall clock as the phase, with the work inside it in thread-seconds and
a ratio against that wall -- which are additive, comparable, and show the
difference between working and waiting.

### What that left, on sixteen cores

    total                        65.403 s   from 310
    user                        370       s   so 5.7 cores busy

    prediction                   25.982 s   39.7%
    decompressing                22.716 s   34.7%
    writing                       8.203 s   12.5%
    background and summation      5.461 s    8.4%
    opening shoeboxes             4.792 s    7.3%
    the profile model             4.541 s    6.9%
    filling shoeboxes             1.825 s    2.8%

Prediction became the largest item by being the only phase still serial, so it
is threaded too: one unit of work per `h`, each collecting into its own vector,
joined in `h` order. The order matters as much as the content -- everything
downstream is indexed by position in that list, so a prediction list that
reordered itself with the thread count would make every comparison between two
runs meaningless. Checked with a checksum that depends on position: identical
at one thread and four.

The report also says how many frames a shoebox wanted and how many times each
was read, because `27366 frames read` of a 36000 image sweep is not
interpretable on its own -- it could be windows re-reading, or frames nothing
needed, or reads that failed. Reads that find no frame are counted separately:
an unallocated chunk is a frame the writer never received, and a shoebox
spanning one is missing a slice and will integrate low.

Three things the pixels settled that no synthetic test would have:Three things the pixels settled that no synthetic test would have:

* **The bad-pixel marker is excluded from both sums, not counted as zero.** It
  is 5.8 per cent of a frame, module gaps and dead pixels, and 470000 voxels of
  the boxes here. Counting them as zero would drag the background down wherever
  a gap crosses a shoebox.
* **`d` is `1 / |A h|` from the STATIC cell.** Three nearly equal numbers were
  candidates: the static cell agrees with the column to 1.2e-15 for every
  reflection, while `1/|s1 - s0|` and the scan-varying `A` are both 2.6e-4 out.
  A tolerance of a part in a thousand would have accepted any of them.
* **`partiality` is the fraction of the rocking curve inside the BOX**, not
  inside the scan. What was summed is what is in the box, and a box spanning
  three sigma holds 0.9973 of the curve rather than all of it. Against DIALS:
  0.99887 against 0.99914, with ninety per cent of differences below 0.0014,
  and 208 reflections that both call genuinely partial. Computing it against
  the scan instead gives 1.0000 for everything away from the ends, which is the
  kind of small error a scale factor multiplies through a whole dataset.
* **`zeta`, `partial_id`, `num_pixels.background_used` and
  `xyzobs.mm.variance`** are written too. What is still missing is the three
  profile-fitting columns, which is the next step rather than an oversight.
* **`xyzobs.px.value` is the centre of mass of the foreground**, background
  subtracted, which `dials.scale` requires. Over the whole box instead it is
  0.15 to 0.19 pixels from DIALS'; over the foreground it is 0.07 in fast and
  0.10 in slow, and `xyzobs - xyzcal` comes out at [-0.001, -0.007, -0.076]
  against DIALS' [-0.002, -0.002, -0.076] -- including that frame offset, which
  is the part that would show a convention error.

  `xyzobs.px.variance` does **not** reproduce DIALS'. The second moment alone is
  about twenty times too large and the variance of the mean about thirty times
  too small, and neither is a constant factor away, so DIALS is computing
  something else; its values sit near 0.1 square pixels on every axis, which
  looks like a quantisation term rather than anything that scales with
  intensity. What is written is the variance of the mean, which is the standard
  quantity, and this is recorded rather than tuned to match.
* **`lp` is `L / P`**, with `L = |s1 . (m2 x s0)| / (|s1||s0|)` and
  `P = (1-p) + (2p-1)(u.n)^2 + p(u.s0hat)^2`. Recovered from the oracle rather
  than recalled: the Lorentz part was identifiable as the candidate with the
  least scatter against the column, and the three coefficients of P were then
  solved for by least squares and came back as 0.001000, 0.998000 and 0.999000
  with a residual of 2e-16 -- which is `(1-p)`, `(2p-1)` and `p` for the
  `p = 0.999` in the file, not three free numbers. At `p = 0.5` it collapses to
  the unpolarized `(1 + cos^2 2theta)/2`, which is the check that the form is
  physics; there is a test for that limit. Against the column it agrees to
  4e-5, the residual being our predicted `s1` rather than the formula.
* **`qe` is `1 - exp(-mu t / cos theta)`**, the fraction of photons the sensor
  stops rather than passes. Checked against the `qe` column: identical for 100
  per cent of reflections to 2e-16. It is stored rather than applied, as DIALS
  does, so the scaler applies it later and the two tables stay comparable.
* **The flag bits were read from the oracle, not guessed**: 256 is
  integrated_sum and 512 integrated_prf, from flag value 769 splitting 21972
  and 20688 ways.

`--save-shoeboxes` keeps the pixels, the mask and the fitted background in the
output, at about a megabyte per hundred reflections. It is the difference
between "the intensity is wrong" and "the intensity is wrong because the mask
is here and the spot is there".

### Profile fitting: implemented, and not yet right

`src/reference.h`, Kabsch sections 3.3 and 3.4 with Leslie section 6. Two
passes over the images, which is the slow way round and the right one: the
reference profiles are learned from the reflections themselves, and a profile
learned from part of a scan and applied to the rest would be a different
algorithm whose errors would be hard to attribute.

Both learning and fitting happen on a cube in the Kabsch frame. That is what
the frame is for: a reflection's shape is the same there wherever it sits on
the detector, so profiles from different reflections can be averaged, which on
the detector they cannot. A pixel maps to a region rather than a point, so each
is subdivided five ways an axis as Kabsch does, and one image covers a range of
`eps3` which is shared between grid planes in proportion to the overlap.

**What works.** Against a DIALS `integrated.refl` on 19360 profile-fitted
reflections:

    intensity.prf.value   ours 153.1   DIALS 159.7   ratio 0.9813 on I/sig > 10
    profile.correlation   ours 0.8141  DIALS 0.8709

and the unit tests hold both of Leslie's limits: fitting a profile to itself
returns the intensity exactly, and the fit reduces to the sum when the
background is negligible, which is section 6.4 and is a check on the weighting
rather than a remark.

**What does not.** The per-reflection correlation with DIALS is 0.4414, against
0.9453 for the summed intensities from the same run. So the scale is right and
the weighting is not. And the variance is too small: on weak reflections it is
52.7 against the summed 185.4, a factor of 3.5 where Leslie section 6.6
predicts about 2 for a typical profile.

Both point the same way. The grid is 729 points and a single pixel is spread
over several of them, by the subdivision on the detector and by the `eps3`
overlap across planes. Neighbouring grid points are therefore correlated, and
`1 / sum(P^2 / v)` treats them as independent observations -- which overcounts
the information, understates the variance, and mis-weights the fit. The pixels
are the independent measurements, not the grid points.

So the next step is to fit against the pixels with the profile transformed onto
them, rather than against the grid with the pixels transformed onto it -- or to
carry the covariance the transform induces. Recorded here rather than tuned
away, because a factor of 3.5 where the theory says 2 is the kind of
disagreement that a fudge factor would hide permanently.

The timing, on the thirty degree sweep at one thread: learning 7.8 s, fitting
20.1 s of 35.1 s total. It was 65 and 150 before the subdivided face was
computed once per box instead of once per image -- the same optimisation, and
the same mistake, as the mask.

### The steps

1. **Bounding boxes and masks.** Done, `mxi_mask`, and checked above.
2. **The Kabsch transform and the grid.** Done, `mxi_grid`.
3. **Background.** Parkhurst et al. (2016) is in hand: a robust GLM with a
   Poisson link, Huber weights at c = 1.345, and the constant-background case
   simplified in its Appendix B to a scalar iteration. The oracle is
   `background.mean` in a DIALS `integrated.refl`, which is present.
4. **Summation.** Leslie (1999) is in hand: the intensity and, more to the
   point, the variance, `G(Is + Ibg + (m/n) Ibg)`, which says the background
   dominates the error for weak reflections. The oracle is
   `intensity.sum.value` and `.variance`.
5. **Reference profiles and profile fitting.** Kabsch sections 3.3 and 3.4 and
   Leslie section 6. The profile work says why a learned profile is necessary
   rather than merely traditional: the real spot is more peaked than a
   Gaussian.
6. **Predicted shoeboxes need images.** Everything above is checkable on the
   strong spots' own pixels or against the oracle columns. Integration proper
   is not: predicted reflections are a superset of the strong ones.

### What this will cost, before it is written

Integration touches every pixel of every shoebox of every predicted reflection.
The performance work was practice for it and three lessons apply directly.

**Measure before optimising, and check the instrument first.** Two conclusions
here were drawn from timings that were wrong: a clock inside a lambda called
four hundred and fifty million times, and a benchmark passing the wrong kind of
memory to a device.

**Do not materialise what can be consumed.** Refinement builds a Jacobian of
parameters by reflections by three doubles, reallocated every iteration, and
threading its inner loop gave 1.29x on sixteen cores because the allocation is
serial. A shoebox does not have to exist as an array to be summed.

**The layout is expensive to change afterwards.** Decide whether a shoebox is
pixel-major or reflection-major, and whether the profile grid is dense or
sparse, before writing the loop that walks it.

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

## mxeq trend: where two integrations disagree, and against what

    mxeq trend ours.refl theirs.refl
    mxeq trend ours.refl theirs.refl --value intensity.prf.value --bins 5

A single correlation says two columns disagree. It does not say whether the
disagreement is with resolution, with intensity, with position on the detector,
with how near a reflection sits to the rotation axis, or with how well the
profile describes it -- and those point at different causes. This bins the
ratio against each in turn.

Bins hold equal populations rather than equal widths: equal widths on a
quantity like I/sigma put almost everything in the first bin and say nothing.

On the profile-fitted intensities it localised the problem straight away. The
ratio is flat -- 0.966 to 0.989 across every variable -- so the scale is right
everywhere. What moves is the agreement, and it moves the wrong way round:

    intensity.prf.value against I/sigma
              from         to        n       ours     theirs    ratio    corr
            -3.729      2.539     3846      11.97      12.32   0.9720  0.3718
             6.238      11.72     3845      150.5      155.7   0.9704  0.9008
             23.61      357.2     3846       1693       1783   0.9873  0.2288

    intensity.sum.value against I/sigma
            -3.729      2.465     3872      13.23      12.41   0.9593 -0.0099
             6.151      11.64     3872      154.6      158.6   0.9714  0.1163
             23.48      357.2     3872       1731       1735   0.9994  0.9751

Summation agrees with DIALS almost perfectly on the strongest reflections
(0.9751) and hardly at all within the weak bins, which is what noise looks
like and is unremarkable. Profile fitting does the opposite: it agrees best in
the middle and WORST on the strongest (0.2288), where there is the least noise
to blame.

So the fit misbehaves where the signal is largest, which is the opposite of
what a weighting error alone would do to weak data, and it is what drags the
overall correlation to 0.42. Strong reflections are where the profile's tails
carry the most counts and where a saturated or mismodelled pixel has the most
leverage. That is the next thing to look at, and this tool is how.

Note that within-bin correlation is not a measure of agreement on its own: a
narrow bin of a noisy quantity correlates poorly however well two programs
agree, which is why the summation rows read as they do. The ratio column is
what to read for agreement and the correlation column for outliers.
