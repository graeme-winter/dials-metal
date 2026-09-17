# Integration

Not started. This is the plan, what it is bound by, and what is missing before
parts of it can be written at all.

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

## What is missing

**The papers with the formulas in them.** Winter et al. (2018) defers to others
for everything numerical, and none of them are here:

| | for |
| --- | --- |
| Kabsch (2010a), Acta Cryst. D66, 133–144 | `sigma_D` and `sigma_M` estimation; profile-fitted intensity and its error |
| Kabsch (1988b), J. Appl. Cryst. 21, 916–924 | the local reciprocal-space coordinate system |
| Parkhurst et al. (2016) | the Poisson GLM background |
| Leslie (1999) | summation error estimates |

Kabsch (2010a) is the one that blocks the most. Writing any of it from memory
would be exactly the kind of convention guessed rather than read that has cost
this project five separate format bugs.

**Pixels.** This is the first stage that needs them, and this container has
none: no image data, no HDF5, and therefore no way to build the spot finder or
read NXmx. Everything up to here worked on reflection tables alone.

That splits the work in two, and the split is useful rather than merely
unfortunate:

* **Needs no pixels** -- profile parameter estimation from indexed strong
  spots, bounding boxes, the foreground and background masks, the coordinate
  transform, the polygon clipping. All testable here against synthetic
  shoeboxes and known geometry.
* **Needs pixels** -- shoebox extraction, background fitting, summation,
  reference profiles, profile fitting. Developed against synthetic shoeboxes
  here; only ever validated by running against real images elsewhere.

## Order of work

1. `sigma_D` and `sigma_M` from the indexed strong spots. No pixels, and
   directly checkable: `dials.integrate` writes them into the `profile` block
   of its output experiment list, so a DIALS-written `integrated.expt` is the
   oracle.
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
