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
| Parkhurst et al. (2016) | the Poisson GLM background |
| Leslie (1999) | summation error estimates |

**Kabsch (2010a) is here now**, and it carries most of what was missing:

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

## sigma_D: implemented, and 4.2 per cent from DIALS

To reproduce:

```sh
mxi_profile refined.expt refined.refl
```

The table has to still have its shoeboxes -- the estimate is made from the
pixels -- which `dials.find_spots` writes and which `mxi_index` and
`mxi_refine` now preserve. A stripped table is refused by name rather than
silently producing a number from nothing.

`src/profile_model.h`. Kabsch section 3.1, as written: for each strong spot,
the counts-weighted variance of the angles between its foreground pixels'
diffracted-beam directions and its own `s1`, background subtracted first; then
`sigma_D` is the root mean of those variances.

On the 1800-image insulin, over all 78618 reflections:

    ours              0.030355856 degrees
    dials.integrate   0.031697888 degrees      -4.23 per cent

**That difference is not explained and the two are not interchangeable until it
is.** What has been ruled out by measurement rather than argument:

    pixel centres at +0.0 instead of +0.5      +9.08 per cent   (wrong the other way)
    dividing by w instead of (w - 1)           -4.63
    variance about the centroid, not s1       -15.28
    parallax-corrected pixel directions       -15.34
    foreground mask instead of all valid       identical, since every valid
                                               pixel in these boxes is foreground

The parallax result is the interesting one. Undoing the correction moves the
spread 15 per cent the wrong way, which says DIALS is not doing it -- and makes
sense: the correction describes where a ray of a given direction is recorded,
so applying its inverse to a pixel asks where the ray came from, which is a
different question.

Four per cent is small enough to be tempting to chase by turning knobs until
the number matches. There are enough knobs in this recipe that a wrong
estimator could be tuned onto the right answer, so the tests plant a known
angular spread and check it comes back, rather than checking agreement with
DIALS. Settling the remainder needs DIALS' own source, which is not here.

Candidates worth testing when it is: which reflections DIALS selects (it may
exclude by `zeta`, by resolution, or by a minimum count); whether it weights by
counts or by counts minus background; and whether its mask has already been
narrowed from the spot finder's.

## sigma_M: implemented, and a factor of three from DIALS

    mxi_profile refined.expt refined.refl

Kabsch section 3.1: maximise the likelihood of the observed offsets under
`R(Delta, sigma_M/zeta)`, the fraction of a reflection's intensity recorded on
an image whose centre is `Delta` from its Bragg angle.

**One sample per image, not one per reflection.** With one per reflection the
offset is bounded by half an oscillation width by construction, the likelihood
is maximised by driving sigma to zero, and the estimate is meaningless -- which
is what the first attempt produced, 0.00006 degrees. A mosaic crystal puts a
spot on images well away from its Bragg angle, and that spread is the entire
signal.

Two rejections, both from the paper rather than invented:

* `|zeta| >= 0.05`. A reflection near the rotation axis has a modelled range of
  `sigma_M/|zeta|`, which diverges; such samples carry no information and would
  decide the answer. Measured on synthetic data: without the cut they pull the
  estimate DOWN by a third, not up, because they make every sigma look equally
  bad and flatten the likelihood.
* Kabsch step (vii), rejecting a spot whose observed centroid is far from its
  predicted angle. Without it the estimate is 0.506 rather than 0.293. The
  paper says "deviates too much" without a number; measured, the answer is flat
  at 0.293309 for any cut between one and ten images, so within that plateau it
  is not a knob.

On the 1800-image insulin, 283242 samples from 77153 spots:

    ours              0.293309 degrees
    dials.integrate   0.097667 degrees      a factor of 3.003

**Unexplained, like the 4.2 per cent on sigma_D.** The factor is close enough
to `n_sigma = 3` to be suspicious, and that is a coincidence worth testing
rather than a conclusion. What has been ruled out: it is not the zeta
correction, since dropping zeta entirely gives 0.309 rather than 0.098 and the
median zeta over the samples used is 0.64, not a third.

The estimator itself is checked against planted values rather than against
DIALS: samples drawn from the model with a known sigma come back within five
per cent at 0.05, 0.1 and 0.3 degrees, and the density is verified to integrate
to one over Delta, without which the product of these is not a likelihood at
all.

One trap on the way, the same one as twice before: rows with no prediction
carry uninitialised `xyzcal` -- denormals, not zeros -- and feeding them to the
likelihood gave 1.94 degrees. `has_prediction` now exists on the C++ side too.

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
