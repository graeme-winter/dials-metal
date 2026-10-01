# mxeq

Equivalence checks between two runs of an MX processing pipeline, boundary by
boundary, without cctbx -- and, for integration, tools that find and explain
where two runs differ.

```sh
mxeq check strong      dials/strong.refl     gpu/strong.refl
mxeq check indexed     dials/indexed.refl    gpu/indexed.refl -e dials/indexed.expt
mxeq check refined     dials/refined.expt    gpu/refined.expt
mxeq check integrated  dials/integrated.refl gpu/integrated.refl -e dials/refined.expt
mxeq check scaled      dials/scaled.refl     gpu/scaled.refl     -e dials/scaled.expt
```

The target is **equivalence, not identity**. Identity is achievable for
thresholding and stops being achievable somewhere around integration, and a
checker that demands it would report every difference in summation order as a
failure. So each boundary is compared on quantities that mean the same thing on
both sides, and the tool reports distributions rather than verdicts.

## Comparing two integrations

`check` says whether two runs are equivalent at a boundary. These say where and
why two integrations differ:

```sh
mxeq explain  ours.refl theirs.refl                 # why rows went unpartnered
mxeq trend    ours.refl theirs.refl --value intensity.prf.value
mxeq html     ours.refl theirs.refl -o comparison.html
mxeq disagree ours.refl theirs.refl -o disagree.refl --sigma 5
mxeq residuals integrated.refl -o residuals.html    # one file: predicted positions
mxeq profiles profiles.txt -o profiles.png          # mxi_integrate --save-profiles
```

**The ratio is always the first file over the second**, so the order of the
arguments sets which way a ratio above one points. Put the one being judged
first.

**Rows are paired as observations, not reflections.** Miller index and entering
flag group them and the frame separates them, because a sweep of several turns
records each reflection many times; `--radius` (5 images) is a sanity check and
not a discriminator. Every report states how many matched. **That count should
be close to the smaller table's rows**, and when it is not, nothing after it is
worth reading -- run `explain`, which says for each unpartnered observation
whether the other program did not predict it, put it in another turn, disagreed
about its side of the Ewald sphere, or placed it further than the radius.

Bins hold equal populations, since equal widths on a quantity like I/sigma put
nearly everything in the first bin. The band on a ratio is the robust spread,
not the standard error.

`trend` and `html` bin against everything a disagreement might follow:
resolution, I/sigma in the second file, |zeta|, partiality, and -- from the
first -- profile correlation, foreground pixels and background, then distance
from the beam centre and frame. `disagree` writes only the disagreeing rows, ours
with the other side's value alongside, for `dials.image_viewer`; its three
criteria find different things -- `--factor` proportional error, `--difference`
absolute counts, `--sigma` the only one that knows whether a difference is
larger than the measurement -- and each one given must be exceeded.

`residuals` reads `xyzres.px`, the observed centre less the predicted, and
shows it against image, resolution and detector position, with a summary that
separates the systematic offset, the counting noise and the prediction error.
See `docs/integration.md` for what it cannot tell you in z.

## Bias against the data itself

`mxeq equivalents scaled.expt scaled.refl` compares every observation with the
symmetry equivalents of its reflection, so that a bias in one kind of
observation shows without a second program or a reference data set. The
intensities are corrected as scaling corrects them -- times `lp`, over `qe` and
the partiality, over `inverse_scale_factor` -- and each is set against the
weighted mean of its reflection's CLEAN equivalents, itself left out: fully
recorded, nothing of the profile on masked pixels, away from the scan's ends.
Against every equivalent the bias would leak into the reference: with a fifth of
the observations partials 10 per cent high, partials read low and the clean ones
high. Against clean equivalents a planted bias reads as itself
(`tests/test_equivalents.py`).

It bins the relative difference by partiality, by `profile.measured` -- the
fraction of the profile on valid pixels, below 1 where a box crosses a module
gap -- by whether the foreground reached a masked pixel, by the images from each
end of the scan, and by resolution and I/sigma as controls. By default each
table counts only the observations clean in every other respect, so that it
measures its own cause; `--all` counts them all. Read the median and the ratio
of the sums: the mean of a ratio is pulled up by noise. Profile fitted and
summed intensities are reported apart (`--intensity`), and `--anomalous` keeps
Friedel mates apart.

On the 300 image insulin sweep, profile fitted: partials of partiality 0.8 to
0.9 +9.5 per cent, 0.9 to 0.99 +2.4, fully recorded +0.2; 80 to 95 per cent of the
profile on valid pixels -3 to -4; and no fully recorded, unmasked observation in
the first or last two images at all -- the excess at the scan's ends is the
partials'. `mxeq trend`, comparing two integrations, now bins by
`profile.measured` too.

## No thresholds

Version one applies no pass/fail criteria at all. Every check prints
measurements; a human decides what is acceptable. This is deliberate: choosing
a tolerance before seeing the distribution on real data means choosing it by
guessing. The `--json` output carries the same numbers as the text, so
thresholds can be written against it once there is something to write them
against.

Exit status is zero whenever the comparison ran. Non-zero means it could not
run, never that the two disagreed.

## Two ways to run the pipeline

Both are needed, and they answer different questions.

**Pinned.** Each stage is given DIALS' upstream output as its input. This
isolates one stage: a difference in the integrated intensities is the
integrator's, because both integrators were handed the same experiment model
and the same predictions.

**Cascade.** Each stage is given the previous stage's own output. This is the
number that matters, because it is what a user of the pipeline gets. It is also
the one that is hard to read without the pinned run beside it: a hundredth of a
degree of orientation difference at refinement moves every shoebox by a
fraction of a pixel, and the integrator then looks wrong when it is not.

Run pinned first. When pinned agrees and cascade does not, the compounding is
the finding.

## Install

```sh
pip install -e .
```

Dependencies are numpy, scipy, msgpack and gemmi. Notably **not** cctbx or
DIALS: `.refl` is msgpack and `.expt` is JSON, and both are read directly, so
the checks run in a container with no crystallographic software in it at all.

## When a file will not read

```sh
mxeq inspect strong.refl
```

prints the columns and their C++ types, and if the file is not the shape this
expects, falls back to walking the msgpack document without assuming any of the
key names. That output is the bug report.

## The boundaries

| boundary | join | headline metric |
| --- | --- | --- |
| `strong` | spatial, mutual nearest neighbour | matched fraction, and what did *not* match |
| `indexed` | spatial, then reindexing operator, then Miller index | agreement after reindexing |
| `refined` | none; models compared directly | cell, misorientation, detector distance |
| `integrated` | (Miller index, entering, id) | relative difference and pull, in resolution shells |
| `scaled` | as integrated, plus a merge into the ASU | CC-half of each, and CC between them |

`docs/boundaries.md` explains why each join is the way it is.

## Reindexing

Two indexing runs on the same images can produce the same lattice in a
different basis. `mxeq check indexed` finds the operator from the data --
spatially match the spots first, then search the candidate operators for the
one that maps A's indices onto B's -- and prints it. Pass it to the later
boundaries:

```sh
mxeq check integrated a.refl b.refl -e a.expt --operator 0,1,0,0,0,1,1,0,0
```

Without it the keyed join fails, which the tool says plainly rather than
reporting a near-total disagreement in intensity.

## Status

The `.refl` format is validated against real `dials.find_spots` and
`mxi_find` output, and `tests/data` holds a 48-row cut of it so
that is checked on every run. `.expt` is exercised only against synthetic
files so far.

The first version of the reader was wrong about the format in two ways and its
entire test suite passed, because the only format check was a round-trip
against its own writer. `CLAUDE.md` records that; it is the reason
`tests/test_format.py` asserts on raw bytes.
