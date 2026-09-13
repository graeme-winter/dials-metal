# mxeq

Equivalence checks between two runs of an MX processing pipeline, boundary by
boundary, without cctbx.

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
`dials-metal-find-spots` output, and `tests/data` holds a 48-row cut of it so
that is checked on every run. `.expt` is exercised only against synthetic
files so far.

The first version of the reader was wrong about the format in two ways and its
entire test suite passed, because the only format check was a round-trip
against its own writer. `CLAUDE.md` records that; it is the reason
`tests/test_format.py` asserts on raw bytes.
