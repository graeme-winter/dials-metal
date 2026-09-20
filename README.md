# dials-metal

An independent implementation of the MX data-reduction chain downstream of the
images: spot finding, indexing, refinement, and the beginnings of integration.
C++ with no dependencies outside the spot finder, plus `mxeq`, a checker that
compares its output against DIALS.

The aim is to **stand in for DIALS**, not to improve on it. A pipeline built on
its own ideas would be fast and would agree with nothing anyone runs. Where it
departs from DIALS deliberately, the departure is a run-time option that
defaults off and the reason is written down.

## Layout

| | |
| --- | --- |
| `spotfinder/` | the Metal and CUDA spot finder; needs HDF5 |
| `src/`, `apps/` | indexing, refinement, prediction, profile model |
| `python/` | everything Python: the checker, the plots, the fixtures |
| `docs/` | notes on integration and on a device port |
| `CLAUDE.md` | working notes: what was got wrong, and how it was found |

`mxeq` is the referee and must stay independent of what it judges. Nothing in
`src/` may import it and it must never import them; it reads files. Within it,
the checker must not import `mxeq.plots` or `mxeq.fixtures`, so it keeps
working where matplotlib and h5py are not installed.

## Building

```sh
git submodule update --init --recursive    # for the spot finder
cmake -S . -B build && cmake --build build
./build/mxi_tests
```

The spot finder is built only when HDF5 and the submodule are both present;
otherwise `cmake` says so and builds the rest, so a machine that cannot build
the whole tree still builds and tests most of it.

## The chain

```sh
dials-metal-find-spots  master.h5                          # -> strong.refl
mxi_index               imported.expt strong.refl          # -> indexed.*
mxi_refine              indexed.expt  indexed.refl --analytic
mxeq check indexed      indexed.refl  dials/indexed.refl -e indexed.expt
```

Output is interchangeable with DIALS at every boundary: `dials.refine`,
`dials.integrate` and `dials.export` all read it, and `dials.image_viewer`
draws it.

## The tools

| | |
| --- | --- |
| `mxi_index` | FFT indexing with assign, refine and reassign macrocycles |
| `mxi_refine` | scan-static and scan-varying refinement, analytical derivatives |
| `mxi_residuals` | residuals by resolution, by detector module, by anything |
| `mxi_profile` | the Gaussian profile model, and how much of a spot it holds |
| `mxi_mask` | shoeboxes marking the integration region, for the image viewer |
| `mxi_grid` | spot density in Kabsch space, and spot widths across the face |
| `mxi_forward` | the model rendered onto the pixels, against the data |
| `mxeq check` | two pipelines compared at a boundary, over a common set |

## Where it stands

**Indexing and refinement are done and agree with DIALS.** On 1800 images of
insulin, given the same reflections, refinement reproduces `dials.refine`'s
detector distance to two microns and its cell to six thousandths of an
Angstrom. Analytical derivatives for crystal, detector and beam are validated
against finite differences and are six times faster.

**Integration is not started.** `docs/integration.md` has the plan and what it
is bound by. The profile model is implemented and does not yet agree with
DIALS; `mxi_profile` prints both values and the disagreement with every answer,
which is the only place either number should be read from.

**A device port is designed but not written.** `docs/gpu.md`. The target
evaluation is 88 to 106 per cent of refinement time, the work is one
independent thread per reflection and parameter, and float32 holds seven digits
in an analytical derivative and none at all in a finite difference -- which is
why the analytical ones exist.

## Two things this repository is really about

**A round-trip test cannot detect a wrong convention.** A file written and read
by the same code agrees with itself whatever it means. Five separate format
bugs here were invisible until DIALS read the output: the `.refl` container
shape, the geometry conventions, `0` written where `0.0` was required, blocks
dropped on write, and flags never set. Everything that matters is checked
against a file somebody else wrote.

**When what you are measuring is at the sampling scale, render the model into a
measurement rather than comparing a measurement to a model.** Spots here are
two pixels across. Comparing their moments with a Gaussian in closed form
suggested the model was wrong by factors and that the spots were anisotropic;
integrating the same model over the same pixels inside the same mask put the
positions within a hundredth of a pixel, the widths within four per cent, and
predicted the anisotropy. Most of the disagreement had been the comparison.
