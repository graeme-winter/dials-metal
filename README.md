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
| `src/`, `apps/` | indexing, refinement, prediction, profile model |
| `src/spots/` | the Metal and CUDA spot finder; needs HDF5 |
| `tests/` | the C++ tests, `tests/spots/` for the spot finder's |
| `python/` | everything Python: the checker, the plots, the fixtures |
| `docs/` | notes on integration and on a device port |
| `CLAUDE.md` | working notes: what was got wrong, and how it was found |

`mxeq` is the referee and must stay independent of what it judges. Nothing in
`src/` may import it and it must never import them; it reads files. Within it,
the checker must not import `mxeq.plots` or `mxeq.fixtures`, so it keeps
working where matplotlib and h5py are not installed.

## Building

The pipeline itself has no dependencies:

```sh
cmake -S . -B build && cmake --build build
ctest --test-dir build
```

The spot finder needs HDF5 and the bitshuffle submodule, and is left out
silently if either is missing -- `cmake` prints which. To build it:

```sh
sudo apt-get install libhdf5-dev        # Debian and Ubuntu
brew install hdf5                       # macOS
# or: conda install -c conda-forge hdf5

git submodule update --init --recursive
cmake -S . -B build && cmake --build build
```

`cmake` says `spotfinder: building` when it has both, and
`spotfinder: HDF5 not found, skipping` or `bitshuffle submodule not checked
out` when it does not. If the spot finder is what you want, check that line
rather than the absence of an error: leaving it out is not a failure.

The threshold kernels are CPU by default. One device backend at a time, and
`cmake` refuses both at once:

```sh
cmake -S . -B build -DSPOTFINDER_CUDA=ON       # NVIDIA

cmake -S . -B build -DSPOTFINDER_METAL=ON \
      -DMETAL_CPP_DIR=~/third-party-git/metal-cpp   # Apple
```

Metal needs `METAL_CPP_DIR`. metal-cpp is Apple's header-only wrapper over the
Metal API and is distributed as a zip from
<https://developer.apple.com/metal/cpp/> rather than through any package
manager, so it has to be pointed at rather than found. `METAL_CPP_DIR` also
works as an environment variable, and the directory wanted is the one holding
`Metal/Metal.hpp`. Without it the configure stops and says so.

CUDA needs nothing pointed at, but it does guess what to build for:
`CMAKE_CUDA_ARCHITECTURES` defaults to `native` on CMake 3.24 and later and to
`70` before that, which is a guess rather than an answer. Set it to the cards
the binary will actually run on -- `-DCMAKE_CUDA_ARCHITECTURES="80;90"` -- if
that is not this machine.

## Making indexing faster

The transform is the largest phase of indexing on a fast machine. FFTW is about
three times quicker than the built-in radix-2 and is opt-in, because it is GPL
and the licence of anything linking it has to accommodate that:

```sh
sudo apt-get install libfftw3-dev       # Debian and Ubuntu
brew install fftw                       # macOS
# or: conda install -c conda-forge fftw

cmake -S . -B build -DMXI_FFTW=ON
```

`cmake` prints `FFT: FFTW` with how it found it, or `FFT: the built-in
radix-2`. If it cannot be found -- a Homebrew prefix that is not searched, no
pkg-config -- point at it:

```sh
cmake -S . -B build -DMXI_FFTW=ON \
      -DFFTW3_LIBRARY=/opt/homebrew/lib/libfftw3.dylib \
      -DFFTW3_INCLUDE_DIR=/opt/homebrew/include
```
 The built-in is
always compiled whichever is used, and there is a test that the two agree:
a library can be planned wrongly, normalised differently or linked against
another precision, and none of that shows in a result that still looks like a
lattice.

The Jacobian of the refinement target is the largest remaining phase of
indexing -- 36 per cent of it on a fast machine -- and is built on one thread
per reflection. `--jacobian-threads N` sets the count, 0 being one per core and
1 none:

```sh
mxi_index imported.expt strong.refl --jacobian-threads 1 --timing
```

The threaded and serial results are bit identical, and there is a test that
says so: every thread writes only its own reflections' entries, so there is no
shared accumulator and no reordered sum.

On x86, `MXI_SSE41` is on where the compiler takes it, for the rounding
instruction. Without it `std::rint` compiles to a sequence and the basis search
runs about three times slower. `-DMXI_SSE41=OFF` for a CPU older than 2009. ARM
needs nothing.

`SPOTFINDER_AVX2` is on where the compiler takes it, which means the binary
will not run on a CPU without AVX2. `-DSPOTFINDER_AVX2=OFF` gives a portable
SSE2 build.

`ctest` runs everything: this project's tests and the spot finder's, when the
spot finder was built.

## The Python

```sh
pip install -e python                   # the checker
pip install -e "python[plots,fixtures]" # and the plots and test-data makers
pytest python/tests
```

The checker deliberately needs neither matplotlib nor h5py, so it runs where
those are not installed.

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
