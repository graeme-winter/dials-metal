# dials-metal

**Reviewing this code?** Start with `docs/review.md`: what this is, how it differs
from DIALS and why, and where to find things.

An independent implementation of the MX data-reduction chain for one rotation
sweep, from the images to scaled intensities: spot finding -- on the CPU, or the
GPU through Metal or CUDA -- indexing, refinement, integration, symmetry
determination and scaling, in C++; plus `mxeq`, a checker that compares its
output against DIALS' and explains the differences. It reads and writes DIALS'
own files, so any step can be exchanged for DIALS' at its boundary.

The aim is to **stand in for DIALS**, not to improve on it. A pipeline built on
its own ideas would be fast and would agree with nothing anyone runs. Where it
departs from DIALS the departure is deliberate, measured against DIALS on the
same data, and written down with the reason -- some behind an option, some,
where the measurement said so, by default. `docs/review.md` lists them all.

## Layout

| | |
| --- | --- |
| `src/`, `apps/` | indexing, refinement, prediction, the profile model, integration, symmetry, scaling |
| `src/spots/` | the spot finder, on the CPU or through Metal or CUDA, and the NXmx image reader |
| `tests/` | the C++ tests, `tests/spots/` for the spot finder's |
| `tools/` | measurement harnesses kept for rerunning, such as planted spots |
| `python/` | everything Python: the checker, the comparisons, the plots |
| `docs/` | the references for integration, symmetry, scaling and the spot finder; the guide for reviewers; what is outstanding; `docs/README.md` maps them |
| `CLAUDE.md` | working notes: what was got wrong, and how it was found |

`mxeq` is the referee and must stay independent of what it judges. Nothing in
`src/` may import it and it must never import them; it reads files. Within it,
the checker must not import `mxeq.plots` or `mxeq.fixtures`, so it keeps
working where matplotlib and h5py are not installed.

## Building

Everything needs the gemmi submodule, for space groups; indexing and
refinement need nothing else:

```sh
git submodule update --init --recursive
cmake -S . -B build && cmake --build build
ctest --test-dir build
```

gemmi is pinned at a release (v0.7.5), and only its symmetry source is
compiled. A checkout without it stops at `cmake`, naming the command above,
rather than building part of the pipeline.

**Anything that reads images needs HDF5 and the bitshuffle submodule**: the spot
finder `mxi_find` and the integrator `mxi_integrate` both. Without them both
are left out *silently* -- `cmake` prints which is missing, and the build
otherwise succeeds. To build them:

```sh
sudo apt-get install libhdf5-dev        # Debian and Ubuntu
brew install hdf5                       # macOS
# or: conda install -c conda-forge hdf5

git submodule update --init --recursive
cmake -S . -B build && cmake --build build
```

`cmake` says `spotfinder: building` when it has both, and
`spotfinder: HDF5 not found, skipping` or `bitshuffle submodule not checked
out` when it does not. If you want the spot finder or the integrator, check
that line rather than the absence of an error: leaving them out is not a
failure.

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

In `mxi_refine`, `--jacobian-threads` and `--normal-threads` do the same for
refinement. The first is exact; the second is a reduction and a threaded run
differs from a serial one in the last bits, by about 1e-12 relative. Use
`--normal-threads 1` where that matters.

`--fft-threads N` uses FFTW's own threading where FFTW was built with it: 0 is
one per core, 1 none. The configure line says whether it found the threaded
flavour.

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
dials.import  master.nxs                                        # -> imported.expt
mxi_find      -e imported.expt -j 16 --gpu -o strong.refl
mxi_index     imported.expt strong.refl                         # -> indexed.*
mxi_refine    indexed.expt indexed.refl --analytic              # -> refined.*
              # scan-varying by default, one control point per 10 degrees --
              # nearly every real crystal moves; --scan-varying N for N,
              # --static for a crystal that does not; static under 10 degrees
mxi_integrate refined.expt refined.refl                         # -> integrated.*
              # --gpu: profile fitting on the GPU, in single precision
# or, refining against the centres integration measures and integrating again:
mxi_integrate refined.expt refined.refl --postrefine
mxi_symmetry  integrated.expt integrated.refl                   # -> symmetrized.*
mxi_scale     symmetrized.expt symmetrized.refl                 # -> scaled.*
```

`dials.symmetry` and `dials.scale` take `integrated.*` as readily, at either of the
last two steps.

On the 3600 images of an EIGER2 XE 16M sweep -- `ins10_1.nxs` of
https://zenodo.org/records/8376818 -- on an M4 Max MacBook, the chain from the
images to scaled data, `mxi_find -j 16 --gpu`, `mxi_index`, `mxi_refine --beam
--analytic`, `mxi_integrate --gpu`, `mxi_symmetry` and `mxi_scale --d-min-auto`,
takes 36.6 s of wall time and 3m42 of CPU (30 September 2026).

Each program prints its report to standard output and writes the same report
to `mxi_<program>.log` where it runs -- `mxi_find.log`, `mxi_integrate.log` --
as DIALS writes `dials.find_spots.log`. Warnings and errors go to standard error
and into the log too, in order, so the log of a run that failed says why. A run
that only asks for `--help` or `--version` writes no log, rather than overwrite
the last real run's.

Every program takes `--timing`, and ends with one table of where its time went:
phases in seconds and per cent of the run, and where work runs in parallel --
the spot finder's reading, decompressing and thresholding, the integrator's
frame reading -- time summed across the threads against the time they had.

Output is interchangeable with DIALS at every boundary: any stage can be
swapped for DIALS' own, `dials.scale` and `dials.export` read the integrated
table, and `dials.image_viewer` draws every stage's output. To compare a stage
against DIALS, see `mxeq` below and `python/README.md`.

## The tools

| | |
| --- | --- |
The pipeline:

| | |
| --- | --- |
| `mxi_find` | spot finding on the CPU, CUDA or Metal, from NXmx HDF5 |
| `mxi_index` | FFT indexing with assign, refine and reassign macrocycles |
| `mxi_refine` | scan-static and scan-varying refinement, analytical derivatives |
| `mxi_integrate` | summation and profile fitting; see `docs/integration.md` |
| `mxi_symmetry` | the Laue group and space group, and the data reindexed into them; see `docs/symmetry.md` |
| `mxi_scale` | scaling one sweep, the error model, merging statistics; see `docs/scaling.md` |

For looking inside it:

| | |
| --- | --- |
| `mxi_residuals` | refinement residuals by resolution, by detector module, by anything |
| `mxi_profile` | the Gaussian profile model, and how much of a spot it holds |
| `mxi_mask` | shoeboxes marking the integration region, for the image viewer |
| `mxi_grid` | spot density in Kabsch space, and spot widths across the face |
| `mxi_forward` | the model rendered onto the pixels, against the data |
| `mxi_background` | the robust background fitted to pixel values from a file |

And `mxeq`, which judges the output: `check` compares two pipelines at a
boundary; `trend`, `html`, `explain` and `disagree` find and explain where two
integrations differ; `residuals` shows how well positions were predicted; and
`profiles` draws the learned reference profiles. Every program answers
`--help`.

## Where it stands

**Indexing and refinement are done and agree with DIALS.** On 1800 images of
insulin, given the same reflections, refinement reproduces `dials.refine`'s
detector distance to two microns and its cell to six thousandths of an
Angstrom. Analytical derivatives for crystal, detector and beam are validated
against finite differences and are six times faster.

**Integration works and agrees with DIALS through scaling.** On 1800 images
of insulin `dials.scale` gives Rmerge 0.038 and Rpim 0.009 on this output and on
DIALS' alike; on a 3600 image Eiger 16M sweep the merging statistics are close
to DIALS'. Judged by their own symmetry equivalents, reflections crossing a
module gap come out better here than in DIALS. The 16M sweep -- 1.08 million
reflections -- integrates in 37 seconds on 16 threads. `docs/integration.md` is
the reference, and lists what is still open.

**Symmetry and scaling are written for one sweep.** `mxi_symmetry`
(`docs/symmetry.md`) chooses the Laue group as dials.symmetry does and the space
group from the absences, and reindexes; `mxi_scale` (`docs/scaling.md`) scales
after Beilsten-Edmands et al. (2020), with cubic B-splines for the scale and
decay and spherical harmonics for absorption. Each is compared with DIALS on the
same sweeps; the error model's differences from dials.scale's are findings
about dials.scale, listed in `docs/outstanding.md`.

```sh
mxi_symmetry integrated.expt integrated.refl     # symmetrized.expt, .refl
mxi_scale symmetrized.expt symmetrized.refl      # scaled.expt, .refl
```

**The spot finder runs on the GPU, Metal and CUDA, and the two agree exactly.**
On 3600 frames of 16M pixels Metal takes 13.9 s on a MacBook and CUDA 32 s on an
RTX 4060, whose kernels are its limit; a fused CUDA kernel for that is written
and verified on the CPU, not yet run (`docs/spots.md`, and item 35 of
`docs/outstanding.md`).

**A device port is designed but not written.** `docs/gpu.md`. The target
evaluation is 88 to 106 per cent of refinement time, the work is one
independent thread per reflection and parameter, and float32 holds seven digits
in an analytical derivative and none at all in a finite difference -- which is
why the analytical ones exist.

**What is open is one list**, `docs/outstanding.md`, grouped by where the work
is and marked where a number is known to be wrong; `docs/README.md` says what
every document is for.

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
