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

## Usage

**What you need:** the images, as an NXmx HDF5 master file and its data files.
`mxi_import` reads the master and writes `imported.expt`, as `dials.import`
does -- `docs/import.md`, and `dials.import`'s own output serves as well.

### From images to scaled data

```sh
mxi_import    master.nxs                                    # -> imported.expt
mxi_find      imported.expt -o strong.refl                  # spots
mxi_index     imported.expt strong.refl                     # -> indexed.expt, .refl
mxi_refine    indexed.expt indexed.refl --analytic          # -> refined.expt, .refl
mxi_integrate refined.expt refined.refl                     # -> integrated.expt, .refl
mxi_symmetry  integrated.expt integrated.refl               # -> symmetrized.expt, .refl
mxi_scale     symmetrized.expt symmetrized.refl --d-min-auto  # -> scaled.expt, .refl
```

Each step reads the previous one's pair of files and writes its own, under the
names shown unless `--output-expt` and `--output-refl` (or `-o`) say otherwise.
`dials.merge` takes `scaled.refl` to make a merged MTZ file.

On the 3600 images of an EIGER2 XE 16M sweep -- `ins10_1.nxs` of
https://zenodo.org/records/8376818 -- on an M4 Max MacBook, that chain with
`--gpu` for finding and integrating and `--beam` in refinement takes 36.6 s of
wall time and 3m42 of CPU (30 September 2026).

### What each step does, and the options most often wanted

| program | does | options worth knowing |
| --- | --- | --- |
| `mxi_find` | finds spots on every frame, of the images an experiment list names or a master given itself | `-j N` threads (every core), `-g`/`--gpu` the threshold on the GPU, `--no-shoeboxes` a much smaller table, `--min-spot-size N` |
| `mxi_index` | indexes by 3D FFT, reduces the cell, refines | `--d-min D`, `--max-cell A`, `--verbose` the search |
| `mxi_refine` | refines beam, detector and crystal; scan-varying by default, one control point per 10 degrees | `--analytic` analytical derivatives (recommended), `--beam` the beam direction too, `--static` one crystal setting for the scan, `--scan-varying N` N control points |
| `mxi_integrate` | predicts, integrates by summation and profile fitting, in one pass over the images | `-g`/`--gpu` profile fitting on the GPU, `--threads N` (default every core), `--d-min D`, `--postrefine` refine against integration's own centres and integrate again, `--summation-only`, `--save-shoeboxes` |
| `mxi_symmetry` | determines the Laue group and reindexes | `--max-delta D` the lattice's obliquity tolerance, `--threads N` |
| `mxi_scale` | scales, with an error model, and reports merging statistics | `--d-min-auto` cut where CC half falls to 0.3, `--d-min D`, `--space-group NAME`, `--l-max L` absorption surface degree, `--anomalous` Friedel mates apart for a strong anomalous signal, `--threads N` |

Every program takes `--help` for the rest, and `--version`.

### Threads and the GPU

Where a program works in parallel it uses every core by default; `-j` or
`--threads` sets fewer. Built with
`-DSPOTFINDER_METAL=ON` (Apple silicon) or `-DSPOTFINDER_CUDA=ON` (NVIDIA), two
steps can use the GPU. `mxi_find --gpu` runs the threshold there, with the same
spots as the CPU (16-bit images only under Metal, which has no double
precision; for 32-bit images whose counts never reach 0xFFFD -- `mxi_max` says
-- `--gpu-force` narrows them to 16 bits, the same spots). `mxi_integrate --gpu` fits profiles there, in single precision
where the CPU's fit is double -- the same arithmetic run on the CPU puts its
intensities a median of 9e-6 sigma from the double fit's, none above 0.1 sigma,
on a 300 image sweep -- and on the 16M sweep above took integration from 20.0 s
to 14.6 on an M4 Max. A
build or machine without a GPU says so and runs on the CPU. `docs/gpu.md` has the
measurements.

### Logs and timing

Each program prints its report to standard output and writes the same report
to `mxi_<program>.log` where it runs -- `mxi_find.log`, `mxi_integrate.log` --
as DIALS writes `dials.find_spots.log`. Warnings and errors go to standard error
and into the log too, in order, so the log of a run that failed says why. A run
that only asks for `--help` or `--version` writes no log, rather than overwrite
the last real run's.

Every program takes `--timing`, and ends with one table of where its time went:
phases in seconds and per cent of the run, and where work runs in parallel --
the spot finder's reading, decompressing and thresholding, the integrator's
frame reading -- time summed across the threads against the time they had. To
time a whole chain step by step, wall clock and CPU, bash's own `time` will do:

```sh
TIMEFORMAT='%R s wall, %U s user, %S s system'
time mxi_find imported.expt -o strong.refl > find.out
```

### Mixing with DIALS

Output is interchangeable with DIALS at every boundary: any step can be swapped
for DIALS' own -- `dials.symmetry` and `dials.scale` take `integrated.*` as
readily as `mxi_symmetry` and `mxi_scale` do -- `dials.export` reads the
integrated table, and `dials.image_viewer` draws every step's output.
To compare a step against DIALS, see `mxeq` below and `python/README.md`.

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

## The tools

The pipeline:

| | |
| --- | --- |
| `mxi_import` | an experiment list from an NXmx master, as dials.import writes it; see `docs/import.md` |
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
| `mxi_max` | the largest count on a valid pixel of a series, and whether it fits in 16 bits -- for `mxi_find --gpu-force` |
| `mxi_readtest` | how fast the images read and decompress, apart from the programs: `--direct-chunk` through HDF5 as they read, or `--pread` outside its lock, a checksum to show both read the same |

And `mxeq`, which judges the output: `check` compares two pipelines at a
boundary; `trend`, `html`, `explain` and `disagree` find and explain where two
integrations differ; `residuals` shows how well positions were predicted; and
`profiles` draws the learned reference profiles. Every program answers
`--help`.

## Where it stands

**The chain needs nothing of DIALS.** `mxi_import` (`docs/import.md`) writes the
experiment list from an NXmx master as `dials.import` does -- on insulin the same
in every model but the exposure time -- and `mxeq compare-expt` checks it against
`dials.import` on each new data set.

**Indexing and refinement are done and agree with DIALS.** On 1800 images of
insulin, given the same reflections, refinement reproduces `dials.refine`'s
detector distance to two microns and its cell to six thousandths of an
Angstrom. Analytical derivatives for crystal, detector and beam are validated
against finite differences and are six times faster.

**Integration works and agrees with DIALS through scaling.** On 1800 images
of insulin `dials.scale` gives Rmerge 0.038 and Rpim 0.009 on this output and on
DIALS' alike; on a 3600 image Eiger 16M sweep the merging statistics are close
to DIALS'. The 16M sweep -- 1.08 million
reflections -- integrates in 19.4 s on a MacBook's 16 threads, 14.6 with
profile fitting on its GPU. `docs/integration.md` is the reference, and lists
what is still open -- among it two biases `mxeq equivalents` measures against the
data's own equivalents: partials high, likely because our sigma_m is too wide
(the partials suggest DIALS's value), and reflections crossing a module gap 3 to
5 per cent low in profile fitting.

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

**Profile fitting runs on the GPU too**, `mxi_integrate --gpu`, in single
precision, overlapped with reading frames: on the MacBook's Metal the 16M sweep's
integration takes 14.6 s against the CPU's 20.0; on CUDA it is a smaller gain,
the cost moving the boxes to the device. A device port of the refinement target
is designed and not written. `docs/gpu.md` has the measurements and the plan.

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
