# The spot finder in this tree

It came in by `git subtree` with its history and sat beside this project rather
than in it: its own repository-level files, its own project() and C++ standard,
its own JSON parser, its own reflection-table writer, and its tests invisible to
the top-level `ctest`.

## Where it is now

    src/spots/          its sources, beside src/ rather than under a project
    tests/spots/        its tests, run by the same ctest
    cmake/spots.cmake   its build, included by the top-level CMakeLists
    third_party/        the bitshuffle submodule
    vcpkg.json          at the root, checked against project() as before

There is one `project()`, one C++ standard and one build type, all the
parent's. The build is `include()`d rather than `add_subdirectory()`d, because
a subdirectory that is not a project has no reason to be one.

It is still guarded: the top-level decides whether HDF5 and the submodule are
both present before including anything, because the spot finder's build calls
`find_package(HDF5 REQUIRED)` and that would abort the whole configure rather
than skip one part.

## What is shared now

**The JSON parser.** The spot finder's expt.cc, then spotfinder/src/expt.cc and
now `src/spots/expt.cc`, was 461 lines, of which about 400
were a second parser for a format this repository already parses. It is now 152
and uses `src/json.cc`.

**It does NOT share the experiment reader**, and that is deliberate. The shared
reader refuses a scan with no oscillation, because assuming zero would make one
compare equal to anything -- a strictness that caught a real bug. The spot
finder needs the image range, the panel size and where the images are, and none
of those is the oscillation. Making it refuse an `.expt` it can work from
perfectly well would be a regression dressed up as tidying. So the parsing is
shared and the interpretation is not.

**Repository-level files.** `.clang-format` and `LICENSE` were inside
`spotfinder/` and now sit at the root, where they cover the whole tree.

**Tests.** `enable_testing()` is now called before the subdirectory is entered.
It was called after, so `ctest` at the top ran this project's tests and
silently none of the spot finder's.

## What is still duplicated

`src/spots/refl.cc` is 353 lines that write a reflection table, which
`src/refl.cc` also does and does more generally -- opaque columns, every dtype,
round-tripped against real DIALS files.

The duplication has already cost something: msgpack's 4 GB limit on a binary
was guarded in this writer and not in the other, and a 4.87 GB shoebox column
went through the unguarded one and was written corrupt.

It has not been consolidated because the two writers must produce the same
bytes for the same spots and that has not been demonstrated. The acceptance
test is written and is in `tests/spots/test_refl_golden.cc`: it holds the
current output of a fixed set of spots, so a rewrite can be shown to change
nothing before it is believed. **Do the consolidation against that test, not
against a reading of the code.**

The difference in field order between the two writers means byte identity may
not be achievable, in which case the test should compare decoded contents --
every column, every dtype, every shoebox byte -- rather than be weakened to
pass.

## The Python

There is one package, `python/`, and everything Python in this repository is in
it. It was in five places: the checker in `mxeq/`, two plotters in `docs/`,
three scripts in `spotfinder/tests/`, and two data makers in `tests/`.

    mxeq.checks     the equivalence checks, including the spot finder's two
    mxeq.plots      drawing what the C++ tools wrote; needs matplotlib
    mxeq.fixtures   making test data; needs h5py

The checker needs neither matplotlib nor h5py and must keep working where they
are absent, so nothing in it may import those two subpackages. They are extras
in `pyproject.toml` and there is a check that `mxeq.cli` imports with both
blocked.

`mxeq.checks.spotfinder_refl` carried its own sixty-line msgpack decoder,
because it began as a script beside the spot finder where no dependency was
allowed. Inside the package `msgpack` is already required and the hand-rolled
one is gone.

Consolidating found a real bug: `mxeq.plots.anisotropy` unpacked eleven columns
positionally from a file that had grown to eleven from nine. It had simply
stopped running, and nothing said so. It now reads the columns by name from the
header and refuses a file whose header and data disagree.

## Building it is the only check that works

For most of this work the spot finder could not be built here and the moves
were checked statically: that every include resolved against the include path
the build sets, and that every path named in `cmake/spots.cmake` existed. Both
checks passed. Both missed real breakage, because both were checking against my
description of what the build does rather than against the build:

* the version guard compares `vcpkg.json` against `project()`, and the
  top-level `project()` had no `VERSION`, so every configure aborted;
* the substitution that repointed `src/` at `src/spots/` required a trailing
  slash, so ten `target_include_directories(... PUBLIC src)` and two
  `PRIVATE tests` were left pointing at directories that no longer held the
  headers.

`apt-get install libhdf5-dev` fixes all of that, and should have been the first
thing tried rather than the last. The whole tree now builds here, including
`mxi_find`, and `ctest` runs eight suites.

What still cannot be checked here is anything needing a device -- the Metal and
CUDA kernels compile only on a machine that has them -- and anything needing
real images.

## Reading only the scan's frames, and giving them the scan's z

Spot finding 1800 images of a 36000 image run printed

    warning: imported.expt covers 1800 images and the series has 36000; z
    will be wrong for anything outside the scan

and then thresholded all 36000 -- twenty times the work, producing 1115233
spots of which a twentieth belonged to the experiment asked about. The warning
described the bug instead of preventing it.

`Series::restrict_frames` now limits what a series offers, and the spot finder
asks for the scan's frames. For NXmx the key is the file's array index, so the
restriction is a filter on it and costs nothing. A series that cannot restrict
itself is an error rather than a silent read of everything, since spots outside
the scan land at z values the experiment does not have. A scan claiming more
images than the file holds is an error too, where it used to be the same
warning.

**Restricting it exposed a second bug, which had been there all along.** z was
`frame number + (first_image - 1)`, on the assumption that a sliced import's
frames are numbered from zero. That holds for a file containing only the slice;
it is false for a master file covering the whole run, where image 6 already
arrives as frame 5. So a scan of images 6 to 10 put its spots at z 10.5 to
14.5 instead of 5.5 to 9.5. Reading every frame had hidden this, because the
frames outside the scan were exactly where the doubled offset sent things.

The `.expt` says which file index each scan image is -- `single_file_indices`
-- so the offset is `(first_image - 1) - first_index`: nothing when the two
line up, five for a slice file starting at zero. Without `single_file_indices`
the file index of image n is taken as n - 1, which is the same assumption the
frame restriction makes, so the frames read and the z they are given cannot
disagree.

The tests build both layouts and run the real binary. On the old code the
whole-run test fails and the slice-file test passes, which is the point of
having both: one catches this, the other catches breaking the case that
already worked.


## The name

The spot finder is `mxi_find`, like the rest of the pipeline's programs. It was
`dials-metal-find-spots`, from when it was a separate project. Nothing writes
the program's name into a file it produces -- it appears only in `--version`
-- so tables written under the old name read exactly as before. A binary built
under the old name is not removed by rebuilding and should be deleted by hand,
or it will go on being found first on a PATH.
