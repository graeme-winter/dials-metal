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

**The JSON parser.** `spotfinder/src/expt.cc` was 461 lines, of which about 400
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

`spotfinder/src/refl.cc` is 353 lines that write a reflection table, which
`src/refl.cc` also does and does more generally -- opaque columns, every dtype,
round-tripped against real DIALS files.

It has not been consolidated because the two writers must produce the same
bytes for the same spots and that has not been demonstrated. The acceptance
test is written and is in `spotfinder/tests/test_refl_golden.cc`: it holds the
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

## What cannot be checked here

This container has no HDF5, so `spotfinder/` does not build in it and nothing
that touches images can be run. `src/expt.cc` and `src/refl.cc` have no HDF5
dependency and their tests build and run standalone, which is why those two
were the parts consolidated: they are the parts that can be held to a test.
