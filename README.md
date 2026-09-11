# dials-metal-index

Indexing, refinement and prediction for the standalone Metal pipeline, in C++
with no dependencies. Companion to `dials-metal-find-spots`: it takes a
`strong.refl` and produces an `indexed.expt` and `indexed.refl`, so the chain
runs end to end without DIALS in the middle.

Name is provisional. Most of this is not Metal work and does not need to be.

## What is here

| | |
| --- | --- |
| `src/linalg.h`, `.cc` | Vec3, Mat3, rotations, a Cholesky solver |
| `src/geometry.h`, `.cc` | Beam, detector, goniometer, scan, crystal; the map from a spot to a reciprocal lattice point |
| `src/predict.h`, `.cc` | Reflection prediction on a rotation scan |

## What is not here yet

Indexing and refinement. Prediction is first on purpose — see below.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/mxi_tests
```

No third-party libraries. C++20, one static library, one test binary.

## Why prediction first

Prediction is the oracle for everything after it. Given a crystal it generates
a reflection list with known `hkl` and known centroids; throw the indices away,
feed the centroids to the indexer, and the indexer must recover the cell and
orientation it started from. That closes a loop against ground truth, needs no
DIALS, and runs in CI.

Build indexing first and there is nothing to test it with except real data that
can only be checked by eye.

It is also the stage most worth putting on a GPU. Each reflection needs one
`atan2` and one `acos` with no iteration, and nothing in the body of the loop
touches anything outside it — one thread per candidate `hkl`, fixed work, and
the only shared state is the output.

## The thing to be careful about

Every geometry convention in `src/geometry.h` is a **belief** about DIALS and
has not been checked against a real `.expt`. The closed-loop test cannot catch
a convention error, because an error present in both the forward and reverse
map cancels exactly.

`docs/conventions.md` lists all seven, what each would break, and the test that
would falsify it. Read it before trusting any output against DIALS.
