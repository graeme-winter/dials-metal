// Predicting where reflections appear on a rotation scan.
//
// The problem for one reflection: a reciprocal lattice point r0 = A h is
// carried round by the goniometer, and diffracts at the rotation angles where
// it crosses the Ewald sphere. Writing R(phi) = S Rot(m2, phi) F, the
// condition |s0 + R(phi) r0| = |s0| reduces, because |R r0| does not depend on
// phi, to
//
//     |r0|^2 + 2 s0 . R(phi) r0 = 0
//
// and decomposing r0 about the rotation axis turns the right-hand term into
// B cos(phi) + C sin(phi) + A0. So each reflection needs one atan2 and one
// acos, with two roots -- the point enters the sphere and later leaves it --
// and no iteration anywhere. That closed form is why this is the stage worth
// putting on the GPU first: one thread per candidate hkl, fixed work, no
// divergence beyond the reflections that never intersect at all.
//
// This implementation is the scalar reference the device version has to agree
// with, so it is written for clarity and checkability, not for speed.

#pragma once

#include <cstdint>
#include <vector>

#include "geometry.hh"

namespace mxi {

struct Prediction {
  int h = 0, k = 0, l = 0;
  bool entering = false;
  std::size_t panel = 0;
  double phi = 0.0;                       // radians
  double px_fast = 0.0, px_slow = 0.0, z = 0.0;
  Vec3 s1;
};

struct PredictOptions {
  // Highest resolution to predict to, in Angstrom. Zero means the limit set by
  // the wavelength, where the Ewald sphere can no longer be reached.
  double d_min = 0.0;
  // Predict outside the scan range. Off by default: a prediction at an angle
  // the scan never visited is not an observable thing, and including them
  // quietly inflates every completeness figure computed downstream.
  bool allow_outside_scan = false;
  // Cap on the Miller index box, as a guard against a wildly wrong cell
  // turning a prediction into an out-of-memory.
  int max_index = 500;
  //: Threads to predict on: 0 for one per core, 1 for none.
  //:
  //: One independent unit of work per h, and nothing in the body touches
  //: anything outside it. The results are concatenated in h order, so the
  //: output is identical whatever the thread count -- which matters because
  //: everything downstream is indexed by position in this vector.
  std::size_t threads = 0;
};

// The two rotation angles at which a reciprocal lattice point meets the Ewald
// sphere, if it meets it at all. Exposed separately from `predict` because it
// is the part with the mathematics in it and therefore the part worth testing
// on its own.
struct Intersections {
  bool any = false;
  double phi[2] = {0.0, 0.0};
  bool entering[2] = {false, false};
};
Intersections ewald_intersections(const Experiment &e, const Vec3 &r0);

// Predict every reflection of `crystal` observable on `e`'s scan and detector.
std::vector<Prediction> predict(const Experiment &e,
                                const PredictOptions &options = {});

// Predict a given list of Miller indices, whether or not they are observable.
// Used by the tests, and by anything that wants to ask "where would this one
// have been" without enumerating a whole sphere.
std::vector<Prediction> predict_indices(
    const Experiment &e, const std::vector<std::array<int, 3>> &indices,
    const PredictOptions &options = {});

// The Miller index box that covers a given resolution, from the real-space
// cell: |h_i| <= |a_i| / d_min, because h_i is the dot product of the i-th
// real-space basis vector with the scattering vector.
std::array<int, 3> index_bounds(const Crystal &crystal, double d_min,
                                int max_index);

}  // namespace mxi
