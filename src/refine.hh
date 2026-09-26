// Refining the experimental model against observed spot centroids.
//
// The target is the one DIALS uses: the difference between observed and
// predicted position in x and y on the detector and in the rotation angle,
// weighted by the inverse of the centroid variances so the sum is
// dimensionless. Gauss-Newton on the normal equations, with Levenberg
// damping when a step makes things worse.
//
// DERIVATIVES ARE NUMERICAL, and that is a decision rather than a shortcut.
// Fifteen parameters over thirteen thousand reflections costs sixteen
// predictions per iteration, which is a fraction of a second. Analytical
// derivatives would be perhaps twenty times faster and would introduce the
// single most dangerous class of bug a refinement engine can have: a wrong
// derivative does not crash, it converges smoothly to the wrong answer and
// reports a small residual while doing it. If profiling ever says this matters,
// the numerical version stays as the thing any analytical version is checked
// against.
//
// PARAMETERISATION
//
//   crystal   the nine elements of A. For a triclinic crystal that is exactly
//             the degrees of freedom available -- three of orientation and six
//             of cell -- so nothing is constrained and nothing is redundant.
//             Imposing a Bravais setting is a separate job, done afterwards.
//   detector  three for the panel origin, three for a small rotation of the
//             panel about its own centre. The rotation is applied to the fast
//             and slow axes rather than being accumulated into them, so the
//             axes stay orthonormal however many iterations run.
//   beam      two, for a small tilt of the incident direction. Off by default:
//             beam direction and detector position are strongly correlated on
//             a single sweep and refining both at once mostly moves the pair
//             around a valley.
//
// ON MULTIPLE SWEEPS
//
// `shared_crystal` decides whether the sweeps refine one crystal or one each.
// Indexing has to assume one, because that is what puts every sweep in a
// common basis. Refinement is where that constraint should be broken: no real
// goniometer returns to precisely the same place, so a single matrix cannot
// fit several sweeps, and on l-cysteine that costs a factor of a hundred in
// residual. Refine shared first so the sweeps stay in one basis, then split.

#pragma once

#include <cstddef>
#include <vector>

#include "expt.hh"
#include "geometry.hh"
#include "refl.hh"

namespace mxi {

struct RefineOptions {
  bool crystal = true;
  bool detector = true;
  bool beam = false;
  //: One crystal for every experiment, or one each.
  bool shared_crystal = true;
  //: Control points in A across each scan. One or zero means a static
  //: crystal. A scan-varying model absorbs the orientation and cell drift a
  //: real goniometer and a real crystal produce over a sweep, and on
  //: l-cysteine it is worth a factor of nearly three in residual -- so a
  //: comparison of detectors made with static crystals is measuring the
  //: crystals.
  //:
  //: Scan-varying refinement must follow a static one, never replace it: the
  //: control points start from the static answer, and starting them from an
  //: unrefined model lets them absorb errors that belong to the detector.
  std::size_t scan_points = 1;
  int max_iterations = 30;
  //: Reject reflections whose residual exceeds this many times the robust
  //: spread, between macrocycles. Zero disables rejection entirely.
  double outlier_sigma = 4.0;
  int macrocycles = 3;
  //: Stop when the weighted residual improves by less than this fraction.
  double convergence = 1e-6;
  //: Use analytical derivatives instead of finite differences.
  //:
  //: Same answer, and that is the point of having both: the numerical path is
  //: the oracle the analytical one is checked against, and it stays. The
  //: analytical path exists for single precision, where a finite difference
  //: spends most of its significance on the cancellation.
  bool analytic = false;
  bool verbose = false;
  //: Ignore the centroid variances and weight every residual equally.
  //:
  //: On a photon-counting detector the centroid variances are nearly constant
  //: -- dominated by the 1/12 pixel quantisation floor -- so this changes very
  //: little, which is itself worth being able to demonstrate rather than
  //: assume.
  bool unit_weights = false;
  //: Discard reflections whose rotation angle is not determined by the data.
  //:
  //: The denominator of Waterman eqn (40) is the volume of the parallelepiped
  //: formed by the rotation axis, the reciprocal lattice vector and the beam.
  //: It goes to zero for reflections near the rotation axis, where the angle
  //: at which they diffract is arbitrarily sensitive to the model -- and the
  //: Lorentz factor has the same asymptote, so their observed angular
  //: centroids are poorly determined as well. Keeping them puts noise into the
  //: rotation-angle residual that no model can fit.
  //:
  //: DIALS discards below 0.05 by default and this did not, which is most of
  //: why its rotation-angle residual was worse.
  double min_volume = 0.05;

  //: Build the model from reflections at or above the median strength.
  //:
  //: Weak, marginally indexed reflections bias a refinement and cannot be
  //: rejected as outliers, because they are internally consistent -- the model
  //: moves until they fit. On insulin, refining on everything puts the
  //: detector 0.27 mm further out than DIALS; on the stronger half, 0.01 mm.
  //: Off by default because discarding half the data should be asked for.
  bool strong_only = false;
  //: Multiply the weight on the rotation-angle residual. The balance between
  //: the positional and angular parts of the target is a choice, not a
  //: measurement, and two programs can make it differently.
  double z_weight = 1.0;
};

struct RefineResult {
  //: Reflections dropped because their rotation angle is not determined by the
  //: data. Reported because it changes which reflections the residual is
  //: averaged over, and a residual quoted over a different set is not
  //: comparable with anything.
  std::size_t n_ill_conditioned = 0;
  //: Rows of the reflection table the fit actually used, after outlier
  //: rejection and after the ill-conditioned ones were dropped. Reported so
  //: that `used_in_refinement` can be set from it: which reflections a residual
  //: was averaged over is not a detail, and this is where DIALS records it.
  std::vector<std::size_t> rows_used;
  //: Rows drawn into the fit and then rejected as outliers. Separate from the
  //: ill-conditioned ones, which were never candidates.
  std::vector<std::size_t> rows_rejected;
  std::size_t n_used = 0;
  std::size_t n_rejected = 0;
  //: RMS of observed minus calculated, in pixels, pixels and images.
  double rmsd_x = 0.0;
  double rmsd_y = 0.0;
  double rmsd_z = 0.0;
  int iterations = 0;
  bool converged = false;
};

// Observed minus calculated for one reflection, in pixels and images, plus
// whether a prediction could be made at all.
struct Residual {
  bool valid = false;
  double dx = 0.0, dy = 0.0, dz = 0.0;
};

// Predict where a reflection with the given index would be seen, choosing the
// Ewald root nearer the observed rotation angle. Exposed because it is the
// heart of the target function and worth testing on its own.
Residual centroid_residual(const Experiment &e, std::size_t panel, int h, int k,
                           int l, double px_fast, double px_slow, double z);

RefineResult refine(ExperimentList &experiments, const Table &reflections,
                    const RefineOptions &options = {});

//: How far the analytical Jacobian sits from the numerical one, over the
//: entries large enough to compare. Exposed so a test can check them against
//: each other on any model -- a single sweep, several sweeps, static or
//: scan-varying -- rather than only inferring agreement from a refined result,
//: which can hide a wrong derivative behind a flat minimum.
struct JacobianComparison {
  std::size_t compared = 0;
  double median_relative = 0.0;
  double percentile_99 = 0.0;
  double percentile_999 = 0.0;
  double worst_relative = 0.0;
  //: Entries where the two disagree by more than half, which for a derivative
  //: means they do not even agree on the sign.
  std::size_t grossly_different = 0;
};
JacobianComparison compare_jacobians(const ExperimentList &experiments,
                                     const Table &reflections,
                                     const RefineOptions &options);

// Write xyzcal.px and xyzcal.mm into the table from the current models, so the
// result can be compared against DIALS' own predictions.
void update_predictions(const ExperimentList &experiments, Table &reflections);

//: Add the columns dials.index produces alongside the Miller indices, which
//: are what dials.refine and dials.integrate then expect to find.
//:
//: `xyzobs.mm.value` above all: DIALS measures its refinement residual in
//: millimetres and radians, so this is the observation it minimises against,
//: and without it dials.refine stops at
//:
//:     The supplied reflection table does not have the required data column:
//:     xyzobs.mm.value
//:
//: **Written once, at indexing, and never recomputed.** It is derived from the
//: detector model, so refining the detector makes it stale -- and DIALS does
//: not recompute it either, which is why it must never be used as a join key
//: between two processing runs. Recomputing it here would be a difference from
//: DIALS dressed up as a correction.
//: Set the `indexed` bit on reflections that have a Miller index, and clear it
//: on those that do not.
//:
//: dials.* filters on the flags rather than on the indices, so a table with
//: correct Miller indices and empty flags processes perfectly and is then
//: invisible to every selection downstream.
void set_indexed_flags(Table &reflections);

//: Set the `used_in_refinement` bit on the reflections the last refinement
//: actually fitted, and clear it on the rest.
//:
//: Which reflections a residual was averaged over is not a detail, and this is
//: where DIALS records it. Without it there is no way to ask afterwards what
//: the number was computed from.
void set_refinement_flags(const RefineResult &result, Table &reflections);

void add_observed_columns(const ExperimentList &experiments,
                          Table &reflections);

//: s1, rlp, entering and imageset_id, which DIALS recomputes from the CURRENT
//: model rather than freezing at import.
//:
//: The split is not arbitrary and was measured, not assumed. Computing
//: xyzobs.mm from the model as imported reproduces DIALS bit for bit, while
//: computing it from the refined model does not; s1 and rlp are the other way
//: round. dials.index leaves the millimetre centroids exactly as
//: dials.find_spots wrote them and calls map_centroids_to_reciprocal_space
//: after refining, so the two groups of columns describe the same observations
//: through different models. Anything joining or comparing them has to know
//: which.
void add_reciprocal_columns(const ExperimentList &experiments,
                            Table &reflections);

//: Seconds spent building the Jacobian and accumulating the normal equations,
//: summed over every refinement since the process started. For deciding what
//: to move to a device: the two have entirely different shapes, one being a
//: pass over reflections and the other a reduction into a small matrix.
//: Threads used to build the ANALYTICAL Jacobian: 0 for hardware_concurrency,
//: 1 for none. It does nothing on the finite-difference path, which is the
//: default and which builds its columns one residual evaluation at a time.
//: Seconds in each part of refinement, summed since the process started.
//: Globals for the same reason the Jacobian's are: refine() is on a public
//: header and its signature is not worth changing to answer a question about
//: where the time goes.
extern double g_observations_seconds; //: building the target rows
extern double g_solve_seconds;        //: the damped solve and the step
extern double g_outlier_seconds;      //: rejecting outliers between cycles
extern double g_residual_seconds;     //: residuals outside the Jacobian

extern std::size_t g_jacobian_threads;

//: Threads for the normal equations: 0 for one per core, 1 for none.
//:
//: Unlike the Jacobian, this one is a reduction, and a reduction split across
//: threads sums its terms in a different order. Two runs at the same setting
//: agree bit for bit; a threaded run and a serial one agree to about 1e-13
//: relative and not exactly. Set it to 1 where that matters.
extern std::size_t g_normal_threads;
extern double g_jacobian_seconds;
extern double g_normal_seconds;

} // namespace mxi
