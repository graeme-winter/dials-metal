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

#include "expt.h"
#include "geometry.h"
#include "refl.h"

namespace mxi {

struct RefineOptions {
  bool crystal = true;
  bool detector = true;
  bool beam = false;
  //: One crystal for every experiment, or one each.
  bool shared_crystal = true;
  int max_iterations = 30;
  //: Reject reflections whose residual exceeds this many times the robust
  //: spread, between macrocycles. Zero disables rejection entirely.
  double outlier_sigma = 4.0;
  int macrocycles = 3;
  //: Stop when the weighted residual improves by less than this fraction.
  double convergence = 1e-6;
  bool verbose = false;
};

struct RefineResult {
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

// Write xyzcal.px and xyzcal.mm into the table from the current models, so the
// result can be compared against DIALS' own predictions.
void update_predictions(const ExperimentList &experiments, Table &reflections);

}  // namespace mxi
