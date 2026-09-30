#pragma once

// Post-refinement: refining the models against the centres INTEGRATION
// measured, rather than the spot finder's, and integrating again with them.
//
// Refinement is ordinarily fitted to the spot finder's centres, which depend on
// strength: for a weak spot the finder sees only the pixels above its
// threshold, the peak of a rocking curve that has a tail toward earlier images,
// and for a strong one the whole curve. The fitted model sits between, and
// every strong spot is then early by the difference -- about 0.1 images, with a
// period of 180 degrees of rotation on a 3600 image sweep. The integrator's
// centre is over the whole profile, and a model fitted to it removes that: on a
// 300 image sweep, from -0.115 images to -0.008, the drift along the scan gone
// flat.

#include <cstddef>
#include <vector>

#include "expt.hh"
#include "refine.hh"
#include "refl.hh"

namespace mxi {

//: The rows of an integrated table fit to refine against: integrated by
//: summation, with a centre of mass (a finite xyzres.px.value -- where there is
//: none, xyzobs.px.value is the prediction, and a residual of exactly zero
//: would pull the refinement back to the model it started from), and indexed. A
//: reflection whose foreground reaches a masked pixel is not summed, and its
//: centre, pulled off by the gap, is left out with it.
std::vector<std::size_t> rows_for_postrefinement(const Table &integrated);

//: Control points for a scan-varying crystal by default: one per 10 degrees of
//: the scan, as the reference profiles have a scan block per 10 degrees, and at
//: least five -- 36 on a full turn, five on 30 degrees. mxi_refine's
//: --scan-varying with no number and --postrefine both use it.
std::size_t scan_varying_points(const Scan &scan);
std::size_t postrefinement_points(const Scan &scan);

struct PostrefineResult {
  std::size_t candidates = 0; //: rows of the integrated table
  std::size_t selected = 0;   //: of those, fit to refine against
  TwoPassRefinement refinement;
};

//: Refine `experiments` against the centres in `integrated`, as mxi_refine
//: would: scan-static, then scan-varying with `points` control points (0 for
//: postrefinement_points of the first experiment's scan), the detector held.
PostrefineResult postrefine(ExperimentList &experiments,
                            const Table &integrated, std::size_t points = 0);

} // namespace mxi
