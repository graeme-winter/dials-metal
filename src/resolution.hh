#pragma once

// A resolution limit from CC half, as dials.estimate_resolution finds it
// (dials/util/resolution_analysis.py): CC half in 50 bins of equal count, a
// tanh in d*^2 fitted to it, and where the fit crosses 0.3; and where CC half
// stops being significant, from a logistic curve fitted to which bins are.

#include <cstddef>
#include <vector>

#include "scale.hh"

namespace mxi {

struct ResolutionBin {
  double d_min = 0.0;     //: the bin's high-resolution limit
  double d_star_sq = 0.0; //: 1 / d_min^2
  double cc_half = 0.0;
  std::size_t n = 0;     //: reflections with two or more observations
  double critical = 0.0; //: CC half above this is significant at the level
  bool significant = false;
};

//: CC half in bins of equal count of unique reflections, sorted by d, as
//: iotbx's counting_sorted binning; each reflection's observations split at
//: random into halves, their plain means correlated. d is the crystal's at the
//: unique index. Reflections whose normalised merged intensity E^2 is 16 or
//: more
//: -- Wilson outliers -- are left out, as dials.estimate_resolution leaves
//: them.
std::vector<ResolutionBin> cc_half_bins(const ScaleData &data,
                                        const std::vector<double> &g,
                                        const Crystal &crystal, int bins = 50,
                                        std::size_t min_per_bin = 10,
                                        double significance_level = 0.1,
                                        std::size_t *wilson_outliers = nullptr);

struct ResolutionEstimate {
  //: 0 where there is none: the fit crossed nothing, or nothing was
  //: significant.
  double d_min_cc_half = 0.0, d_min_significance = 0.0;
  double r = 0.0, s0 = 0.0; //: the tanh: CC half = (1 - tanh((x - s0) / r)) / 2
  bool fitted = false;
};

//: dials.estimate_resolution's two limits. The tanh is fitted by weighted least
//: squares, sigma 1/sqrt(n - 3), from r = 0.2 and s0 = 0.4; the limit is where
//: the FITTED values cross `limit`, by linear interpolation between bins; or,
//: where CC half is above it in every bin with pairs enough to fit (n > 3), or
//: the fit never comes down to it, the last such bin. dials.estimate_resolution
//: counts every bin, and a bin with no pairs -- CC half zero -- left it none.
//: The significance limit is 1/sqrt(res) for 1 - expit(r (x - res)) fitted to
//: which bins are significant.
ResolutionEstimate estimate_resolution(const std::vector<ResolutionBin> &bins,
                                       double limit = 0.3);

//: The upper p quantile of Student's t with nu degrees of freedom, by the
//: Cornish-Fisher expansion about the normal: better than 1e-4 from nu = 10,
//: and the bins here have hundreds.
double student_t_quantile(double p_upper, double nu);

} // namespace mxi
