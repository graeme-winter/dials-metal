#pragma once

// What an integration produced, in the tables a user reads to judge it: one
// row per resolution shell, and overall, lowest and highest. Modelled on the
// summaries dials.integrate prints, limited to what this program measures --
// there is no overload or ice-ring column, because nothing here detects
// either, and a column of zeroes would say that it had looked.

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace mxi {

//: The columns a summary is made from, one entry per reflection. Every vector
//: has the same length; res_fast and res_slow are NaN where no centre of mass
//: was found, as xyzres.px.value is.
struct SummaryInput {
  std::vector<double> d;
  std::vector<std::int64_t> flags;
  std::vector<double> intensity_sum, variance_sum;
  std::vector<double> intensity_prf, variance_prf;
  std::vector<double> profile_correlation;
  std::vector<double> background;
  std::vector<double> partiality;
  std::vector<double> res_fast, res_slow;
};

//: One shell, or the whole of the data.
struct SummaryRow {
  double d_min = 0.0, d_max = 0.0;
  std::size_t n = 0;
  std::size_t full = 0, partial = 0; //: partiality of at least 0.99, or less
  std::size_t summed = 0, fitted = 0;
  std::size_t bad_foreground = 0, bad_background = 0;
  //: Means over the reflections they apply to: background and I/sigma(sum)
  //: over the summed, I/sigma(prf) and CC prf over the fitted, RMSD XY over
  //: those with a centre of mass.
  double background = 0.0;
  double isigma_sum = 0.0, isigma_prf = 0.0;
  double cc_prf = 0.0;
  double rmsd_xy = 0.0;
};

struct IntegrationSummary {
  //: Highest resolution first, as dials.integrate shows them.
  std::vector<SummaryRow> shells;
  SummaryRow overall;
  //: Between the summed and the fitted intensity, over reflections with both:
  //: a sanity check that the two methods measured the same thing.
  double pearson_sum_prf = 0.0, spearman_sum_prf = 0.0;
  std::size_t both = 0;
};

//: Shells of equal volume in reciprocal space, equal steps in 1/d^3, between
//: the lowest and highest resolution present. A reflection without a finite
//: positive d is left out of every row.
IntegrationSummary summarise(const SummaryInput &in, int shells = 10);

//: The two tables, as text.
void print_summary(std::FILE *out, const IntegrationSummary &s);

} // namespace mxi
