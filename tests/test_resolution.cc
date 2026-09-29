#include <cmath>

#include "../src/resolution.hh"
#include "check.hh"

namespace mxi {

TEST(students_t_quantiles_are_the_tabulated_ones) {
  // From the tables: t at the upper 0.1 with 10, 30 and 100 degrees of
  // freedom, and at 0.05 with 20.
  check::close(student_t_quantile(0.1, 10), 1.3722, 2e-4, "t(0.9, 10)");
  check::close(student_t_quantile(0.1, 30), 1.3104, 2e-4, "t(0.9, 30)");
  check::close(student_t_quantile(0.1, 100), 1.2901, 2e-4, "t(0.9, 100)");
  check::close(student_t_quantile(0.05, 20), 1.7247, 2e-4, "t(0.95, 20)");
}

namespace {

//: Fifty bins from 1/d^2 = 0.02 to 0.6, CC half on the tanh given, n each.
std::vector<ResolutionBin> planted_bins(double r, double s0, std::size_t n) {
  std::vector<ResolutionBin> bins;
  for (int k = 0; k < 50; ++k) {
    ResolutionBin b;
    b.d_star_sq = 0.02 + 0.58 * k / 49.0;
    b.d_min = 1.0 / std::sqrt(b.d_star_sq);
    b.cc_half = 0.5 * (1.0 - std::tanh((b.d_star_sq - s0) / r));
    b.n = n;
    const double t = student_t_quantile(0.1, static_cast<double>(n - 2));
    b.critical = t / std::sqrt(static_cast<double>(n) - 2.0 + t * t);
    b.significant = b.cc_half > b.critical;
    bins.push_back(b);
  }
  return bins;
}

} // namespace

TEST(the_tanh_fitted_crosses_the_limit_where_the_curve_does) {
  // A planted tanh, r 0.05 and s0 0.4: the fit recovers it, and the limit is
  // where the curve is 0.3, s0 + r atanh(0.4), to within the interpolation
  // between bins.
  const auto bins = planted_bins(0.05, 0.4, 500);
  const ResolutionEstimate e = estimate_resolution(bins);
  check::is_true(e.fitted, "fitted");
  check::close(e.r, 0.05, 1e-4, "r");
  check::close(e.s0, 0.4, 1e-4, "s0");
  const double expected = 1.0 / std::sqrt(0.4 + 0.05 * std::atanh(0.4));
  check::close(e.d_min_cc_half, expected, 2e-3, "the limit at CC half 0.3");
  // Significance: the bins switch from significant to not where the curve
  // crosses the critical value of about 0.06; the logistic finds it.
  const double critical = bins.front().critical;
  const double switch_at =
      1.0 / std::sqrt(0.4 + 0.05 * std::atanh(1.0 - 2.0 * critical));
  check::close(e.d_min_significance, switch_at, 0.02, "the significance limit");
}

TEST(a_curve_above_the_limit_everywhere_gives_the_last_bin) {
  const auto bins = planted_bins(0.05, 5.0, 500); // CC half near one throughout
  const ResolutionEstimate e = estimate_resolution(bins);
  check::close(e.d_min_cc_half, bins.back().d_min, 1e-12,
               "the highest resolution measured");
  check::close(e.d_min_significance, bins.back().d_min, 1e-12,
               "and all significant");
}

TEST(bins_with_no_pairs_do_not_decide_the_limit) {
  // A sweep with CC half near one to the edge, and two last bins of
  // reflections seen once -- a detector's corners -- with no pairs and CC half
  // reported as zero. The limit is the last bin with pairs, not none.
  auto bins = planted_bins(0.05, 5.0, 500);
  for (std::size_t k = bins.size() - 2; k < bins.size(); ++k) {
    bins[k].n = 0;
    bins[k].cc_half = 0.0;
    bins[k].significant = false;
  }
  const ResolutionEstimate e = estimate_resolution(bins);
  check::close(e.d_min_cc_half, bins[bins.size() - 3].d_min, 1e-12,
               "the last bin with pairs, not none");
}

} // namespace mxi
