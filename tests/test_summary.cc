#include <cmath>
#include <limits>

#include "../src/refl.hh"
#include "../src/spots/histogram.hh"
#include "../src/summary.hh"
#include "check.hh"

namespace mxi {

namespace {

//: Four reflections chosen so every number in the summary can be worked out by
//: hand: two at low resolution, two at high, with known flags and intensities.
SummaryInput hand_worked() {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  SummaryInput in;
  // d of 4 and 2: 1/d^3 is 1/64 and 1/8, so two shells split at their mean.
  in.d = {4.0, 4.0, 2.0, 2.0};
  in.flags = {
      flag::kPredicted | flag::kIntegratedSum | flag::kIntegratedPrf,
      flag::kPredicted | flag::kIntegratedSum,
      flag::kPredicted | flag::kIntegratedSum | flag::kIntegratedPrf,
      // a gap-crossing reflection: not summed, and says why
      flag::kPredicted | flag::kForegroundIncludesBadPixels |
          flag::kFailedDuringSummation | flag::kIntegratedPrf,
  };
  in.intensity_sum = {400.0, 90.0, 16.0, 5.0};
  in.variance_sum = {100.0, 9.0, 4.0, 1.0}; // I/sigma 40, 30, 8
  in.intensity_prf = {410.0, 0.0, 18.0, 30.0};
  in.variance_prf = {100.0, 1.0, 9.0, 25.0}; // I/sigma 41, -, 6, 6
  in.profile_correlation = {0.9, 0.0, 0.5, 0.7};
  in.background = {1.0, 3.0, 0.5, 9.0};
  in.partiality = {1.0, 0.5, 1.0, 1.0};
  // the second has no centre of mass: it must not count as a perfect one
  in.res_fast = {0.3, nan, 0.6, 0.0};
  in.res_slow = {0.4, nan, 0.8, 0.0};
  return in;
}

} // namespace

TEST(the_integration_summary_is_what_the_hand_says) {
  const IntegrationSummary s = summarise(hand_worked(), 2);
  check::equal(static_cast<long long>(s.shells.size()), 2, "two shells");

  // Highest resolution first, as dials.integrate shows them.
  const SummaryRow &high = s.shells.front();
  const SummaryRow &low = s.shells.back();
  check::close(high.d_min, 2.0, 1e-12, "the high shell is the d = 2 pair");
  check::close(low.d_min, 4.0, 1e-12, "the low shell is the d = 4 pair");

  check::equal(static_cast<long long>(low.n), 2, "two in the low shell");
  check::equal(static_cast<long long>(low.full), 1, "one fully recorded");
  check::equal(static_cast<long long>(low.partial), 1, "one partial");
  check::equal(static_cast<long long>(low.summed), 2, "both summed");
  check::equal(static_cast<long long>(low.fitted), 1, "one fitted");
  check::close(low.background, 2.0, 1e-12,
               "background over the summed: (1 + 3) / 2");
  check::close(low.isigma_sum, 35.0, 1e-12, "I/sigma sum: (40 + 30) / 2");
  check::close(low.isigma_prf, 41.0, 1e-12, "I/sigma prf over the one fitted");
  check::close(low.cc_prf, 0.9, 1e-12, "CC over the fitted only");
  // One centre, residual 0.5 px; the reflection without one is left out rather
  // than counted as a perfect prediction.
  check::close(low.rmsd_xy, 0.5, 1e-12, "RMSD over reflections with a centre");

  check::equal(static_cast<long long>(high.summed), 1,
               "the gap one is not summed");
  check::equal(static_cast<long long>(high.fitted), 2, "but both are fitted");
  check::equal(static_cast<long long>(high.bad_foreground), 1,
               "and it is counted");
  check::close(high.background, 0.5, 1e-12, "background only over the summed");
  check::close(high.isigma_prf, 6.0, 1e-12, "I/sigma prf: (6 + 6) / 2");
  check::close(high.rmsd_xy, std::sqrt((1.0 + 0.0) / 2.0), 1e-12,
               "RMSD: residuals 1.0 and 0.0 px");

  check::equal(static_cast<long long>(s.overall.n), 4, "four overall");
  check::equal(static_cast<long long>(s.overall.summed), 3,
               "three summed overall");
  check::equal(static_cast<long long>(s.both), 2,
               "two have both a sum and a fit, for the correlation");
}

TEST(a_summary_leaves_out_what_has_no_resolution) {
  SummaryInput in = hand_worked();
  in.d[1] = std::numeric_limits<double>::quiet_NaN();
  in.d[2] = 0.0;
  const IntegrationSummary s = summarise(in, 2);
  check::equal(static_cast<long long>(s.overall.n), 2,
               "a NaN or zero d is in no row");
}

TEST(the_spot_histogram_is_what_the_hand_says) {
  // Three images with 0, 5 and 10 spots, two rows high: the top row is
  // starred where a column exceeds half the largest, the bottom wherever there
  // are any spots at all.
  const std::vector<std::string> h =
      spots::spot_histogram({0, 5, 10}, 1, 60, 2);
  check::equal(static_cast<long long>(h.size()), 4,
               "a header, two rows, an axis");
  check::is_true(h[0] == "15 spots found on 3 images (max 10 / bin)", h[0]);
  check::is_true(h[1] == "  *", "top row: only the column at the maximum");
  check::is_true(h[2] == " **", "bottom row: every column with spots");
  check::is_true(h[3] == "1 3", "first image left, last right");
}

TEST(the_spot_histogram_sums_runs_of_images_into_a_column) {
  // A hundred images of one spot each, ten columns: each column is ten images
  // and holds ten spots, so every row is full.
  const std::vector<std::string> h =
      spots::spot_histogram(std::vector<std::size_t>(100, 1), 1, 10, 3);
  check::is_true(h[0] == "100 spots found on 100 images (max 10 / bin)", h[0]);
  for (std::size_t row = 1; row <= 3; ++row)
    check::is_true(h[row] == "**********", "a flat scan is a full block");
  check::is_true(h[4].find("1") == 0 && h[4].rfind("100") == 7, h[4]);
}

TEST(a_column_under_a_twentieth_of_the_tallest_shows_no_star) {
  // Heights are rounded, as dials.find_spots rounds them: a column of 1
  // against a tallest of 25 scales to 0.4 of a row and draws nothing, where
  // rounding up would have drawn it on the bottom row.
  const std::vector<std::string> h = spots::spot_histogram({1, 25}, 1, 60, 10);
  check::is_true(h[10] == " *", "the bottom row: the small column is blank");
  check::is_true(h[1] == " *", "the top row: the tallest reaches it");
}

} // namespace mxi
