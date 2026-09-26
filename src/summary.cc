#include "summary.hh"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "refl.hh"

namespace mxi {

namespace {

//: Sums that become a row's means once every reflection has been seen.
struct Accumulator {
  SummaryRow row;
  double background = 0.0, isigma_sum = 0.0, isigma_prf = 0.0, cc = 0.0;
  double residual2 = 0.0;
  std::size_t n_isigma_sum = 0, n_isigma_prf = 0, n_cc = 0, n_residual = 0;
  double d_min = HUGE_VAL, d_max = 0.0;

  void add(const SummaryInput &in, std::size_t i) {
    const std::int64_t f = in.flags[i];
    const double d = in.d[i];
    ++row.n;
    d_min = std::min(d_min, d);
    d_max = std::max(d_max, d);
    if (in.partiality[i] >= 0.99)
      ++row.full;
    else
      ++row.partial;
    if (f & flag::kForegroundIncludesBadPixels)
      ++row.bad_foreground;
    if (f & flag::kBackgroundIncludesBadPixels)
      ++row.bad_background;
    if (f & flag::kIntegratedSum) {
      ++row.summed;
      background += in.background[i];
      if (in.variance_sum[i] > 0.0) {
        isigma_sum += in.intensity_sum[i] / std::sqrt(in.variance_sum[i]);
        ++n_isigma_sum;
      }
    }
    if (f & flag::kIntegratedPrf) {
      ++row.fitted;
      if (in.variance_prf[i] > 0.0) {
        isigma_prf += in.intensity_prf[i] / std::sqrt(in.variance_prf[i]);
        ++n_isigma_prf;
      }
      if (std::isfinite(in.profile_correlation[i])) {
        cc += in.profile_correlation[i];
        ++n_cc;
      }
    }
    if (std::isfinite(in.res_fast[i]) && std::isfinite(in.res_slow[i])) {
      residual2 +=
          in.res_fast[i] * in.res_fast[i] + in.res_slow[i] * in.res_slow[i];
      ++n_residual;
    }
  }

  SummaryRow finish() const {
    SummaryRow out = row;
    out.d_min = row.n ? d_min : 0.0;
    out.d_max = row.n ? d_max : 0.0;
    const auto mean = [](double sum, std::size_t n) {
      return n ? sum / static_cast<double>(n) : 0.0;
    };
    out.background = mean(background, row.summed);
    out.isigma_sum = mean(isigma_sum, n_isigma_sum);
    out.isigma_prf = mean(isigma_prf, n_isigma_prf);
    out.cc_prf = mean(cc, n_cc);
    out.rmsd_xy = std::sqrt(mean(residual2, n_residual));
    return out;
  }
};

double pearson(const std::vector<double> &a, const std::vector<double> &b) {
  const std::size_t n = a.size();
  if (n < 2)
    return 0.0;
  const double ma =
      std::accumulate(a.begin(), a.end(), 0.0) / static_cast<double>(n);
  const double mb =
      std::accumulate(b.begin(), b.end(), 0.0) / static_cast<double>(n);
  double sab = 0.0, saa = 0.0, sbb = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    sab += (a[i] - ma) * (b[i] - mb);
    saa += (a[i] - ma) * (a[i] - ma);
    sbb += (b[i] - mb) * (b[i] - mb);
  }
  return saa > 0.0 && sbb > 0.0 ? sab / std::sqrt(saa * sbb) : 0.0;
}

//: Ranks, ties sharing the mean of the ranks they span.
std::vector<double> ranks(const std::vector<double> &v) {
  std::vector<std::size_t> order(v.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(),
            [&](std::size_t x, std::size_t y) { return v[x] < v[y]; });
  std::vector<double> out(v.size());
  for (std::size_t i = 0; i < order.size();) {
    std::size_t j = i;
    while (j + 1 < order.size() && v[order[j + 1]] == v[order[i]])
      ++j;
    const double rank = 0.5 * static_cast<double>(i + j);
    for (std::size_t k = i; k <= j; ++k)
      out[order[k]] = rank;
    i = j + 1;
  }
  return out;
}

} // namespace

IntegrationSummary summarise(const SummaryInput &in, int shells) {
  IntegrationSummary out;
  const std::size_t n = in.d.size();
  shells = std::max(shells, 1);

  // Equal volume in reciprocal space: equal steps in 1/d^3.
  double lo = HUGE_VAL, hi = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    if (!(in.d[i] > 0.0) || !std::isfinite(in.d[i]))
      continue;
    const double v = 1.0 / (in.d[i] * in.d[i] * in.d[i]);
    lo = std::min(lo, v);
    hi = std::max(hi, v);
  }
  if (!(hi > 0.0))
    return out;

  std::vector<Accumulator> bins(static_cast<std::size_t>(shells));
  Accumulator overall;
  std::vector<double> both_sum, both_prf;
  for (std::size_t i = 0; i < n; ++i) {
    if (!(in.d[i] > 0.0) || !std::isfinite(in.d[i]))
      continue;
    const double v = 1.0 / (in.d[i] * in.d[i] * in.d[i]);
    // Index 0 is the lowest resolution; the highest-resolution reflection
    // belongs in the last shell rather than one past it.
    std::size_t b = hi > lo
                        ? static_cast<std::size_t>((v - lo) / (hi - lo) *
                                                   static_cast<double>(shells))
                        : 0;
    b = std::min(b, static_cast<std::size_t>(shells) - 1);
    bins[b].add(in, i);
    overall.add(in, i);
    const std::int64_t f = in.flags[i];
    if ((f & flag::kIntegratedSum) && (f & flag::kIntegratedPrf)) {
      both_sum.push_back(in.intensity_sum[i]);
      both_prf.push_back(in.intensity_prf[i]);
    }
  }
  for (auto it = bins.rbegin(); it != bins.rend(); ++it) {
    if (it->row.n > 0)
      out.shells.push_back(it->finish());
  }
  out.overall = overall.finish();
  out.both = both_sum.size();
  out.pearson_sum_prf = pearson(both_sum, both_prf);
  out.spearman_sum_prf = pearson(ranks(both_sum), ranks(both_prf));
  return out;
}

void print_summary(std::FILE *out, const IntegrationSummary &s) {
  if (s.shells.empty()) {
    std::fprintf(out, "\nno reflections with a resolution to summarise\n");
    return;
  }
  std::fprintf(out, "\nSummary against resolution\n");
  std::fprintf(out, "  %6s %7s %7s %7s %7s %7s %8s %8s %6s %7s\n", "d min",
               "# full", "# part", "# sum", "# prf", "Ibg", "I/sigI", "I/sigI",
               "CC", "RMSD XY");
  std::fprintf(out, "  %6s %7s %7s %7s %7s %7s %8s %8s %6s %7s\n", "(A)", "",
               "", "", "", "", "(sum)", "(prf)", "prf", "(px)");
  for (const SummaryRow &r : s.shells) {
    std::fprintf(out,
                 "  %6.2f %7zu %7zu %7zu %7zu %7.2f %8.2f %8.2f %6.3f %7.2f\n",
                 r.d_min, r.full, r.partial, r.summed, r.fitted, r.background,
                 r.isigma_sum, r.isigma_prf, r.cc_prf, r.rmsd_xy);
  }

  const SummaryRow &o = s.overall;
  const SummaryRow &low = s.shells.back();
  const SummaryRow &high = s.shells.front();
  std::fprintf(out, "\nSummary\n");
  std::fprintf(out, "  %-44s %9s %9s %9s\n", "", "overall", "low", "high");
  const auto count = [&](const char *what, std::size_t SummaryRow::*m) {
    std::fprintf(out, "  %-44s %9zu %9zu %9zu\n", what, o.*m, low.*m, high.*m);
  };
  const auto value = [&](const char *what, double SummaryRow::*m,
                         const char *fmt) {
    char a[32], b[32], c[32];
    std::snprintf(a, sizeof a, fmt, o.*m);
    std::snprintf(b, sizeof b, fmt, low.*m);
    std::snprintf(c, sizeof c, fmt, high.*m);
    std::fprintf(out, "  %-44s %9s %9s %9s\n", what, a, b, c);
  };
  value("d min (A)", &SummaryRow::d_min, "%.2f");
  value("d max (A)", &SummaryRow::d_max, "%.2f");
  count("reflections", &SummaryRow::n);
  count("  fully recorded", &SummaryRow::full);
  count("  partially recorded", &SummaryRow::partial);
  count("integrated by summation", &SummaryRow::summed);
  count("integrated by profile fitting", &SummaryRow::fitted);
  count("with the foreground reaching a masked pixel",
        &SummaryRow::bad_foreground);
  count("with the background reaching a masked pixel",
        &SummaryRow::bad_background);
  value("mean background (counts/pixel)", &SummaryRow::background, "%.3f");
  value("mean I/sigma (summation)", &SummaryRow::isigma_sum, "%.2f");
  value("mean I/sigma (profile fitting)", &SummaryRow::isigma_prf, "%.2f");
  value("mean profile correlation", &SummaryRow::cc_prf, "%.3f");
  value("RMSD of the centre from the prediction (px)", &SummaryRow::rmsd_xy,
        "%.3f");
  std::fprintf(out,
               "  summed against fitted intensity, over %zu reflections: "
               "Pearson %.4f, Spearman %.4f\n",
               s.both, s.pearson_sum_prf, s.spearman_sum_prf);
}

} // namespace mxi
