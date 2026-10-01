#include "resolution.hh"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>

namespace mxi {

namespace {

//: The upper p quantile of the standard normal, by bisection on erfc: slow
//: and exact, which is what a handful of calls wants.
double normal_upper(double p) {
  double lo = -40.0, hi = 40.0;
  for (int it = 0; it < 200; ++it) {
    const double mid = 0.5 * (lo + hi);
    (0.5 * std::erfc(mid / std::sqrt(2.0)) > p ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

double pearson(const std::vector<double> &a, const std::vector<double> &b) {
  const double n = static_cast<double>(a.size());
  if (a.size() < 3)
    return 0.0;
  const double ma = std::accumulate(a.begin(), a.end(), 0.0) / n;
  const double mb = std::accumulate(b.begin(), b.end(), 0.0) / n;
  double sab = 0.0, saa = 0.0, sbb = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    sab += (a[i] - ma) * (b[i] - mb);
    saa += (a[i] - ma) * (a[i] - ma);
    sbb += (b[i] - mb) * (b[i] - mb);
  }
  return saa > 0.0 && sbb > 0.0 ? sab / std::sqrt(saa * sbb) : 0.0;
}

//: Levenberg-Marquardt on two parameters, residuals r_i(p) / sigma_i, with a
//: numerical Jacobian: enough for a curve through fifty points.
template <typename Model>
bool fit_two(Model model, const std::vector<double> &x,
             const std::vector<double> &y, const std::vector<double> &sigma,
             double p[2]) {
  const auto cost = [&](const double q[2]) {
    double c = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      const double r = (y[i] - model(x[i], q)) / sigma[i];
      c += r * r;
    }
    return c;
  };
  double c = cost(p), mu = 1e-3;
  for (int it = 0; it < 200; ++it) {
    double JtJ[2][2] = {{0, 0}, {0, 0}}, Jtr[2] = {0, 0};
    for (std::size_t i = 0; i < x.size(); ++i) {
      const double f = model(x[i], p);
      double j[2];
      for (int k = 0; k < 2; ++k) {
        double q[2] = {p[0], p[1]};
        const double h = 1e-7 * std::fmax(1.0, std::abs(p[k]));
        q[k] += h;
        j[k] = (model(x[i], q) - f) / h / sigma[i];
      }
      const double r = (y[i] - f) / sigma[i];
      for (int a = 0; a < 2; ++a) {
        Jtr[a] += j[a] * r;
        for (int b = 0; b < 2; ++b)
          JtJ[a][b] += j[a] * j[b];
      }
    }
    bool better = false;
    for (int tries = 0; tries < 20 && !better; ++tries) {
      const double a = JtJ[0][0] * (1 + mu), d = JtJ[1][1] * (1 + mu),
                   b = JtJ[0][1];
      const double det = a * d - b * b;
      if (!(std::abs(det) > 0.0)) {
        mu *= 10;
        continue;
      }
      const double q[2] = {p[0] + (d * Jtr[0] - b * Jtr[1]) / det,
                           p[1] + (a * Jtr[1] - b * Jtr[0]) / det};
      const double c2 = cost(q);
      if (c2 < c) {
        const bool done = c - c2 < 1e-12 * std::fmax(1.0, c);
        p[0] = q[0];
        p[1] = q[1];
        c = c2;
        mu = std::fmax(mu / 10, 1e-12);
        better = true;
        if (done)
          return true;
      } else {
        mu *= 10;
      }
    }
    if (!better)
      return true;
  }
  return true;
}

} // namespace

double student_t_quantile(double p_upper, double nu) {
  const double z = normal_upper(p_upper);
  const double z3 = z * z * z, z5 = z3 * z * z, z7 = z5 * z * z,
               z9 = z7 * z * z;
  return z + (z3 + z) / (4 * nu) + (5 * z5 + 16 * z3 + 3 * z) / (96 * nu * nu) +
         (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * nu * nu * nu) +
         (79 * z9 + 776 * z7 + 1482 * z5 - 1920 * z3 - 945 * z) /
             (92160 * nu * nu * nu * nu);
}

std::vector<ResolutionBin>
cc_half_bins(const ScaleData &data, const std::vector<double> &g,
             const Crystal &crystal, int nbins, std::size_t min_per_bin,
             double significance_level, std::size_t *wilson_outliers) {
  std::vector<std::vector<std::size_t>> members(data.stats_unique().size());
  for (std::size_t i = 0; i < data.size(); ++i)
    if (!data.outlier[i])
      members[data.stats_group()[i]].push_back(i);
  std::vector<double> d(data.stats_unique().size(), 0.0),
      mean(data.stats_unique().size(), 0.0);
  std::vector<std::size_t> present;
  for (std::size_t h = 0; h < members.size(); ++h) {
    if (members[h].empty())
      continue;
    const Miller &u = data.stats_unique()[h];
    d[h] = crystal.d_spacing(u[0], u[1], u[2]);
    double s = 0.0;
    for (std::size_t i : members[h])
      s += data.intensity[i] / g[i];
    mean[h] = s / static_cast<double>(members[h].size());
    present.push_back(h);
  }
  std::sort(present.begin(), present.end(),
            [&](std::size_t a, std::size_t b) { return d[a] > d[b]; });
  // Wilson outliers: E^2 = <I> over the mean <I> of its neighbours in d.
  std::size_t removed = 0;
  {
    std::vector<bool> keep(members.size(), true);
    const std::size_t shell = 500;
    for (std::size_t a = 0; a < present.size(); a += shell) {
      const std::size_t b = std::min(present.size(), a + shell);
      double s = 0.0;
      for (std::size_t k = a; k < b; ++k)
        s += mean[present[k]];
      const double m = s / static_cast<double>(b - a);
      for (std::size_t k = a; k < b; ++k)
        if (m > 0.0 && mean[present[k]] / m >= 16.0) {
          keep[present[k]] = false;
          ++removed;
        }
    }
    present.erase(std::remove_if(present.begin(), present.end(),
                                 [&](std::size_t h) { return !keep[h]; }),
                  present.end());
  }
  if (wilson_outliers)
    *wilson_outliers = removed;
  std::vector<ResolutionBin> out;
  if (present.empty())
    return out;
  const std::size_t per =
      std::max(min_per_bin, present.size() / static_cast<std::size_t>(nbins));
  std::mt19937 rng(20);
  for (std::size_t a = 0; a < present.size(); a += per) {
    std::size_t b = std::min(present.size(), a + per);
    if (present.size() - b < per / 2)
      b = present.size(); // the remainder joins the last bin
    ResolutionBin bin;
    bin.d_min = d[present[b - 1]];
    bin.d_star_sq = 1.0 / (bin.d_min * bin.d_min);
    std::vector<double> x, y;
    for (std::size_t k = a; k < b; ++k) {
      std::vector<std::size_t> o = members[present[k]];
      if (o.size() < 2)
        continue;
      std::shuffle(o.begin(), o.end(), rng);
      const std::size_t split = o.size() / 2;
      double s1 = 0.0, s2 = 0.0;
      for (std::size_t j = 0; j < o.size(); ++j)
        (j < split ? s1 : s2) += data.intensity[o[j]] / g[o[j]];
      x.push_back(s1 / static_cast<double>(split));
      y.push_back(s2 / static_cast<double>(o.size() - split));
    }
    bin.n = x.size();
    bin.cc_half = pearson(x, y);
    if (bin.n > 2) {
      const double t = student_t_quantile(significance_level,
                                          static_cast<double>(bin.n - 2));
      bin.critical = t / std::sqrt(static_cast<double>(bin.n) - 2.0 + t * t);
      bin.significant = bin.cc_half > bin.critical;
    }
    out.push_back(bin);
    if (b == present.size())
      break;
  }
  return out;
}

ResolutionEstimate estimate_resolution(const std::vector<ResolutionBin> &bins,
                                       double limit) {
  ResolutionEstimate out;
  if (bins.empty())
    return out;
  std::vector<double> x, y, sigma;
  std::size_t enough = 0;
  for (const ResolutionBin &b : bins) {
    x.push_back(b.d_star_sq);
    y.push_back(b.cc_half);
    sigma.push_back(b.n > 4 ? 1.0 / std::sqrt(static_cast<double>(b.n) - 3.0)
                            : 100.0);
    if (b.n > 3)
      ++enough;
  }
  const double highest = 1.0 / std::sqrt(x.back());
  // CC half: the tanh, and where its fitted values cross the limit.
  if (enough >= 3) {
    double p[2] = {0.2, 0.4};
    const auto tanh_model = [](double v, const double q[2]) {
      return 0.5 * (1.0 - std::tanh((v - q[1]) / q[0]));
    };
    fit_two(tanh_model, x, y, sigma, p);
    out.r = p[0];
    out.s0 = p[1];
    out.fitted = true;
    // Only bins with pairs enough to fit decide it -- the n > 3 the fit itself
    // asks. A bin of reflections seen once, as a detector's corners are, has no
    // halves and a CC half of zero, and under dials.estimate_resolution's rule
    // (every bin above the limit, or else the fitted values' crossing) it left
    // a sweep with CC half near one to its last bin with no limit at all.
    double min_informative = HUGE_VAL, last_informative = highest;
    for (const ResolutionBin &b : bins)
      if (b.n > 3) {
        min_informative = std::fmin(min_informative, b.cc_half);
        last_informative = b.d_min;
      }
    std::vector<double> f;
    for (double v : x)
      f.push_back(tanh_model(v, p));
    bool crossed = false;
    if (min_informative <= limit) {
      for (std::size_t j = 1; j < x.size(); ++j)
        if ((f[j - 1] - limit) * (f[j] - limit) < 0.0) {
          const double s = x[j - 1] + (limit - f[j - 1]) * (x[j] - x[j - 1]) /
                                          (f[j] - f[j - 1]);
          out.d_min_cc_half = 1.0 / std::sqrt(s);
          crossed = true;
          out.reached = true;
          break;
        }
    }
    // Above the limit in every bin that says anything, or a fit that never
    // comes down to it: the data reach it nowhere, and the limit is the last
    // bin with pairs.
    if (!crossed && (min_informative > limit || f.back() > limit))
      out.d_min_cc_half = last_informative;
  }
  // Significance: a logistic through which bins are significant.
  bool all = true, none = true;
  for (const ResolutionBin &b : bins) {
    all = all && b.significant;
    none = none && !b.significant;
  }
  if (all) {
    out.d_min_significance = highest;
  } else if (!none) {
    double start = x.front();
    for (std::size_t j = 0; j < bins.size(); ++j)
      if (!bins[j].significant)
        start = x[j]; // the highest d* at which it is not significant
    std::vector<double> flags, ones(bins.size(), 1.0);
    for (const ResolutionBin &b : bins)
      flags.push_back(b.significant ? 1.0 : 0.0);
    double p[2] = {100.0, start};
    const auto logistic = [](double v, const double q[2]) {
      return 1.0 - 1.0 / (1.0 + std::exp(-q[0] * (v - q[1])));
    };
    fit_two(logistic, x, flags, ones, p);
    if (p[0] > 0.0 && p[1] > 0.0)
      out.d_min_significance = 1.0 / std::sqrt(p[1]);
  }
  return out;
}

} // namespace mxi
