#include "laue.hh"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <random>

namespace mxi {

namespace {

const Rotation kIdentityRotation{1, 0, 0, 0, 1, 0, 0, 0, 1};

Rotation times(const Rotation &a, const Rotation &b) {
  Rotation out{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      int v = 0;
      for (int k = 0; k < 3; ++k)
        v += a[static_cast<std::size_t>(i * 3 + k)] *
             b[static_cast<std::size_t>(k * 3 + j)];
      out[static_cast<std::size_t>(i * 3 + j)] = v;
    }
  return out;
}

//: The P -1 representative: of h and -h, the one first in order.
Miller p1_unique(const Miller &h) {
  const Miller n{-h[0], -h[1], -h[2]};
  return n < h ? n : h;
}

double cauchy_cdf(double x, double loc, double scale) {
  return 0.5 + std::atan((x - loc) / scale) / std::acos(-1.0);
}

double trunc_cauchy_pdf(double x, double a, double b, double loc,
                        double scale) {
  const double pdf =
      scale / (std::acos(-1.0) * ((x - loc) * (x - loc) + scale * scale));
  return pdf / (cauchy_cdf(b, loc, scale) - cauchy_cdf(a, loc, scale));
}

struct Accumulator {
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  std::size_t n = 0;
  void add(double x, double y) {
    sx += x;
    sy += y;
    sxx += x * x;
    syy += y * y;
    sxy += x * y;
    ++n;
  }
  double cc() const {
    if (n < 2)
      return 0.0;
    const double N = static_cast<double>(n);
    const double v = (sxx - sx * sx / N) * (syy - sy * sy / N);
    return v > 0.0 ? (sxy - sx * sy / N) / std::sqrt(v) : 0.0;
  }
};

//: Quartiles by linear interpolation.
double quantile_of(std::vector<double> v, double q) {
  std::sort(v.begin(), v.end());
  const double at = q * static_cast<double>(v.size() - 1);
  const auto lo = static_cast<std::size_t>(at);
  const std::size_t hi = std::min(lo + 1, v.size() - 1);
  return v[lo] + (at - static_cast<double>(lo)) * (v[hi] - v[lo]);
}

} // namespace

int rotation_order(const Rotation &r) {
  Rotation p = r;
  for (int n = 1; n <= 6; ++n) {
    if (p == kIdentityRotation)
      return n;
    p = times(p, r);
  }
  return 0;
}

std::vector<Rotation> symmetry_elements(const std::vector<Rotation> &lattice) {
  std::vector<Rotation> out;
  for (const Rotation &r : lattice) {
    if (std::find(out.begin(), out.end(), r) != out.end())
      continue;
    const int n = rotation_order(r);
    // A rotation of order above two stands with its inverse, r^(n-1).
    bool inverse_listed = false;
    if (n > 2) {
      Rotation inv = r;
      for (int k = 1; k < n - 1; ++k)
        inv = times(inv, r);
      inverse_listed = std::find(out.begin(), out.end(), inv) != out.end();
    }
    if (!inverse_listed)
      out.push_back(r);
  }
  std::stable_sort(out.begin(), out.end(),
                   [](const Rotation &a, const Rotation &b) {
                     return rotation_order(a) == 1 && rotation_order(b) != 1;
                   });
  return out;
}

void normalise(P1Intensities &data, std::size_t per_shell) {
  std::vector<std::size_t> order(data.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return data.d[a] > data.d[b];
  });
  for (std::size_t start = 0; start < order.size(); start += per_shell) {
    std::size_t end = std::min(order.size(), start + per_shell);
    if (order.size() - end < per_shell / 2)
      end = order.size(); // no shell of a handful at the end
    double mean = 0.0;
    for (std::size_t k = start; k < end; ++k)
      mean += data.i[order[k]] / static_cast<double>(end - start);
    if (mean > 0.0)
      for (std::size_t k = start; k < end; ++k) {
        data.i[order[k]] /= mean;
        data.sigma[order[k]] /= mean;
      }
    if (end == order.size())
      break;
  }
}

double p_cc_given_present(double cc, double sigma_cc, double expected) {
  return trunc_cauchy_pdf(cc, -1.0, 1.0, expected, sigma_cc);
}

double p_cc_given_absent(double cc, double sigma_cc) {
  // integral over x in [0, 1] of pdf(cc; x) (1 - x^2)^(1/2), over that of
  // (1 - x^2)^(1/2), which is pi/4; by the midpoint rule on 2000 steps.
  const int steps = 2000;
  double sum = 0.0;
  for (int k = 0; k < steps; ++k) {
    const double x = (k + 0.5) / steps;
    sum +=
        trunc_cauchy_pdf(cc, -1.0, 1.0, x, sigma_cc) * std::sqrt(1.0 - x * x);
  }
  return (sum / steps) / (std::acos(-1.0) / 4.0);
}

LaueScores score_laue_groups(const P1Intensities &data,
                             const std::vector<Rotation> &lattice,
                             unsigned seed) {
  LaueScores out;
  std::map<Miller, std::size_t> where;
  for (std::size_t k = 0; k < data.size(); ++k)
    where[p1_unique(data.hkl[k])] = k;

  // sigma(CC) as a function of sample size, from pairs of unrelated
  // reflections at similar resolution: rms CC against 1/sqrt(n), fitted.
  std::mt19937 rng(seed);
  {
    std::vector<std::size_t> by_d(data.size());
    std::iota(by_d.begin(), by_d.end(), std::size_t{0});
    std::sort(by_d.begin(), by_d.end(), [&](std::size_t a, std::size_t b) {
      return data.d[a] > data.d[b];
    });
    std::vector<double> a, b;
    const std::size_t bin = std::max<std::size_t>(200, data.size() / 500 + 1);
    for (std::size_t s = 0; s + 1 < by_d.size(); s += bin) {
      std::vector<std::size_t> in(
          by_d.begin() + static_cast<long>(s),
          by_d.begin() + static_cast<long>(std::min(by_d.size(), s + bin)));
      std::shuffle(in.begin(), in.end(), rng);
      for (std::size_t k = 0; k + 1 < in.size(); k += 2) {
        a.push_back(data.i[in[k]]);
        b.push_back(data.i[in[k + 1]]);
      }
    }
    const std::size_t pairs = a.size();
    const std::size_t max_n = std::min<std::size_t>(pairs / 10, 200),
                      min_n = std::min<std::size_t>(5, max_n);
    if (max_n >= min_n + 4) {
      double sx = 0, sy = 0, sxx = 0, sxy = 0, k = 0;
      std::uniform_int_distribution<std::size_t> pick(0, pairs - 1);
      for (std::size_t n = min_n; n <= max_n; ++n) {
        double ms = 0.0;
        for (int rep = 0; rep < 200; ++rep) {
          Accumulator acc;
          for (std::size_t j = 0; j < n; ++j) {
            const std::size_t p = pick(rng);
            acc.add(a[p], b[p]);
          }
          ms += acc.cc() * acc.cc() / 200.0;
        }
        const double x = 1.0 / std::sqrt(static_cast<double>(n)),
                     y = std::sqrt(ms);
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
        k += 1.0;
      }
      out.cc_sig_fac = (sxy - sx * sy / k) / (sxx - sx * sx / k);
    }
  }
  // E(CC; S): the intensities' variance over it plus the mean variance of the
  // errors, combined as dials.symmetry does with CC of the identity, which on
  // intensities with Friedel mates merged is one.
  {
    double mean = 0.0, s2 = 0.0;
    for (std::size_t k = 0; k < data.size(); ++k) {
      mean += data.i[k] / static_cast<double>(data.size());
      s2 += data.sigma[k] * data.sigma[k] / static_cast<double>(data.size());
    }
    double var = 0.0;
    for (std::size_t k = 0; k < data.size(); ++k)
      var += (data.i[k] - mean) * (data.i[k] - mean) /
             static_cast<double>(data.size() - 1);
    out.e_cc_true = var / (var + s2);
    out.cc_identity = 1.0;
    const double sigma_1 = std::fmax(0.05, out.cc_sig_fac / std::sqrt(200.0));
    const double sigma_2 = std::fmax(
        0.05, out.cc_sig_fac / std::sqrt(static_cast<double>(data.size())));
    const double w1 = sigma_1 > 1e-4 ? 1.0 / (sigma_1 * sigma_1) : 0.0;
    const double w2 = data.size() > 10 ? 1.0 / (sigma_2 * sigma_2) : 0.0;
    out.cc_true = (w1 * out.e_cc_true + w2 * out.cc_identity) / (w1 + w2);
  }
  // Each element: CC over the reflections it relates, excluding those it
  // leaves in place (on the axis), with dials.symmetry's generous outlier cut.
  for (const Rotation &r : symmetry_elements(lattice)) {
    ElementScore e;
    e.rotation = r;
    e.order = rotation_order(r);
    std::vector<Rotation> ops{r};
    if (e.order > 2) {
      Rotation inv = r;
      for (int k = 1; k < e.order - 1; ++k)
        inv = times(inv, r);
      ops.push_back(inv);
    }
    std::vector<double> xs, ys;
    for (const Rotation &op : ops) {
      for (std::size_t k = 0; k < data.size(); ++k) {
        const Miller h = data.hkl[k];
        const Miller image = p1_unique(
            e.order == 1 ? Miller{-h[0], -h[1], -h[2]} : apply(op, h));
        if (e.order > 1 && image == p1_unique(h))
          continue; // on the axis: epsilon > 1
        const auto it = where.find(image);
        if (it == where.end())
          continue;
        xs.push_back(data.i[k]);
        ys.push_back(data.i[it->second]);
      }
    }
    std::vector<bool> keep(xs.size(), true);
    for (const std::vector<double> *col : {&xs, &ys}) {
      if (col->empty())
        break;
      const double q1 = quantile_of(*col, 0.25), q3 = quantile_of(*col, 0.75);
      const double cut = 20.0 * (q3 - q1);
      for (std::size_t k = 0; k < col->size(); ++k)
        if ((*col)[k] > q3 + cut || (*col)[k] < q1 - cut)
          keep[k] = false;
    }
    Accumulator acc;
    for (std::size_t k = 0; k < xs.size(); ++k)
      if (keep[k])
        acc.add(xs[k], ys[k]);
    e.pairs = acc.n;
    e.sx = acc.sx;
    e.sy = acc.sy;
    e.sxx = acc.sxx;
    e.syy = acc.syy;
    e.sxy = acc.sxy;
    if (e.pairs > 0) {
      e.cc = acc.cc();
      e.sigma_cc = std::fmax(0.1, out.cc_sig_fac /
                                      std::sqrt(static_cast<double>(e.pairs)));
      e.z = e.cc / e.sigma_cc;
      e.p_given_present = p_cc_given_present(e.cc, e.sigma_cc, out.cc_true);
      e.p_given_absent = p_cc_given_absent(e.cc, e.sigma_cc);
      e.likelihood = e.p_given_present / (e.p_given_present + e.p_given_absent);
    }
    out.elements.push_back(e);
  }
  // Each subgroup: Evans (2011) A2.
  double total = 0.0;
  for (const std::vector<Rotation> &g : subgroups(lattice)) {
    GroupScore s;
    s.rotations = g;
    Accumulator in, out_of;
    double pl = 0.0, zf = 0.0, za = 0.0;
    std::size_t nf = 0, na = 0;
    for (const ElementScore &e : out.elements) {
      const bool has = std::find(g.begin(), g.end(), e.rotation) != g.end();
      s.contains.push_back(has);
      Accumulator &acc = has ? in : out_of;
      acc.sx += e.sx;
      acc.sy += e.sy;
      acc.sxx += e.sxx;
      acc.syy += e.syy;
      acc.sxy += e.sxy;
      acc.n += e.pairs;
      if (e.pairs <= 2)
        continue;
      if (has) {
        zf += e.z * e.z;
        ++nf;
        pl += std::log(e.p_given_present);
      } else {
        za += e.z * e.z;
        ++na;
        pl += std::log(e.p_given_absent);
      }
    }
    s.likelihood = std::exp(pl);
    s.z_for = nf ? std::sqrt(zf / nf) : 0.0;
    s.z_against = na ? std::sqrt(za / na) : 0.0;
    s.z_net = s.z_for - s.z_against;
    s.cc_for = in.cc();
    s.cc_against = out_of.cc();
    total += s.likelihood;
    out.groups.push_back(std::move(s));
  }
  for (GroupScore &s : out.groups)
    s.likelihood = total > 0.0 ? s.likelihood / total : 0.0;
  std::stable_sort(out.groups.begin(), out.groups.end(),
                   [](const GroupScore &a, const GroupScore &b) {
                     return a.likelihood > b.likelihood;
                   });
  return out;
}

} // namespace mxi
