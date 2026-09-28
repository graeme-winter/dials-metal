#include "scale.hh"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>

namespace mxi {

ScaleData build_scale_data(const ExperimentList &experiments,
                           const Table &reflections, const SpaceGroup &group,
                           const ScaleModelShape &shape,
                           const ScaleDataOptions &options) {
  if (experiments.size() == 0)
    throw std::runtime_error("scaling needs an experiment");
  const Experiment &e = experiments[0];
  const bool profile = options.use == IntensityChoice::profile;
  const Column &miller = reflections.at("miller_index");
  const Column &flags = reflections.at("flags");
  const Column &value =
      reflections.at(profile ? "intensity.prf.value" : "intensity.sum.value");
  const Column &var = reflections.at(profile ? "intensity.prf.variance"
                                             : "intensity.sum.variance");
  const Column &dcol = reflections.at("d");
  const Column &cal = reflections.at("xyzcal.px");
  const bool has_lp = reflections.has("lp"), has_qe = reflections.has("qe"),
             has_part = reflections.has("partiality"),
             has_s1 = reflections.has("s1"), has_id = reflections.has("id");
  const std::int64_t wanted =
      profile ? flag::kIntegratedPrf : flag::kIntegratedSum;
  const double images = static_cast<double>(e.scan.num_images());

  ScaleData data;
  std::map<Miller, std::size_t> groups;
  std::vector<double> harmonics_s1, harmonics_s0;
  for (std::size_t i = 0; i < reflections.nrows; ++i) {
    if ((flags.ints[i] & wanted) == 0)
      continue;
    if (has_id && reflections.at("id").ints[i] != 0)
      continue; // one sweep, for now
    const double v = var.reals[i];
    const double d = dcol.reals[i];
    if (!(v > 0.0) || !std::isfinite(value.reals[i]) || !(d > 0.0) ||
        !std::isfinite(d) || d < options.d_min)
      continue;
    const double part = has_part ? reflections.at("partiality").reals[i] : 1.0;
    if (!(part >= options.partiality_cutoff))
      continue;
    const Miller h{static_cast<int>(miller.ints[i * 3]),
                   static_cast<int>(miller.ints[i * 3 + 1]),
                   static_cast<int>(miller.ints[i * 3 + 2])};
    if ((h[0] == 0 && h[1] == 0 && h[2] == 0) || group.absent(h))
      continue;
    double factor = 1.0 / part;
    if (has_lp)
      factor *= reflections.at("lp").reals[i];
    if (has_qe && reflections.at("qe").reals[i] > 0.0)
      factor /= reflections.at("qe").reals[i];

    const Miller u = group.unique(h);
    const auto [it, fresh] = groups.emplace(u, data.unique.size());
    if (fresh)
      data.unique.push_back(u);

    ScaleObservation o;
    const double z = cal.reals[i * 3 + 2];
    o.rotation = images > 0.0 ? z / images : 0.0;
    o.time = o.rotation;
    o.inv_2d2 = 1.0 / (2.0 * d * d);
    if (shape.lmax > 0 && has_s1) {
      // The crystal frame: undo the goniometer's rotation at this angle.
      const Mat3 r = e.goniometer.setting *
                     rotation(e.goniometer.axis, e.scan.phi_from_z(z)) *
                     e.goniometer.fixed;
      const Mat3 back = r.transpose();
      const Column &s1 = reflections.at("s1");
      const Vec3 s1c = (back * Vec3{s1.reals[i * 3], s1.reals[i * 3 + 1],
                                    s1.reals[i * 3 + 2]})
                           .normalized();
      const Vec3 s0c =
          (back * e.beam.direction).normalized(); // the reverse incident beam
      real_spherical_harmonics(shape.lmax, s1c, &harmonics_s1);
      real_spherical_harmonics(shape.lmax, s0c, &harmonics_s0);
      o.absorption.resize(harmonics_s1.size());
      for (std::size_t k = 0; k < harmonics_s1.size(); ++k)
        o.absorption[k] = 0.5 * (harmonics_s1[k] + harmonics_s0[k]);
    }
    data.intensity.push_back(value.reals[i] * factor);
    data.variance.push_back(v * factor * factor);
    data.variance_before.push_back(v * factor * factor);
    data.scale_term.push_back(0.0);
    const Column &pv = reflections.at("intensity.prf.value");
    const Column &pvar = reflections.at("intensity.prf.variance");
    const Column &sv = reflections.at("intensity.sum.value");
    const Column &svar = reflections.at("intensity.sum.variance");
    const bool summed =
        (flags.ints[i] & flag::kIntegratedSum) != 0 && svar.reals[i] > 0.0;
    data.prf.push_back(pv.reals[i] * factor);
    data.prf_variance.push_back(pvar.reals[i] * factor * factor);
    data.sum.push_back(sv.reals[i] * factor);
    data.sum_variance.push_back(svar.reals[i] * factor * factor);
    data.has_sum.push_back(summed);
    data.observation.push_back(std::move(o));
    data.group.push_back(it->second);
    data.row.push_back(i);
    data.d.push_back(d);
    data.outlier.push_back(false);
  }
  return data;
}

std::vector<double> inverse_scales(const ScaleModel &model,
                                   const ScaleData &data) {
  std::vector<double> g(data.size());
  for (std::size_t i = 0; i < data.size(); ++i)
    g[i] = model.inverse_scale(data.observation[i]);
  return g;
}

std::vector<double> merged_intensities(const ScaleData &data,
                                       const std::vector<double> &g) {
  std::vector<double> num(data.unique.size(), 0.0),
      den(data.unique.size(), 0.0);
  for (std::size_t i = 0; i < data.size(); ++i) {
    if (data.outlier[i])
      continue;
    const double w = 1.0 / data.variance[i];
    num[data.group[i]] += w * g[i] * data.intensity[i];
    den[data.group[i]] += w * g[i] * g[i];
  }
  for (std::size_t h = 0; h < num.size(); ++h)
    num[h] = den[h] > 0.0 ? num[h] / den[h] : 0.0;
  return num;
}

std::vector<std::size_t> select_for_fitting(const ScaleData &data,
                                            std::size_t min_groups,
                                            std::size_t min_observations,
                                            unsigned seed) {
  std::vector<std::vector<std::size_t>> members(data.unique.size());
  for (std::size_t i = 0; i < data.size(); ++i)
    if (!data.outlier[i])
      members[data.group[i]].push_back(i);
  std::vector<std::size_t> order(members.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::mt19937 rng(seed);
  std::shuffle(order.begin(), order.end(), rng);
  std::vector<std::size_t> out;
  std::size_t groups = 0;
  for (std::size_t h : order) {
    if (groups >= min_groups && out.size() >= min_observations)
      break;
    if (members[h].size() < 2)
      continue; // a single observation says nothing about scale
    out.insert(out.end(), members[h].begin(), members[h].end());
    ++groups;
  }
  std::sort(out.begin(), out.end());
  return out;
}

namespace {

//: The target over `use`, with <I_h> from the same observations, plus the
//: restraints.
double target(const ScaleModel &model, const ScaleData &data,
              const std::vector<std::size_t> &use,
              const ScaleFitOptions &options, std::vector<double> *g_out,
              std::vector<double> *merged_out) {
  std::vector<double> g(data.size(), 0.0);
  for (std::size_t i : use)
    g[i] = model.inverse_scale(data.observation[i]);
  // Merged over `use` only, which is what is being fitted.
  std::vector<double> num(data.unique.size(), 0.0),
      den(data.unique.size(), 0.0);
  for (std::size_t i : use) {
    const double w = 1.0 / data.variance[i];
    num[data.group[i]] += w * g[i] * data.intensity[i];
    den[data.group[i]] += w * g[i] * g[i];
  }
  for (std::size_t h = 0; h < num.size(); ++h)
    num[h] = den[h] > 0.0 ? num[h] / den[h] : 0.0;
  double phi = 0.0;
  for (std::size_t i : use) {
    const double r = data.intensity[i] - g[i] * num[data.group[i]];
    phi += r * r / data.variance[i];
  }
  const ScaleModelShape &s = model.shape();
  for (std::size_t k = 0; k < s.decay_points; ++k)
    phi += options.decay_restraint * model.parameters[model.first_decay() + k] *
           model.parameters[model.first_decay() + k];
  for (std::size_t k = 0; k < harmonic_count(s.lmax); ++k)
    phi += options.absorption_restraint *
           model.parameters[model.first_absorption() + k] *
           model.parameters[model.first_absorption() + k];
  if (g_out)
    *g_out = std::move(g);
  if (merged_out)
    *merged_out = std::move(num);
  return phi;
}

} // namespace

namespace {

//: The normal equations of the fit at the model's current parameters, over the
//: groups in `members`: J^T J and -J^T r, with J the variable-projection
//: Jacobian of the residuals sqrt(w)(I - g <I>), and the restraints. `g` and
//: `merged` are the model's scales and merged intensities over the same
//: observations. The fit steps with these; the covariance inverts them.
void normal_equations(const ScaleModel &model, const ScaleData &data,
                      const std::vector<std::vector<std::size_t>> &members,
                      const std::vector<double> &g,
                      const std::vector<double> &merged,
                      const ScaleFitOptions &options, std::vector<double> &N,
                      std::vector<double> &b) {
  const std::size_t n = model.size();
  N.assign(n * n, 0.0);
  b.assign(n, 0.0);
  std::vector<std::pair<std::size_t, double>> grad;
  std::vector<double> G(n, 0.0), row(n, 0.0);
  std::vector<std::vector<std::pair<std::size_t, double>>> grads;
  std::vector<std::size_t> touched;
  for (const std::vector<std::size_t> &obs : members) {
    if (obs.empty())
      continue;
    const std::size_t h = data.group[obs.front()];
    const double m = merged[h];
    grads.resize(obs.size());
    std::fill(G.begin(), G.end(), 0.0);
    touched.clear();
    double den = 0.0;
    for (std::size_t k = 0; k < obs.size(); ++k) {
      const std::size_t i = obs[k];
      const double w = 1.0 / data.variance[i];
      model.inverse_scale(data.observation[i], &grads[k]);
      den += w * g[i] * g[i];
      for (const auto &[a, da] : grads[k]) {
        if (G[a] == 0.0)
          touched.push_back(a);
        G[a] += w * g[i] * da;
      }
    }
    if (!(den > 0.0))
      continue;
    std::sort(touched.begin(), touched.end());
    touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
    for (std::size_t a : touched)
      G[a] /= den;
    for (std::size_t k = 0; k < obs.size(); ++k) {
      const std::size_t i = obs[k];
      const double w = 1.0 / data.variance[i];
      for (std::size_t a : touched)
        row[a] = -g[i] * G[a];
      for (const auto &[a, da] : grads[k])
        row[a] += da;
      const double r = data.intensity[i] - g[i] * m;
      // J = -sqrt(w) m row; accumulate J^T J and -J^T r.
      for (std::size_t a : touched) {
        if (row[a] == 0.0)
          continue;
        b[a] += w * m * row[a] * r;
        for (std::size_t c : touched)
          N[a * n + c] += w * m * m * row[a] * row[c];
      }
      for (std::size_t a : touched)
        row[a] = 0.0;
    }
    for (std::size_t a : touched)
      G[a] = 0.0;
  }
  const ScaleModelShape &s = model.shape();
  for (std::size_t k = 0; k < s.decay_points; ++k) {
    const std::size_t a = model.first_decay() + k;
    N[a * n + a] += options.decay_restraint;
    b[a] -= options.decay_restraint * model.parameters[a];
  }
  for (std::size_t k = 0; k < harmonic_count(s.lmax); ++k) {
    const std::size_t a = model.first_absorption() + k;
    N[a * n + a] += options.absorption_restraint;
    b[a] -= options.absorption_restraint * model.parameters[a];
  }
}

} // namespace

ScaleFitResult fit_scale_model(ScaleModel &model, const ScaleData &data,
                               const ScaleFitOptions &options,
                               const std::vector<std::size_t> &given) {
  std::vector<std::size_t> use = given;
  if (use.empty())
    for (std::size_t i = 0; i < data.size(); ++i)
      if (!data.outlier[i])
        use.push_back(i);
  ScaleFitResult result;
  result.observations = use.size();
  const std::size_t n = model.size();
  if (use.empty() || n == 0)
    return result;

  std::vector<double> g, merged;
  double phi = target(model, data, use, options, &g, &merged);
  result.target_start = phi;
  double mu = 1e-3;
  // The observations of each group being fitted.
  std::vector<std::vector<std::size_t>> members(data.unique.size());
  for (std::size_t i : use)
    members[data.group[i]].push_back(i);
  for (int it = 0; it < options.max_iterations; ++it) {
    result.iterations = it + 1;
    // Normal equations of the residuals sqrt(w)(I - g <I>), differentiated
    // with <I> as a function of the parameters (variable projection, in
    // Kaufman's form): d<I>/dp = -<I> G_h, G_h = sum w g dg/dp / sum w g^2.
    // Holding <I> instead gives the right gradient but overstates the
    // curvature in every direction <I> can partly follow, so the steps come
    // out short: three fits of 50 steps on a 300 image sweep ended still
    // creeping. The row for observation i is dg_i/dp - g_i G_h.
    std::vector<double> N, b;
    normal_equations(model, data, members, g, merged, options, N, b);
    bool improved = false;
    for (int tries = 0; tries < 12 && !improved; ++tries) {
      std::vector<double> A = N, step = b;
      for (std::size_t a = 0; a < n; ++a)
        A[a * n + a] += mu * std::fmax(N[a * n + a], 1e-12);
      if (!solve_spd(A.data(), step.data(), n)) {
        mu *= 10.0;
        continue;
      }
      ScaleModel trial = model;
      for (std::size_t a = 0; a < n; ++a)
        trial.parameters[a] += step[a];
      trial.normalise();
      std::vector<double> g2, m2;
      const double phi2 = target(trial, data, use, options, &g2, &m2);
      if (phi2 < phi) {
        const double gain = phi - phi2;
        model = trial;
        g = std::move(g2);
        merged = std::move(m2);
        phi = phi2;
        mu = std::fmax(mu / 10.0, 1e-9);
        improved = true;
        if (gain < 1e-8 * std::fmax(1.0, phi)) {
          result.converged = true;
          result.target_end = phi;
          return result;
        }
      } else {
        mu *= 10.0;
      }
    }
    if (!improved) {
      result.converged = true; // no step lowers the target: a minimum
      break;
    }
  }
  result.target_end = phi;
  return result;
}

namespace {

//: Weighted sums over a group's observations that are not outliers.
struct GroupSums {
  std::vector<double> num, den;
  std::vector<std::size_t> count;
};

GroupSums group_sums(const ScaleData &data, const std::vector<double> &g) {
  GroupSums s;
  s.num.assign(data.unique.size(), 0.0);
  s.den.assign(data.unique.size(), 0.0);
  s.count.assign(data.unique.size(), 0);
  for (std::size_t i = 0; i < data.size(); ++i) {
    if (data.outlier[i])
      continue;
    const double w = 1.0 / data.variance[i];
    s.num[data.group[i]] += w * g[i] * data.intensity[i];
    s.den[data.group[i]] += w * g[i] * g[i];
    ++s.count[data.group[i]];
  }
  return s;
}

//: The standard normal quantile, by Acklam's rational approximation refined
//: with a Newton step on erfc: better than 1e-12.
double normal_quantile(double p) {
  const double a[6] = {-3.969683028665376e+01, 2.209460984245205e+02,
                       -2.759285104469687e+02, 1.383577518672690e+02,
                       -3.066479806614716e+01, 2.506628277459239e+00};
  const double b[5] = {-5.447609879822406e+01, 1.615858368580409e+02,
                       -1.556989798598866e+02, 6.680131188771972e+01,
                       -1.328068155288572e+01};
  const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                       -2.400758277161838e+00, -2.549732539343734e+00,
                       4.374664141464968e+00,  2.938163982698783e+00};
  const double d[4] = {7.784695709041462e-03, 3.224671290700398e-01,
                       2.445134137142996e+00, 3.754408661907416e+00};
  double x;
  if (p < 0.02425) {
    const double q = std::sqrt(-2.0 * std::log(p));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  } else if (p > 1.0 - 0.02425) {
    const double q = std::sqrt(-2.0 * std::log(1.0 - p));
    x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  } else {
    const double q = p - 0.5, r = q * q;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) *
        q /
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
  }
  const double e = 0.5 * std::erfc(-x / std::sqrt(2.0)) - p;
  const double u = e * std::sqrt(2.0 * std::acos(-1.0)) * std::exp(x * x / 2.0);
  return x - u / (1.0 + x * u / 2.0);
}

} // namespace

std::size_t reject_outliers(ScaleData &data, const std::vector<double> &g,
                            double zmax) {
  std::fill(data.outlier.begin(), data.outlier.end(), false);
  std::vector<std::vector<std::size_t>> members(data.unique.size());
  for (std::size_t i = 0; i < data.size(); ++i)
    members[data.group[i]].push_back(i);
  std::size_t flagged = 0;
  for (const std::vector<std::size_t> &m : members) {
    while (true) {
      double sn = 0.0, sd = 0.0;
      std::size_t live = 0;
      for (std::size_t i : m) {
        if (data.outlier[i])
          continue;
        const double w = 1.0 / data.variance[i];
        sn += w * g[i] * data.intensity[i];
        sd += w * g[i] * g[i];
        ++live;
      }
      if (live < 3)
        break; // with two, which of them is wrong cannot be said
      double worst = 0.0;
      std::size_t at = m.front();
      for (std::size_t i : m) {
        if (data.outlier[i])
          continue;
        const double w = 1.0 / data.variance[i];
        const double num = sn - w * g[i] * data.intensity[i];
        const double den = sd - w * g[i] * g[i];
        const double z = std::abs(data.intensity[i] - g[i] * num / den) /
                         std::sqrt(data.variance[i] + g[i] * g[i] / den);
        if (z > worst) {
          worst = z;
          at = i;
        }
      }
      if (worst <= zmax)
        break;
      data.outlier[at] = true;
      ++flagged;
    }
  }
  return flagged;
}

ErrorModel refine_error_model(const ScaleData &data,
                              const std::vector<double> &g) {
  const GroupSums s = group_sums(data, g);
  // Groups fit to judge an error model by, as dials.scale: <I> above 25 and
  // <I / sigma^2> above 0.85, with at least two observations.
  std::vector<double> ratio(data.unique.size(), 0.0);
  for (std::size_t i = 0; i < data.size(); ++i)
    if (!data.outlier[i])
      ratio[data.group[i]] +=
          data.intensity[i] / (data.variance_before[i] + data.scale_term_at(i));
  std::vector<std::size_t> use;
  for (std::size_t i = 0; i < data.size(); ++i) {
    const std::size_t h = data.group[i];
    if (data.outlier[i] || s.count[h] < 2 || !(s.den[h] > 0.0))
      continue;
    const double mean = s.num[h] / s.den[h];
    if (mean > 25.0 && ratio[h] / static_cast<double>(s.count[h]) > 0.85)
      use.push_back(i);
  }
  ErrorModel em;
  em.used = use.size();
  if (use.size() < 50)
    return em;
  // Normalised deviations of x_i = I_i / g_i from the merged <I_h>, with the
  // EXACT variance of the difference rather than eqn 12's prefactor. With
  // <I_h> = sum c_j x_j, weights held as dials.scale holds them, and v_j the
  // variance of x_j under the candidate error model,
  //     Var(x_i - <I_h>) = v_i (1 - 2 c_i) + sum_j c_j^2 v_j,
  // which is v (n - 1) / n when every weight is equal. dials.scale multiplies
  // by sqrt((n - 1) / n) instead (calc_deltahl), which leaves the deviations'
  // spread at (n - 1) / n and a low by that factor: planted a = 1.3 in groups
  // of eight came back 1.158.
  std::vector<double> c(data.size(), 0.0);
  for (std::size_t i : use) {
    const std::size_t h = data.group[i];
    c[i] = (g[i] * g[i] / data.variance[i]) / s.den[h];
  }
  std::vector<double> spread(data.unique.size(), 0.0), v(data.size(), 0.0);
  const auto deviations = [&](double a, double b, std::vector<double> *out) {
    std::fill(spread.begin(), spread.end(), 0.0);
    for (std::size_t i : use) {
      v[i] = a * a *
             (data.variance_before[i] + data.scale_term_at(i) +
              b * b * data.intensity[i] * data.intensity[i]) /
             (g[i] * g[i]);
      spread[data.group[i]] += c[i] * c[i] * v[i];
    }
    out->clear();
    for (std::size_t i : use) {
      const std::size_t h = data.group[i];
      const double var =
          std::fmax(v[i] * (1.0 - 2.0 * c[i]) + spread[h], 1e-12 * v[i]);
      out->push_back((data.intensity[i] / g[i] - s.num[h] / s.den[h]) /
                     std::sqrt(var));
    }
  };
  // Intensity bins for b: logarithmically spaced over the used intensities.
  double lo = HUGE_VAL, hi = 0.0;
  for (std::size_t i : use) {
    const double v = std::fmax(data.intensity[i], 1e-3);
    lo = std::fmin(lo, v);
    hi = std::fmax(hi, v);
  }
  const int nbins = 20;
  std::vector<int> bin(use.size());
  for (std::size_t k = 0; k < use.size(); ++k) {
    const double v = std::fmax(data.intensity[use[k]], 1e-3);
    bin[k] = hi > lo ? std::min(nbins - 1,
                                static_cast<int>(nbins * std::log(v / lo) /
                                                 std::log(hi / lo)))
                     : 0;
  }
  std::vector<double> dev;
  for (int round = 0; round < 20; ++round) {
    const double a0 = em.a, b0 = em.b;
    // a: the slope of the central normal probability plot, |x| < 1.5.
    deviations(em.a, em.b, &dev);
    std::vector<double> sorted = dev;
    std::sort(sorted.begin(), sorted.end());
    const double n = static_cast<double>(sorted.size());
    double sxx = 0.0, sxy = 0.0;
    for (std::size_t k = 0; k < sorted.size(); ++k) {
      const double x = normal_quantile((static_cast<double>(k) + 0.5) / n);
      if (std::abs(x) < 1.5) {
        sxx += x * x;
        sxy += x * sorted[k];
      }
    }
    if (sxx > 0.0)
      em.a *= sxy / sxx;
    // b: minimise eqn 17, sum w_i [(0.5 - v_i)^2 + 1/v_i], the bin variances
    // v_i of the deviations and w_i proportional to the bin's mean intensity.
    const auto phi = [&](double b) {
      deviations(em.a, b, &dev);
      std::vector<double> sum(nbins, 0.0), sum2(nbins, 0.0), mean_i(nbins, 0.0);
      std::vector<std::size_t> count(nbins, 0);
      for (std::size_t k = 0; k < dev.size(); ++k) {
        sum[bin[k]] += dev[k];
        sum2[bin[k]] += dev[k] * dev[k];
        mean_i[bin[k]] += data.intensity[use[k]];
        ++count[bin[k]];
      }
      double f = 0.0;
      for (int j = 0; j < nbins; ++j) {
        if (count[j] < 10)
          continue;
        const double c = static_cast<double>(count[j]);
        const double v =
            std::fmax(sum2[j] / c - (sum[j] / c) * (sum[j] / c), 1e-6);
        f += (mean_i[j] / c) * ((0.5 - v) * (0.5 - v) + 1.0 / v);
      }
      return f;
    };
    double x0 = 0.0, x1 = 0.5; // golden section over b
    const double r = 0.5 * (std::sqrt(5.0) - 1.0);
    double c1 = x1 - r * (x1 - x0), c2 = x0 + r * (x1 - x0);
    double f1 = phi(c1), f2 = phi(c2);
    for (int it = 0; it < 60; ++it) {
      if (f1 < f2) {
        x1 = c2;
        c2 = c1;
        f2 = f1;
        c1 = x1 - r * (x1 - x0);
        f1 = phi(c1);
      } else {
        x0 = c1;
        c1 = c2;
        f1 = f2;
        c2 = x0 + r * (x1 - x0);
        f2 = phi(c2);
      }
    }
    em.b = 0.5 * (x0 + x1);
    if (std::abs(em.a - a0) < 1e-5 && std::abs(em.b - b0) < 1e-6)
      break;
  }
  return em;
}

void apply_error_model(ScaleData &data, const ErrorModel &model) {
  for (std::size_t i = 0; i < data.size(); ++i)
    data.variance[i] =
        model.a * model.a *
        (data.variance_before[i] + data.scale_term_at(i) +
         model.b * model.b * data.intensity[i] * data.intensity[i]);
}

void combine_intensities(ScaleData &data, double i_mid) {
  for (std::size_t i = 0; i < data.size(); ++i) {
    double w = 1.0;
    if (!data.has_sum[i])
      w = 1.0;
    else if (std::isinf(i_mid))
      w = 0.0;
    else if (i_mid > 0.0)
      w = 1.0 / (1.0 + std::pow(std::fmax(data.sum[i], 0.0) / i_mid, 3.0));
    data.intensity[i] = w * data.prf[i] + (1.0 - w) * data.sum[i];
    const double sp = std::sqrt(data.prf_variance[i]),
                 ss = std::sqrt(data.sum_variance[i]);
    const double sigma = w * sp + (1.0 - w) * (data.has_sum[i] ? ss : sp);
    data.variance[i] = sigma * sigma;
    data.variance_before[i] = sigma * sigma;
  }
}

double rmeas(const ScaleData &data, const std::vector<double> &g) {
  const GroupSums s = group_sums(data, g);
  double top = 0.0, bottom = 0.0;
  for (std::size_t i = 0; i < data.size(); ++i) {
    const std::size_t h = data.group[i];
    if (data.outlier[i] || s.count[h] < 2 || !(s.den[h] > 0.0))
      continue;
    const double n = static_cast<double>(s.count[h]);
    const double mean = s.num[h] / s.den[h];
    top += std::sqrt(n / (n - 1.0)) * std::abs(data.intensity[i] / g[i] - mean);
    bottom += data.intensity[i] / g[i];
  }
  return bottom > 0.0 ? top / bottom : 0.0;
}

double choose_intensity_combination(ScaleData &data,
                                    const std::vector<double> &g) {
  double biggest = 1.0;
  for (std::size_t i = 0; i < data.size(); ++i)
    if (data.has_sum[i])
      biggest = std::fmax(biggest, data.sum[i]);
  std::vector<double> candidates = {0.0, HUGE_VAL};
  for (double m = 10.0; m < biggest; m *= 10.0)
    candidates.push_back(m);
  double best = 0.0, best_r = HUGE_VAL;
  for (double m : candidates) {
    combine_intensities(data, m);
    const double r = rmeas(data, g);
    if (r < best_r) {
      best_r = r;
      best = m;
    }
  }
  combine_intensities(data, best);
  return best;
}

std::vector<MergingShell> merging_statistics(const ScaleData &data,
                                             const std::vector<double> &g,
                                             const SpaceGroup &group,
                                             const Crystal &crystal, int shells,
                                             MergingShell *overall) {
  shells = std::max(shells, 1);
  // The range from the same d as the shells and the possible reflections: the
  // crystal's, at each group's unique index.
  std::vector<bool> seen(data.unique.size(), false);
  for (std::size_t i = 0; i < data.size(); ++i)
    if (!data.outlier[i])
      seen[data.group[i]] = true;
  double lo = HUGE_VAL, hi = 0.0; // 1/d^3
  for (std::size_t h = 0; h < data.unique.size(); ++h) {
    if (!seen[h])
      continue;
    const Miller &u = data.unique[h];
    const double d = crystal.d_spacing(u[0], u[1], u[2]);
    const double v = 1.0 / (d * d * d);
    lo = std::fmin(lo, v);
    hi = std::fmax(hi, v);
  }
  std::vector<MergingShell> out(static_cast<std::size_t>(shells));
  if (!(hi > 0.0))
    return out;
  const auto shell_of = [&](double d) {
    const double v = 1.0 / (d * d * d);
    const auto s =
        hi > lo ? static_cast<int>((v - lo) / (hi - lo) * shells) : 0;
    return static_cast<std::size_t>(std::min(std::max(s, 0), shells - 1));
  };
  for (int s = 0; s < shells; ++s) {
    out[static_cast<std::size_t>(s)].d_max =
        std::cbrt(1.0 / (lo + (hi - lo) * s / shells));
    out[static_cast<std::size_t>(s)].d_min =
        std::cbrt(1.0 / (lo + (hi - lo) * (s + 1) / shells));
  }
  const double d_min = out.back().d_min, d_max = out.front().d_max;

  // Per group: the scaled observations x = I/g, their variances, and a random
  // half for CC half.
  std::vector<std::vector<std::size_t>> members(data.unique.size());
  for (std::size_t i = 0; i < data.size(); ++i)
    if (!data.outlier[i])
      members[data.group[i]].push_back(i);
  std::mt19937 rng(20);
  std::vector<std::vector<std::pair<double, double>>> halves(out.size() + 1);
  const auto accumulate = [&](MergingShell &m, std::size_t h,
                              std::size_t slot) {
    const std::vector<std::size_t> &o = members[h];
    if (o.empty())
      return;
    m.observations += o.size();
    ++m.unique;
    double sw = 0.0, swx = 0.0, sx = 0.0;
    for (std::size_t i : o) {
      const double x = data.intensity[i] / g[i];
      const double w = g[i] * g[i] / data.variance[i];
      sw += w;
      swx += w * x;
      sx += x;
    }
    m.mean_i += swx / sw;
    m.i_over_sigma +=
        (swx / sw) * std::sqrt(sw); // <I> / sigma(<I>), summed for now
    if (o.size() < 2)
      return;
    const double n = static_cast<double>(o.size()), mean = sx / n;
    double dev = 0.0;
    for (std::size_t i : o)
      dev += std::abs(data.intensity[i] / g[i] - mean);
    m.rmerge += dev;
    m.rmeas += std::sqrt(n / (n - 1.0)) * dev;
    m.rpim += std::sqrt(1.0 / (n - 1.0)) * dev;
    m.cc_half += sx; // the denominator, summed for now
    std::vector<std::size_t> shuffled = o;
    std::shuffle(shuffled.begin(), shuffled.end(), rng);
    double a = 0.0, b = 0.0;
    const std::size_t split = shuffled.size() / 2;
    for (std::size_t k = 0; k < shuffled.size(); ++k)
      (k < split ? a : b) += data.intensity[shuffled[k]] / g[shuffled[k]];
    halves[slot].emplace_back(a / static_cast<double>(split),
                              b / static_cast<double>(shuffled.size() - split));
  };
  MergingShell all;
  for (std::size_t h = 0; h < members.size(); ++h) {
    if (members[h].empty())
      continue;
    // The group's d from the crystal and its unique index, the same function
    // the possible reflections are counted with: the d column is from the
    // scan-varying crystal at each observation, and near a shell's edge the two
    // disagreed enough to put one shell at 100.2 per cent.
    const Miller &u = data.unique[h];
    const std::size_t s = shell_of(crystal.d_spacing(u[0], u[1], u[2]));
    accumulate(out[s], h, s);
    accumulate(all, h, out.size());
  }
  // Possible reflections: every index the cell allows to d_min, in the
  // asymmetric unit and not absent.
  std::vector<std::size_t> possible(out.size(), 0);
  const UnitCell cell = crystal.cell();
  const int hmax = static_cast<int>(std::ceil(cell.a / d_min)),
            kmax = static_cast<int>(std::ceil(cell.b / d_min)),
            lmax = static_cast<int>(std::ceil(cell.c / d_min));
  for (int hh = -hmax; hh <= hmax; ++hh)
    for (int kk = -kmax; kk <= kmax; ++kk)
      for (int ll = -lmax; ll <= lmax; ++ll) {
        if (hh == 0 && kk == 0 && ll == 0)
          continue;
        const double d = crystal.d_spacing(hh, kk, ll);
        if (d < d_min || d > d_max)
          continue;
        const Miller m{hh, kk, ll};
        if (group.absent(m) || group.unique(m) != m)
          continue;
        ++possible[shell_of(d)];
        ++all.possible;
      }
  const auto finish = [&](MergingShell &m,
                          const std::vector<std::pair<double, double>> &pairs) {
    const double denominator = m.cc_half;
    m.multiplicity =
        m.unique ? static_cast<double>(m.observations) / m.unique : 0.0;
    m.completeness =
        m.possible ? static_cast<double>(m.unique) / m.possible : 0.0;
    m.i_over_sigma = m.unique ? m.i_over_sigma / m.unique : 0.0;
    m.mean_i = m.unique ? m.mean_i / m.unique : 0.0;
    m.rmerge = denominator > 0.0 ? m.rmerge / denominator : 0.0;
    m.rmeas = denominator > 0.0 ? m.rmeas / denominator : 0.0;
    m.rpim = denominator > 0.0 ? m.rpim / denominator : 0.0;
    // Pearson over the half-dataset means.
    const double n = static_cast<double>(pairs.size());
    double ma = 0.0, mb = 0.0;
    for (const auto &[a, b] : pairs) {
      ma += a / n;
      mb += b / n;
    }
    double sab = 0.0, saa = 0.0, sbb = 0.0;
    for (const auto &[a, b] : pairs) {
      sab += (a - ma) * (b - mb);
      saa += (a - ma) * (a - ma);
      sbb += (b - mb) * (b - mb);
    }
    m.cc_half = pairs.size() > 2 && saa > 0.0 && sbb > 0.0
                    ? sab / std::sqrt(saa * sbb)
                    : 0.0;
  };
  for (std::size_t s = 0; s < out.size(); ++s) {
    out[s].possible = possible[s];
    finish(out[s], halves[s]);
  }
  all.d_max = d_max;
  all.d_min = d_min;
  finish(all, halves[out.size()]);
  if (overall)
    *overall = all;
  return out;
}

ScaleRun scale_sweep(const ExperimentList &experiments,
                     const Table &reflections, const SpaceGroup &group,
                     const ScaleRunOptions &options) {
  ScaleRun run;
  const Experiment &e = experiments[0];
  const double degrees = std::abs(e.scan.osc_width) * e.scan.num_images();
  ScaleModelShape shape = default_shape(degrees);
  if (!options.absorption)
    shape.lmax = 0;
  ScaleDataOptions data_options;
  data_options.d_min = options.d_min;
  run.data =
      build_scale_data(experiments, reflections, group, shape, data_options);
  ScaleData &data = run.data;
  run.model = ScaleModel(shape);
  if (data.size() == 0)
    return run;

  const auto fit = [&]() {
    const std::vector<std::size_t> use = select_for_fitting(data);
    run.fitted_on = use.size();
    run.fits.push_back(fit_scale_model(run.model, data, options.fit, use));
    run.g = inverse_scales(run.model, data);
  };
  run.g.assign(data.size(), 1.0);
  reject_outliers(data, run.g); // on the unscaled intensities
  fit();
  reject_outliers(data, run.g);
  if (options.combine) {
    run.i_mid = choose_intensity_combination(data, run.g);
    reject_outliers(data, run.g);
  }
  run.error_model = refine_error_model(data, run.g);
  apply_error_model(data, run.error_model);
  fit();
  reject_outliers(data, run.g);
  run.error_model = refine_error_model(data, run.g);
  apply_error_model(data, run.error_model);
  fit();
  // The scale's uncertainty into each observation's variance BEFORE the final
  // error model, so that a and b correct what remains after it rather than
  // absorbing it. (dials.scale applies its error model first and then
  // multiplies each variance by 1 + sigma_g / g, a factor linear in the
  // fractional error and the same at every intensity; propagated, the term is
  // I^2 var(g) / g^2.)
  run.covariance = parameter_covariance(run.model, data, options.fit,
                                        select_for_fitting(data));
  run.g_variance = inverse_scale_variances(run.model, data, run.covariance);
  propagate_scale_variances(data, run.g, run.g_variance);
  run.outliers = reject_outliers(data, run.g);
  run.error_model = refine_error_model(data, run.g);
  apply_error_model(data, run.error_model);
  return run;
}

void write_scaling(Table &reflections, const ScaleData &data,
                   const std::vector<double> &g,
                   const std::vector<double> &g_variance) {
  const std::size_t n = reflections.nrows;
  Column &isf = reflections.real_column("inverse_scale_factor", "double", 1);
  Column &isf_var =
      reflections.real_column("inverse_scale_factor_variance", "double", 1);
  Column &value = reflections.real_column("intensity.scale.value", "double", 1);
  Column &variance =
      reflections.real_column("intensity.scale.variance", "double", 1);
  // A copy of the flags, changed and set back. int_column() makes a new column
  // of zeroes: used here, it wiped every flag integration had set.
  Column flags = reflections.at("flags");
  std::vector<bool> scaled(n, false);
  for (std::size_t r = 0; r < n; ++r) {
    isf.reals[r] = 1.0;
    isf_var.reals[r] = 0.0;
    value.reals[r] = 0.0;
    variance.reals[r] = 0.0;
    flags.ints[r] &=
        ~(flag::kOutlierInScaling | flag::kExcludedForScaling | flag::kScaled);
  }
  for (std::size_t i = 0; i < data.size(); ++i) {
    const std::size_t r = data.row[i];
    scaled[r] = true;
    isf.reals[r] = g[i];
    isf_var.reals[r] = i < g_variance.size() ? g_variance[i] : 0.0;
    value.reals[r] = data.intensity[i];
    variance.reals[r] = data.variance[i];
    flags.ints[r] |= data.outlier[i] ? flag::kOutlierInScaling : flag::kScaled;
  }
  for (std::size_t r = 0; r < n; ++r)
    if (!scaled[r])
      flags.ints[r] |= flag::kExcludedForScaling;
  reflections.set("flags", std::move(flags));
}

ParameterCovariance
parameter_covariance(const ScaleModel &model, const ScaleData &data,
                     const ScaleFitOptions &options,
                     const std::vector<std::size_t> &given) {
  ParameterCovariance out;
  std::vector<std::size_t> use = given;
  if (use.empty())
    for (std::size_t i = 0; i < data.size(); ++i)
      if (!data.outlier[i])
        use.push_back(i);
  const std::size_t n = model.size();
  if (use.empty() || n == 0)
    return out;
  std::vector<std::vector<std::size_t>> members(data.unique.size());
  for (std::size_t i : use)
    members[data.group[i]].push_back(i);
  std::vector<double> g, merged;
  const double phi = target(model, data, use, options, &g, &merged);
  std::vector<double> N, b;
  normal_equations(model, data, members, g, merged, options, N, b);

  // Z: the directions keeping the scale's sum and the relative B's sum fixed,
  // e_i - e_last within each block; the absorption terms are free.
  const ScaleModelShape &s = model.shape();
  std::vector<std::vector<double>> Z;
  const auto block = [&](std::size_t first, std::size_t count) {
    for (std::size_t i = 0; i + 1 < count; ++i) {
      std::vector<double> z(n, 0.0);
      z[first + i] = 1.0;
      z[first + count - 1] = -1.0;
      Z.push_back(std::move(z));
    }
  };
  block(0, s.scale_points);
  block(model.first_decay(), s.decay_points);
  for (std::size_t k = 0; k < harmonic_count(s.lmax); ++k) {
    std::vector<double> z(n, 0.0);
    z[model.first_absorption() + k] = 1.0;
    Z.push_back(std::move(z));
  }
  const std::size_t m = Z.size();
  if (m == 0)
    return out;
  std::vector<double> NZ(n * m, 0.0), M(m * m, 0.0);
  for (std::size_t a = 0; a < n; ++a)
    for (std::size_t j = 0; j < m; ++j)
      for (std::size_t c = 0; c < n; ++c)
        NZ[a * m + j] += N[a * n + c] * Z[j][c];
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < m; ++j)
      for (std::size_t a = 0; a < n; ++a)
        M[i * m + j] += Z[i][a] * NZ[a * m + j];
  std::vector<double> Minv(m * m, 0.0);
  for (std::size_t j = 0; j < m; ++j) {
    std::vector<double> A = M, e(m, 0.0);
    e[j] = 1.0;
    if (!solve_spd(A.data(), e.data(), m))
      return out;
    for (std::size_t i = 0; i < m; ++i)
      Minv[i * m + j] = e[i];
  }
  // The degrees of freedom: observations, less the merged intensities, less
  // the free parameters.
  std::size_t groups = 0;
  for (const auto &group_members : members)
    if (!group_members.empty())
      ++groups;
  const std::size_t fitted = groups + m;
  out.degrees_of_freedom = use.size() > fitted ? use.size() - fitted : 1;
  double restraint = 0.0;
  for (std::size_t k = 0; k < s.decay_points; ++k)
    restraint += options.decay_restraint *
                 model.parameters[model.first_decay() + k] *
                 model.parameters[model.first_decay() + k];
  for (std::size_t k = 0; k < harmonic_count(s.lmax); ++k)
    restraint += options.absorption_restraint *
                 model.parameters[model.first_absorption() + k] *
                 model.parameters[model.first_absorption() + k];
  out.goodness_of_fit =
      (phi - restraint) / static_cast<double>(out.degrees_of_freedom);
  out.matrix.assign(n * n, 0.0);
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < m; ++j) {
      const double v = out.goodness_of_fit * Minv[i * m + j];
      if (v == 0.0)
        continue;
      for (std::size_t a = 0; a < n; ++a) {
        if (Z[i][a] == 0.0)
          continue;
        for (std::size_t c = 0; c < n; ++c)
          out.matrix[a * n + c] += Z[i][a] * v * Z[j][c];
      }
    }
  out.ok = true;
  return out;
}

std::vector<double>
inverse_scale_variances(const ScaleModel &model, const ScaleData &data,
                        const ParameterCovariance &covariance) {
  std::vector<double> out(data.size(), 0.0);
  if (!covariance.ok)
    return out;
  const std::size_t n = model.size();
  std::vector<std::pair<std::size_t, double>> grad;
  for (std::size_t i = 0; i < data.size(); ++i) {
    model.inverse_scale(data.observation[i], &grad);
    double v = 0.0;
    for (const auto &[a, da] : grad)
      for (const auto &[c, dc] : grad)
        v += da * covariance.matrix[a * n + c] * dc;
    out[i] = std::fmax(v, 0.0);
  }
  return out;
}

void propagate_scale_variances(ScaleData &data, const std::vector<double> &g,
                               const std::vector<double> &g_variance) {
  data.scale_term.resize(data.size(), 0.0);
  for (std::size_t i = 0; i < data.size(); ++i)
    data.scale_term[i] = g[i] > 0.0 ? data.intensity[i] * data.intensity[i] *
                                          g_variance[i] / (g[i] * g[i])
                                    : 0.0;
}

} // namespace mxi
