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
  std::vector<std::pair<std::size_t, double>> grad;
  for (int it = 0; it < options.max_iterations; ++it) {
    result.iterations = it + 1;
    // Normal equations of the residuals sqrt(w)(I - g <I>), <I> held.
    std::vector<double> N(n * n, 0.0), b(n, 0.0);
    for (std::size_t i : use) {
      const double m = merged[data.group[i]];
      const double w = 1.0 / data.variance[i];
      model.inverse_scale(data.observation[i], &grad);
      const double r = data.intensity[i] - g[i] * m;
      // J = -m dg/dp; accumulate J^T J and -J^T r (the descent direction).
      for (const auto &[a, da] : grad) {
        b[a] += w * m * da * r;
        for (const auto &[c, dc] : grad)
          N[a * n + c] += w * m * m * da * dc;
      }
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
        if (gain < 1e-10 * std::fmax(1.0, phi)) {
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

} // namespace mxi
