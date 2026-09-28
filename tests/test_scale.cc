#include <cmath>
#include <random>

#include "../src/scale.hh"
#include "check.hh"

namespace mxi {

TEST(the_merged_intensity_is_the_weighted_least_squares_one) {
  // Two observations of one reflection: g = (1, 2), I = (10, 22), variances
  // (1, 4). <I> = (1*1*10 + 0.25*2*22) / (1*1 + 0.25*4) = 21 / 2 = 10.5.
  ScaleData data;
  data.intensity = {10.0, 22.0};
  data.variance = {1.0, 4.0};
  data.group = {0, 0};
  data.unique = {Miller{1, 2, 3}};
  data.outlier = {false, false};
  data.observation.resize(2);
  const std::vector<double> m = merged_intensities(data, {1.0, 2.0});
  check::close(m[0], 10.5, 1e-15, "10.5");
  data.outlier[1] = true;
  check::close(merged_intensities(data, {1.0, 2.0})[0], 10.0, 1e-15,
               "and an outlier takes no part");
}

namespace {

//: Observations made from a known model: 1500 reflections seen 6 to 10 times
//: each, across the scan and resolution, with absorption at random directions,
//: and noise of one per cent.
ScaleData planted(const ScaleModel &truth, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  std::exponential_distribution<double> strength(1.0 / 1000.0);
  std::normal_distribution<double> noise(0.0, 1.0);
  ScaleData data;
  std::vector<double> y;
  for (std::size_t h = 0; h < 1500; ++h) {
    data.unique.push_back(Miller{static_cast<int>(h), 0, 0});
    const double inv_d2 =
        1.0 / 64.0 + uni(rng) * (1.0 / 2.25 - 1.0 / 64.0); // 8 to 1.5 A
    const double truth_i = 50.0 + strength(rng);
    const int seen = 6 + static_cast<int>(uni(rng) * 5.0);
    for (int k = 0; k < seen; ++k) {
      ScaleObservation o;
      o.rotation = uni(rng);
      o.time = o.rotation;
      o.inv_2d2 = 0.5 * inv_d2;
      const Vec3 a{noise(rng), noise(rng), noise(rng)},
          b{noise(rng), noise(rng), noise(rng)};
      std::vector<double> ya, yb;
      real_spherical_harmonics(truth.shape().lmax, a.normalized(), &ya);
      real_spherical_harmonics(truth.shape().lmax, b.normalized(), &yb);
      for (std::size_t j = 0; j < ya.size(); ++j)
        o.absorption.push_back(0.5 * (ya[j] + yb[j]));
      const double g = truth.inverse_scale(o);
      const double sigma = 0.01 * g * truth_i;
      data.intensity.push_back(g * truth_i + sigma * noise(rng));
      data.variance.push_back(sigma * sigma);
      data.observation.push_back(std::move(o));
      data.group.push_back(h);
      data.outlier.push_back(false);
    }
  }
  return data;
}

} // namespace

TEST(a_planted_scale_decay_and_absorption_are_recovered) {
  ScaleModel truth({6, 5, 2});
  const double c[6] = {0.9, 1.05, 1.15, 1.0, 0.85, 1.05};
  const double b[5] = {1.0, 0.5, 0.0, -0.5,
                       -1.0}; // mean zero: its offset is unseen
  for (std::size_t i = 0; i < 6; ++i)
    truth.parameters[i] = c[i];
  for (std::size_t i = 0; i < 5; ++i)
    truth.parameters[truth.first_decay() + i] = b[i];
  for (std::size_t k = 0; k < 8; ++k)
    truth.parameters[truth.first_absorption() + k] =
        0.02 * std::cos(1.3 * k + 0.4);
  truth.normalise();
  const ScaleData data = planted(truth, 7);

  ScaleModel fit({6, 5, 2});
  const ScaleFitResult r = fit_scale_model(fit, data);
  check::is_true(r.converged, "converged in " + std::to_string(r.iterations));
  check::is_true(r.target_end < r.target_start, "and lowered the target");
  // The scale's shape, which is all that is seen: both normalised to mean one.
  double worst_c = 0.0, worst_b = 0.0, worst_p = 0.0;
  for (std::size_t i = 0; i < 6; ++i)
    worst_c =
        std::fmax(worst_c, std::abs(fit.parameters[i] - truth.parameters[i]));
  // B up to a common offset, which the merged intensities absorb.
  double off = 0.0;
  for (std::size_t i = 0; i < 5; ++i)
    off += (fit.parameters[fit.first_decay() + i] -
            truth.parameters[truth.first_decay() + i]) /
           5.0;
  for (std::size_t i = 0; i < 5; ++i)
    worst_b = std::fmax(worst_b,
                        std::abs(fit.parameters[fit.first_decay() + i] - off -
                                 truth.parameters[truth.first_decay() + i]));
  for (std::size_t k = 0; k < 8; ++k)
    worst_p = std::fmax(
        worst_p, std::abs(fit.parameters[fit.first_absorption() + k] -
                          truth.parameters[truth.first_absorption() + k]));
  check::is_true(worst_c < 0.003, "the scale, to " + std::to_string(worst_c));
  check::is_true(worst_b < 0.05,
                 "the relative B, to " + std::to_string(worst_b) + " A^2");
  check::is_true(worst_p < 0.003,
                 "the absorption surface, to " + std::to_string(worst_p));
}

TEST(fitting_selects_whole_groups_until_both_thresholds_are_met) {
  ScaleModel truth({6, 0, 0});
  const ScaleData data = planted(truth, 3);
  const std::vector<std::size_t> s = select_for_fitting(data, 200, 3000);
  std::vector<std::size_t> count(data.unique.size(), 0),
      total(data.unique.size(), 0);
  for (std::size_t i : s)
    ++count[data.group[i]];
  for (std::size_t i = 0; i < data.size(); ++i)
    ++total[data.group[i]];
  std::size_t groups = 0;
  bool whole = true;
  for (std::size_t h = 0; h < count.size(); ++h) {
    if (count[h] == 0)
      continue;
    ++groups;
    whole = whole && count[h] == total[h];
  }
  check::is_true(whole, "every group taken is taken whole");
  check::is_true(groups >= 200 && s.size() >= 3000, "both thresholds met");
  check::is_true(s.size() < data.size(),
                 "and not everything, when that is enough");
}

} // namespace mxi
