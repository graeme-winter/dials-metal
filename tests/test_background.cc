// The robust GLM background, Parkhurst et al. (2016).

#include <cmath>
#include <cstdint>
#include <vector>

#include "../src/background.h"
#include "check.h"

namespace mxi {

namespace {

//: Poisson draws, deterministic, so a test cannot pass or fail by luck.
std::vector<double> poisson_sample(double mu, std::size_t n, std::uint64_t seed) {
  std::vector<double> out;
  out.reserve(n);
  std::uint64_t state = seed;
  const auto uniform = [&state]() {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>((state >> 11) & ((1ULL << 53) - 1)) /
           static_cast<double>(1ULL << 53);
  };
  for (std::size_t i = 0; i < n; ++i) {
    // Knuth: multiply uniforms until the product falls below exp(-mu).
    const double limit = std::exp(-mu);
    double product = 1.0;
    int k = -1;
    do {
      product *= std::fmax(uniform(), 1e-300);
      ++k;
    } while (product > limit);
    out.push_back(static_cast<double>(k));
  }
  return out;
}

double mean_of(const std::vector<double> &v) {
  double s = 0.0;
  for (double x : v) s += x;
  return v.empty() ? 0.0 : s / static_cast<double>(v.size());
}

}  // namespace

TEST(with_no_clipping_the_glm_is_the_plain_mean) {
  // The Huber weight with a large tuning constant does nothing, and the
  // estimator must then be the maximum likelihood one, which for a Poisson is
  // the sample mean. If this fails the correction terms are wrong, because
  // they are what has to vanish here: C1 goes to zero and C2 to sqrt(mu).
  const std::vector<double> v = poisson_sample(0.4, 4000, 11);
  BackgroundOptions options;
  options.tuning = 1e6;
  const BackgroundResult r = glm_background(v, options);
  check::is_true(r.valid && r.converged, "converged");
  check::close(r.mean, mean_of(v), 1e-9 * mean_of(v), "and it is the mean");
}

TEST(the_correction_terms_are_what_they_are_defined_to_be) {
  // C1 = E[psi_c(r)] and C2 = E[psi_c(r)(y - mu)], computed here by an
  // independent summation with a different stopping rule, so a mistake in the
  // limit or the recurrence shows up as a disagreement.
  for (double mu : {0.05, 0.3, 1.0, 5.0, 40.0}) {
    const double c = 1.345;
    const double root = std::sqrt(mu);
    double c1 = 0.0, c2 = 0.0;
    // Straight from the definition: factorials in long double, no recurrence.
    for (int j = 0; j < 400; ++j) {
      long double logp = -static_cast<long double>(mu) +
                         static_cast<long double>(j) * std::log((long double)mu);
      for (int k = 2; k <= j; ++k) logp -= std::log((long double)k);
      const double p = static_cast<double>(std::exp(logp));
      const double d = static_cast<double>(j) - mu;
      double psi = d / root;
      psi = psi > c ? c : (psi < -c ? -c : psi);
      c1 += psi * p;
      c2 += psi * d * p;
    }
    check::close(glm_c1(mu, c), c1, 1e-10, "C1 agrees with its definition");
    check::close(glm_c2(mu, c), c2, 1e-10, "C2 agrees with its definition");
  }
}

TEST(with_no_clipping_the_corrections_take_their_known_values) {
  // C1 is E[r], which is zero, and C2 is E[r(y-mu)] = var/sqrt(mu) = sqrt(mu).
  // Known in closed form, so this pins the summation against arithmetic rather
  // than against another summation.
  for (double mu : {0.3, 2.0, 30.0}) {
    check::close(glm_c1(mu, 1e6), 0.0, 1e-9, "C1 is zero without clipping");
    check::close(glm_c2(mu, 1e6), std::sqrt(mu), 1e-7 * std::sqrt(mu),
                 "C2 is sqrt(mu) without clipping");
  }
}

TEST(the_glm_is_unbiased_on_clean_poisson_data) {
  // The point of the Fisher consistency correction: clipping one tail of an
  // asymmetric distribution would otherwise bias the estimate low, which is
  // exactly what the traditional methods do.
  for (double mu : {0.1, 0.3, 1.0, 3.0}) {
    double total = 0.0;
    const int trials = 40;
    for (int t = 0; t < trials; ++t) {
      const std::vector<double> v = poisson_sample(mu, 3000, 1000 + t);
      total += glm_background(v).mean;
    }
    const double average = total / trials;
    check::close(average, mu, 0.04 * mu,
                 "the fitted background is the true one, within four per cent");
  }
}

TEST(an_outlier_moves_the_mean_and_not_the_glm) {
  // A zinger, a hot pixel, or a neighbouring reflection leaking in. The plain
  // mean follows it; the robust estimate should barely notice.
  std::vector<double> v = poisson_sample(0.3, 500, 7);
  const double clean_mean = mean_of(v);
  const double clean_glm = glm_background(v).mean;
  v.push_back(5000.0);

  const double dirty_mean = mean_of(v);
  const double dirty_glm = glm_background(v).mean;

  check::is_true(dirty_mean > 5.0 * clean_mean,
                 "one outlier in five hundred wrecks the mean");
  check::is_true(std::fabs(dirty_glm - clean_glm) < 0.1 * clean_glm,
                 "and moves the glm by less than a tenth");
}

TEST(an_empty_background_is_zero_and_not_a_failure) {
  // Common on a weak dataset: every background pixel of a shoebox is zero.
  // theta is log(mu) and there is no log of zero, so this is answered rather
  // than iterated towards.
  const std::vector<double> v(400, 0.0);
  const BackgroundResult r = glm_background(v);
  check::is_true(r.valid && r.converged, "answered, not failed");
  check::close(r.mean, 0.0, 0.0, "and the answer is zero");
}

TEST(nothing_at_all_is_a_failure) {
  const BackgroundResult r = glm_background({});
  check::is_true(!r.valid, "no pixels is not a background");
}

}  // namespace mxi
