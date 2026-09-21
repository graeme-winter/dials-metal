#include "background.h"

#include <algorithm>
#include <cmath>

namespace mxi {

namespace {

double huber(double r, double c) { return r > c ? c : (r < -c ? -c : r); }

//: Where to stop summing a Poisson. Ten standard deviations past the mean and
//: thirty besides, which covers the small-mean case where the distribution is
//: nearly all at zero and the standard deviation is no guide.
int summation_limit(double mu) {
  const double limit = mu + 10.0 * std::sqrt(std::fmax(mu, 0.0)) + 30.0;
  return static_cast<int>(std::ceil(limit));
}

}  // namespace

double glm_c1(double mu, double tuning) {
  if (!(mu > 0.0)) return 0.0;
  const double root = std::sqrt(mu);
  const int top = summation_limit(mu);
  // Poisson by recurrence from P(0) = exp(-mu): no factorials, no overflow.
  double p = std::exp(-mu);
  double total = 0.0;
  for (int j = 0; j <= top; ++j) {
    total += huber((static_cast<double>(j) - mu) / root, tuning) * p;
    p *= mu / static_cast<double>(j + 1);
  }
  return total;
}

double glm_c2(double mu, double tuning) {
  if (!(mu > 0.0)) return 0.0;
  const double root = std::sqrt(mu);
  const int top = summation_limit(mu);
  double p = std::exp(-mu);
  double total = 0.0;
  for (int j = 0; j <= top; ++j) {
    const double d = static_cast<double>(j) - mu;
    total += huber(d / root, tuning) * d * p;
    p *= mu / static_cast<double>(j + 1);
  }
  return total;
}

BackgroundResult glm_background(const std::vector<double> &values,
                                const BackgroundOptions &options) {
  BackgroundResult out;
  out.n = values.size();
  if (values.empty()) return out;

  double sum = 0.0;
  for (double v : values) sum += v;
  const double n = static_cast<double>(values.size());
  const double mean = sum / n;

  // All zero is a real answer at this background level and the iteration
  // cannot represent it -- theta is log(mu) and there is no log of zero. Said
  // here rather than returned as a failure, because a shoebox whose background
  // pixels are all empty is common on a weak dataset.
  if (!(mean > 0.0)) {
    out.valid = true;
    out.mean = 0.0;
    out.converged = true;
    return out;
  }

  double theta = std::log(mean);
  for (int i = 0; i < options.max_iterations; ++i) {
    const double mu = std::exp(theta);
    const double root = std::sqrt(mu);
    double score = 0.0;
    for (double v : values) score += huber((v - mu) / root, options.tuning);
    const double c1 = glm_c1(mu, options.tuning);
    const double c2 = glm_c2(mu, options.tuning);
    if (!(std::fabs(c2) > 0.0)) break;
    const double step = (score - n * c1) / (n * c2);
    theta += step;
    out.iterations = i + 1;
    if (std::fabs(step) < options.tolerance) {
      out.converged = true;
      break;
    }
    // A background cannot be enormous and the iteration should not wander
    // there looking for one; without this a shoebox of nothing but outliers
    // can take theta off to infinity.
    if (theta > 30.0 || theta < -60.0) break;
  }
  out.valid = true;
  out.mean = std::exp(theta);
  return out;
}

}  // namespace mxi
