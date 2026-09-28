#include "scale_model.hh"

#include <algorithm>
#include <cmath>

#include "derivatives.hh"

namespace mxi {

std::size_t harmonic_count(int lmax) {
  return lmax > 0 ? static_cast<std::size_t>(lmax * (lmax + 2)) : 0;
}

void real_spherical_harmonics(int lmax, const Vec3 &u,
                              std::vector<double> *out) {
  out->assign(harmonic_count(lmax), 0.0);
  if (lmax <= 0)
    return;
  const double pi = std::acos(-1.0);
  const double x = std::fmax(-1.0, std::fmin(1.0, u.z));   // cos theta
  const double s = std::sqrt(std::fmax(0.0, 1.0 - x * x)); // sin theta
  const double phi = std::atan2(u.y, u.x);
  // Associated Legendre functions P_l^m(x) without the (-1)^m, by the stable
  // recurrences: P_m^m, then P_{m+1}^m, then upward in l.
  std::vector<std::vector<double>> P(
      static_cast<std::size_t>(lmax + 1),
      std::vector<double>(static_cast<std::size_t>(lmax + 1), 0.0));
  P[0][0] = 1.0;
  for (int m = 1; m <= lmax; ++m)
    P[m][m] = P[m - 1][m - 1] * static_cast<double>(2 * m - 1) * s;
  for (int m = 0; m < lmax; ++m)
    P[m + 1][m] = x * static_cast<double>(2 * m + 1) * P[m][m];
  for (int m = 0; m <= lmax; ++m)
    for (int l = m + 2; l <= lmax; ++l)
      P[l][m] = (static_cast<double>(2 * l - 1) * x * P[l - 1][m] -
                 static_cast<double>(l + m - 1) * P[l - 2][m]) /
                static_cast<double>(l - m);
  std::size_t k = 0;
  for (int l = 1; l <= lmax; ++l) {
    for (int m = -l; m <= l; ++m, ++k) {
      const int a = std::abs(m);
      // N = sqrt((2l+1)/(4 pi) (l-|m|)!/(l+|m|)!), the ratio by a product.
      double ratio = 1.0;
      for (int j = l - a + 1; j <= l + a; ++j)
        ratio /= static_cast<double>(j);
      const double norm =
          std::sqrt(static_cast<double>(2 * l + 1) / (4.0 * pi) * ratio);
      const double p = P[l][a];
      if (m == 0)
        (*out)[k] = norm * p;
      else if (m > 0)
        (*out)[k] = std::sqrt(2.0) * norm * p * std::cos(a * phi);
      else
        (*out)[k] = std::sqrt(2.0) * norm * p * std::sin(a * phi);
    }
  }
}

ScaleModelShape default_shape(double degrees) {
  ScaleModelShape shape;
  degrees = std::abs(degrees);
  if (degrees < 90.0) {
    // A quarter and a third of the sweep, as DIALS' narrow-sweep defaults:
    // six and five points on 30 degrees, as its log for such a sweep shows.
    shape.scale_points = 6;
    shape.decay_points = 5;
  } else {
    shape.scale_points =
        static_cast<std::size_t>(std::ceil(degrees / 15.0)) + 2;
    shape.decay_points =
        static_cast<std::size_t>(std::ceil(degrees / 20.0)) + 2;
  }
  shape.lmax = degrees >= 60.0 ? 4 : 0;
  return shape;
}

ScaleModel::ScaleModel(const ScaleModelShape &shape) : shape_(shape) {
  shape_.scale_points = std::max<std::size_t>(1, shape_.scale_points);
  parameters.assign(shape_.scale_points + shape_.decay_points +
                        harmonic_count(shape_.lmax),
                    0.0);
  for (std::size_t i = 0; i < shape_.scale_points; ++i)
    parameters[i] = 1.0;
}

double ScaleModel::inverse_scale(
    const ScaleObservation &o,
    std::vector<std::pair<std::size_t, double>> *gradient) const {
  const SplineWeights cw = spline_weights(shape_.scale_points, o.rotation);
  double c = 0.0;
  for (std::size_t k = 0; k < cw.count; ++k)
    c += cw.weight[k] * parameters[cw.index[k]];

  double b = 0.0;
  SplineWeights bw;
  if (shape_.decay_points > 0) {
    bw = spline_weights(shape_.decay_points, o.time);
    for (std::size_t k = 0; k < bw.count; ++k)
      b += bw.weight[k] * parameters[first_decay() + bw.index[k]];
  }
  const double t = std::exp(b * o.inv_2d2);

  double sa = 1.0;
  const std::size_t na = harmonic_count(shape_.lmax);
  for (std::size_t k = 0; k < na && k < o.absorption.size(); ++k)
    sa += parameters[first_absorption() + k] * o.absorption[k];

  const double g = c * t * sa;
  if (gradient) {
    gradient->clear();
    for (std::size_t k = 0; k < cw.count; ++k)
      gradient->emplace_back(cw.index[k], cw.weight[k] * t * sa);
    if (shape_.decay_points > 0)
      for (std::size_t k = 0; k < bw.count; ++k)
        gradient->emplace_back(first_decay() + bw.index[k],
                               g * o.inv_2d2 * bw.weight[k]);
    for (std::size_t k = 0; k < na && k < o.absorption.size(); ++k)
      gradient->emplace_back(first_absorption() + k, c * t * o.absorption[k]);
  }
  return g;
}

void ScaleModel::normalise() {
  double mean = 0.0;
  for (std::size_t i = 0; i < shape_.scale_points; ++i)
    mean += parameters[i];
  mean /= static_cast<double>(shape_.scale_points);
  if (mean > 0.0)
    for (std::size_t i = 0; i < shape_.scale_points; ++i)
      parameters[i] /= mean;
  if (shape_.decay_points > 0) {
    double b = 0.0;
    for (std::size_t i = 0; i < shape_.decay_points; ++i)
      b += parameters[first_decay() + i];
    b /= static_cast<double>(shape_.decay_points);
    for (std::size_t i = 0; i < shape_.decay_points; ++i)
      parameters[first_decay() + i] -= b;
  }
}

} // namespace mxi
