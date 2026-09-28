#include <cmath>

#include "../src/scale_model.hh"
#include "check.hh"

namespace mxi {

namespace {

//: Gauss-Legendre nodes and weights on [-1, 1], by Newton's method on P_n.
void gauss_legendre(int n, std::vector<double> *x, std::vector<double> *w) {
  const double pi = std::acos(-1.0);
  x->assign(static_cast<std::size_t>(n), 0.0);
  w->assign(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    double z = std::cos(pi * (i + 0.75) / (n + 0.5)), pp = 0.0;
    for (int it = 0; it < 100; ++it) {
      double p1 = 1.0, p2 = 0.0;
      for (int j = 1; j <= n; ++j) {
        const double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
      }
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      const double dz = p1 / pp;
      z -= dz;
      if (std::abs(dz) < 1e-15)
        break;
    }
    (*x)[static_cast<std::size_t>(i)] = z;
    (*w)[static_cast<std::size_t>(i)] = 2.0 / ((1.0 - z * z) * pp * pp);
  }
}

} // namespace

TEST(degree_one_harmonics_are_the_coordinates) {
  // Without the Condon-Shortley sign, Y_1 at a unit vector is
  // sqrt(3 / 4 pi) (y, z, x) for m = -1, 0, 1: exact, and independent of the
  // recurrences the code uses.
  const double k = std::sqrt(3.0 / (4.0 * std::acos(-1.0)));
  const Vec3 u = Vec3{0.3, -0.5, 0.81}.normalized();
  std::vector<double> y;
  real_spherical_harmonics(1, u, &y);
  check::equal(static_cast<long long>(y.size()), 3, "three of degree one");
  check::close(y[0], k * u.y, 1e-15, "m = -1 is y");
  check::close(y[1], k * u.z, 1e-15, "m = 0 is z");
  check::close(y[2], k * u.x, 1e-15, "m = 1 is x");
}

TEST(the_harmonics_to_degree_four_are_orthonormal) {
  // Gauss-Legendre in cos(theta), 16 points, and 32 even steps in phi,
  // integrate products of these harmonics exactly: the Gram matrix of all 24
  // must be the identity to rounding. Any wrong normalisation or recurrence
  // term shows as an entry off by far more.
  const int lmax = 4;
  const std::size_t n = harmonic_count(lmax);
  check::equal(static_cast<long long>(n), 24,
               "24 to degree four, as the paper counts");
  std::vector<double> xs, ws;
  gauss_legendre(16, &xs, &ws);
  std::vector<double> gram(n * n, 0.0), y;
  const double pi = std::acos(-1.0);
  for (std::size_t i = 0; i < xs.size(); ++i) {
    const double z = xs[i], s = std::sqrt(1.0 - z * z);
    for (int j = 0; j < 32; ++j) {
      const double phi = 2.0 * pi * j / 32.0;
      real_spherical_harmonics(lmax, {s * std::cos(phi), s * std::sin(phi), z},
                               &y);
      const double w = ws[i] * 2.0 * pi / 32.0;
      for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = 0; b < n; ++b)
          gram[a * n + b] += w * y[a] * y[b];
    }
  }
  double worst = 0.0;
  for (std::size_t a = 0; a < n; ++a)
    for (std::size_t b = 0; b < n; ++b)
      worst =
          std::fmax(worst, std::abs(gram[a * n + b] - (a == b ? 1.0 : 0.0)));
  check::is_true(worst < 1e-12,
                 "orthonormal: worst entry off by " + std::to_string(worst));
}

TEST(the_scale_models_gradient_is_its_derivative) {
  ScaleModel model({6, 5, 4});
  for (std::size_t i = 0; i < model.size(); ++i)
    model.parameters[i] += 0.05 * std::sin(1.7 * static_cast<double>(i) + 0.3);
  ScaleObservation o;
  o.rotation = 0.37;
  o.time = 0.61;
  o.inv_2d2 = 1.0 / (2.0 * 2.1 * 2.1);
  std::vector<double> y0, y1;
  real_spherical_harmonics(4, Vec3{0.2, 0.7, -0.4}.normalized(), &y0);
  real_spherical_harmonics(4, Vec3{-0.1, 0.1, 1.0}.normalized(), &y1);
  for (std::size_t k = 0; k < y0.size(); ++k)
    o.absorption.push_back(0.5 * (y0[k] + y1[k]));
  std::vector<std::pair<std::size_t, double>> grad;
  model.inverse_scale(o, &grad);
  std::vector<double> analytic(model.size(), 0.0);
  for (const auto &[i, v] : grad)
    analytic[i] += v;
  double worst = 0.0;
  for (std::size_t i = 0; i < model.size(); ++i) {
    ScaleModel up = model, down = model;
    const double h = 1e-6;
    up.parameters[i] += h;
    down.parameters[i] -= h;
    const double numeric =
        (up.inverse_scale(o) - down.inverse_scale(o)) / (2.0 * h);
    worst = std::fmax(worst, std::abs(numeric - analytic[i]));
  }
  check::is_true(worst < 1e-8, "every derivative, to " + std::to_string(worst));
}

TEST(the_default_shape_follows_dials) {
  // 30 degrees: six and five, as dials.scale's log for such a sweep shows, and
  // no absorption surface. 360: the paper's 26, 20 and 24 -- 70 parameters.
  const ScaleModelShape narrow = default_shape(30.0);
  check::equal(static_cast<long long>(narrow.scale_points), 6,
               "30 degrees: six scale");
  check::equal(static_cast<long long>(narrow.decay_points), 5, "five decay");
  check::equal(static_cast<long long>(narrow.lmax), 0, "no absorption");
  const ScaleModelShape full = default_shape(360.0);
  check::equal(static_cast<long long>(ScaleModel(full).size()), 70,
               "360 degrees: 70");
  check::equal(static_cast<long long>(full.scale_points), 26, "26 scale");
  check::equal(static_cast<long long>(full.decay_points), 20, "20 decay");
}

TEST(normalising_the_scale_changes_no_ratio_between_observations) {
  ScaleModel model({6, 0, 0});
  for (std::size_t i = 0; i < 6; ++i)
    model.parameters[i] = 1.5 + 0.1 * static_cast<double>(i);
  ScaleObservation a, b;
  a.rotation = 0.1;
  b.rotation = 0.8;
  const double before = model.inverse_scale(a) / model.inverse_scale(b);
  model.normalise();
  double mean = 0.0;
  for (std::size_t i = 0; i < 6; ++i)
    mean += model.parameters[i] / 6.0;
  check::close(mean, 1.0, 1e-15, "the scale's mean is one");
  check::close(model.inverse_scale(a) / model.inverse_scale(b), before, 1e-14,
               "and every ratio, which is all a fit sees, is unchanged");
}

} // namespace mxi
