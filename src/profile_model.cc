#include "profile_model.h"

#include <cmath>

namespace mxi {

bool spot_angular_variance(const Experiment &e, const Shoebox &box, const Vec3 &s1,
                           double *variance) {
  if (box.panel < 0 || static_cast<std::size_t>(box.panel) >= e.detector.size()) {
    return false;
  }
  const Panel &p = e.detector[static_cast<std::size_t>(box.panel)];
  const double length = s1.norm();
  if (!(length > 0.0)) return false;
  const Vec3 direction = s1 / length;

  double weight = 0.0;
  double weighted = 0.0;
  for (std::int32_t z = 0; z < box.nz(); ++z) {
    for (std::int32_t y = 0; y < box.ny(); ++y) {
      for (std::int32_t x = 0; x < box.nx(); ++x) {
        const std::size_t at = box.at(x, y, z);
        if ((box.mask[at] & shoebox_mask::kForeground) == 0) continue;
        const double count = static_cast<double>(box.data[at]) -
                             static_cast<double>(box.background[at]);
        if (!(count > 0.0)) continue;

        // The direction of the ray that would have landed in the middle of
        // this pixel. Millimetres straight from the pixel size, with no
        // parallax: the correction moves where a ray of a given direction is
        // recorded, and undoing it here would be asking where the ray came
        // from, which is a different question and gives a spread 15 per cent
        // smaller.
        const double mm_fast =
            (static_cast<double>(box.bbox[0] + x) + 0.5) * p.pixel_size[0];
        const double mm_slow =
            (static_cast<double>(box.bbox[2] + y) + 0.5) * p.pixel_size[1];
        const Vec3 lab = p.lab_coord_mm(mm_fast, mm_slow);
        const double n = lab.norm();
        if (!(n > 0.0)) continue;

        double cosine = (lab / n).dot(direction);
        cosine = std::fmax(-1.0, std::fmin(1.0, cosine));
        const double angle = std::acos(cosine);
        weight += count;
        weighted += count * angle * angle;
      }
    }
  }
  // One count has no spread, and the sample variance of a single observation
  // is not zero, it is undefined.
  if (!(weight > 1.0)) return false;
  *variance = weighted / (weight - 1.0);
  return true;
}

double recorded_fraction(double delta, double zeta, double sigma,
                         double oscillation) {
  const double spread = sigma / std::fabs(zeta);
  if (!(spread > 0.0) || !(oscillation > 0.0)) return 0.0;
  const double scale = std::sqrt(2.0) * spread;
  const double upper = std::erf((delta + 0.5 * oscillation) / scale);
  const double lower = std::erf((delta - 0.5 * oscillation) / scale);
  // A density in delta: it integrates to one over delta, which is what makes
  // the product of these a likelihood rather than a product of fractions.
  return (upper - lower) / (2.0 * oscillation);
}

std::vector<RangeSample> range_samples(const Experiment &e, const Shoebox &box,
                                       double phi_calculated, double zeta) {
  std::vector<RangeSample> out;
  const double width = Scan::radians(e.scan.osc_width);
  for (std::int32_t z = 0; z < box.nz(); ++z) {
    double counts = 0.0;
    for (std::int32_t y = 0; y < box.ny(); ++y) {
      for (std::int32_t x = 0; x < box.nx(); ++x) {
        const std::size_t at = box.at(x, y, z);
        if ((box.mask[at] & shoebox_mask::kForeground) == 0) continue;
        counts += static_cast<double>(box.data[at]) -
                  static_cast<double>(box.background[at]);
      }
    }
    if (!(counts > 0.0)) continue;
    // The centre of this image, in radians.
    const double image = static_cast<double>(box.bbox[4] + z) + 0.5;
    const double centre = Scan::radians(e.scan.osc_start) +
                          (image - static_cast<double>(e.scan.z_offset)) * width;
    out.push_back({phi_calculated - centre, zeta});
  }
  return out;
}

double reflecting_range(const std::vector<RangeSample> &samples, double oscillation,
                        double min_zeta) {
  std::vector<RangeSample> used;
  for (const RangeSample &s : samples) {
    if (std::fabs(s.zeta) >= min_zeta) used.push_back(s);
  }
  if (used.size() < 2) return 0.0;

  const auto negative_log_likelihood = [&](double sigma) {
    double total = 0.0;
    for (const RangeSample &s : used) {
      const double r = recorded_fraction(s.delta, s.zeta, sigma, oscillation);
      // A sample the model says is impossible would take the whole sum to
      // infinity; floored so that one such sample cannot decide the answer on
      // its own.
      total -= std::log(std::fmax(r, 1e-300));
    }
    return total;
  };

  // Golden section on a one-dimensional problem, which needs no derivative and
  // cannot overshoot. Bracketed from a hundredth of a degree to two degrees:
  // outside that range the answer is not a mosaicity.
  double low = Scan::radians(0.001);
  double high = Scan::radians(2.0);
  const double phi = 0.5 * (std::sqrt(5.0) - 1.0);
  double b = high - phi * (high - low);
  double c = low + phi * (high - low);
  double fb = negative_log_likelihood(b);
  double fc = negative_log_likelihood(c);
  for (int i = 0; i < 200 && (high - low) > 1e-12; ++i) {
    if (fb < fc) {
      high = c;
      c = b;
      fc = fb;
      b = high - phi * (high - low);
      fb = negative_log_likelihood(b);
    } else {
      low = b;
      b = c;
      fb = fc;
      c = low + phi * (high - low);
      fc = negative_log_likelihood(c);
    }
  }
  return Scan::degrees(0.5 * (low + high));
}

double beam_divergence(const Experiment &e, const std::vector<Shoebox> &boxes,
                       const std::vector<Vec3> &s1, std::size_t *n_used) {
  double total = 0.0;
  std::size_t used = 0;
  const std::size_t n = std::min(boxes.size(), s1.size());
  for (std::size_t i = 0; i < n; ++i) {
    double variance = 0.0;
    if (!spot_angular_variance(e, boxes[i], s1[i], &variance)) continue;
    total += variance;
    ++used;
  }
  if (n_used != nullptr) *n_used = used;
  if (used == 0) return 0.0;
  return Scan::degrees(std::sqrt(total / static_cast<double>(used)));
}

}  // namespace mxi
