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
        // Any pixel the spot finder marked, not only the foreground, and the
        // background is NOT subtracted. Kabsch section 3.1 step (v) says to
        // subtract it; DIALS does not, and carries a note saying so. Matching
        // DIALS is the point here, and on a table out of dials.find_spots the
        // background is zero in any case.
        if (box.mask[at] == 0) continue;
        const double count = static_cast<double>(box.data[at]);
        if (!(count > 0.0)) continue;

        // The direction of the ray that would have landed in the middle of
        // this pixel, through the SAME px-to-mm mapping the s1 it is compared
        // against was built with, parallax and all. Mixing the two conventions
        // displaces every pixel radially by about half a pixel, which showed
        // up as a systematic 0.44 sigma offset in eps2 and nowhere else.
        const auto mm = p.px_to_mm(static_cast<double>(box.bbox[0] + x) + 0.5,
                                   static_cast<double>(box.bbox[2] + y) + 0.5);
        const Vec3 lab = p.lab_coord_mm(mm.first, mm.second);
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

double compute_zeta(const Experiment &e, const Vec3 &s1) {
  const Vec3 axis = e.goniometer.lab_axis();
  const Vec3 e1 = s1.cross(e.beam.s0());
  const double length = e1.norm();
  if (!(length > 0.0)) return 0.0;
  return axis.dot(e1 / length);
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
  const std::uint8_t wanted = shoebox_mask::kValid | shoebox_mask::kForeground;
  for (std::int32_t z = 0; z < box.nz(); ++z) {
    // An image counts if it holds any pixel the spot finder marked as valid
    // foreground, whatever its value. Requiring positive counts as well is a
    // different criterion; on these tables it selects the same images, but it
    // is not what DIALS asks.
    bool marked = false;
    for (std::int32_t y = 0; y < box.ny() && !marked; ++y) {
      for (std::int32_t x = 0; x < box.nx() && !marked; ++x) {
        if (box.mask[box.at(x, y, z)] == wanted) marked = true;
      }
    }
    if (!marked) continue;
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

KabschFrame kabsch_frame(const Experiment &e, const Vec3 &s1) {
  KabschFrame out;
  const Vec3 s0 = e.beam.s0();
  const Vec3 cross = s1.cross(s0);
  const double length = cross.norm();
  if (!(length > 0.0) || !(s1.norm() > 0.0)) return out;
  out.e1 = cross / length;
  const Vec3 second = s1.cross(out.e1);
  const double second_length = second.norm();
  if (!(second_length > 0.0)) return out;
  out.e2 = second / second_length;
  const Vec3 third = s1 + s0;
  const double third_length = third.norm();
  if (!(third_length > 0.0)) return out;
  out.e3 = third / third_length;
  out.s1 = s1;
  out.zeta = e.goniometer.lab_axis().dot(out.e1);
  out.valid = true;
  return out;
}

Epsilon epsilon_of(const Experiment &e, const KabschFrame &frame, const Panel &p,
                   double px_fast, double px_slow, double phi_image,
                   double phi_calculated) {
  Epsilon out;
  if (!frame.valid) return out;
  const double length = frame.s1.norm();
  const auto mm = p.px_to_mm(px_fast, px_slow);
  const Vec3 lab = p.lab_coord_mm(mm.first, mm.second);
  const double lab_length = lab.norm();
  if (!(lab_length > 0.0)) return out;
  // The diffracted beam that would have gone through this pixel, of the same
  // length as the reflection's own: the scattering is elastic.
  const Vec3 s_prime = lab * (length / lab_length);
  const Vec3 difference = s_prime - frame.s1;
  // Divided by |s1| to turn a distance in reciprocal space into an angle, and
  // reported in degrees, as Kabsch writes it.
  out.e1 = Scan::degrees(frame.e1.dot(difference) / length);
  out.e2 = Scan::degrees(frame.e2.dot(difference) / length);
  out.e3 = Scan::degrees(frame.zeta * (phi_image - phi_calculated));
  (void)e;
  return out;
}

Capture capture_fractions(const Experiment &e, const std::vector<Shoebox> &boxes,
                          const std::vector<Vec3> &s1,
                          const std::vector<double> &phi_calculated,
                          double sigma_d, double sigma_m) {
  Capture out;
  double inside[4] = {0.0, 0.0, 0.0, 0.0};
  double in_detector[4] = {0.0, 0.0, 0.0, 0.0};
  double in_rotation[4] = {0.0, 0.0, 0.0, 0.0};
  double total = 0.0;
  const double width = Scan::radians(e.scan.osc_width);
  const std::size_t n = std::min(boxes.size(), s1.size());

  for (std::size_t i = 0; i < n && i < phi_calculated.size(); ++i) {
    const Shoebox &box = boxes[i];
    if (box.panel < 0 || static_cast<std::size_t>(box.panel) >= e.detector.size()) {
      continue;
    }
    const Panel &p = e.detector[static_cast<std::size_t>(box.panel)];
    const KabschFrame frame = kabsch_frame(e, s1[i]);
    if (!frame.valid) continue;
    ++out.n_spots;

    for (std::int32_t z = 0; z < box.nz(); ++z) {
      const double image = static_cast<double>(box.bbox[4] + z) + 0.5;
      const double phi = Scan::radians(e.scan.osc_start) +
                         (image - static_cast<double>(e.scan.z_offset)) * width;
      for (std::int32_t y = 0; y < box.ny(); ++y) {
        for (std::int32_t x = 0; x < box.nx(); ++x) {
          const std::size_t at = box.at(x, y, z);
          if (box.mask[at] == 0) continue;
          const double count = static_cast<double>(box.data[at]) -
                               static_cast<double>(box.background[at]);
          if (!(count > 0.0)) continue;
          const Epsilon eps =
              epsilon_of(e, frame, p, static_cast<double>(box.bbox[0] + x) + 0.5,
                         static_cast<double>(box.bbox[2] + y) + 0.5, phi,
                         phi_calculated[i]);
          total += count;
          for (int k = 0; k < 4; ++k) {
            const double limit = static_cast<double>(k + 1);
            const bool on_detector = std::fabs(eps.e1) <= limit * sigma_d &&
                                     std::fabs(eps.e2) <= limit * sigma_d;
            const bool in_rotation_range = std::fabs(eps.e3) <= limit * sigma_m;
            if (on_detector) in_detector[k] += count;
            if (in_rotation_range) in_rotation[k] += count;
            if (on_detector && in_rotation_range) inside[k] += count;
          }
        }
      }
    }
  }
  out.counts = total;
  if (total > 0.0) {
    for (int k = 0; k < 4; ++k) {
      out.fraction[k] = inside[k] / total;
      out.detector[k] = in_detector[k] / total;
      out.rotation[k] = in_rotation[k] / total;
    }
  }
  return out;
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
