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
