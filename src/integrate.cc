#include "integrate.h"

#include <cmath>

namespace mxi {

double quantum_efficiency(const Panel &panel, const Vec3 &s1) {
  if (!(panel.mu > 0.0) || !(panel.thickness > 0.0)) return 1.0;
  const double length = s1.norm();
  if (!(length > 0.0)) return 1.0;
  const Vec3 normal = panel.fast.cross(panel.slow).normalized();
  const double cosine = std::fabs((s1 / length).dot(normal));
  if (!(cosine > 0.0)) return 1.0;
  return 1.0 - std::exp(-panel.mu * panel.thickness / cosine);
}

IntegratedReflection integrate_shoebox(Shoebox *box,
                                       const IntegrateOptions &options) {
  IntegratedReflection out;
  if (box == nullptr || box->data.size() != box->size()) return out;
  if (box->mask.size() != box->size()) return out;

  std::vector<double> background_values;
  double foreground_sum = 0.0;
  std::size_t n_foreground = 0;
  std::size_t n_valid = 0;
  for (std::size_t i = 0; i < box->size(); ++i) {
    const std::uint8_t m = box->mask[i];
    if ((m & shoebox_mask::kValid) == 0) continue;
    ++n_valid;
    if (m & shoebox_mask::kForeground) {
      foreground_sum += static_cast<double>(box->data[i]);
      ++n_foreground;
    } else if (m & shoebox_mask::kBackground) {
      background_values.push_back(static_cast<double>(box->data[i]));
    }
  }
  out.n_foreground = n_foreground;
  out.n_background = background_values.size();
  out.n_valid = n_valid;

  if (background_values.size() < options.min_background) {
    out.too_few_background = true;
    return out;
  }
  const BackgroundResult background =
      glm_background(background_values, options.background);
  if (!background.valid) {
    out.background_failed = true;
    return out;
  }

  const double m = static_cast<double>(n_foreground);
  const double n = static_cast<double>(background_values.size());
  out.background_mean = background.mean;
  out.background_sum = background.mean * m;
  // The uncertainty in the background estimate, propagated to the foreground:
  // the (m/n) I_bg term of Leslie equation 11, written as a variance on the
  // summed background so it can be reported on its own.
  out.background_sum_variance = options.gain * (m / n) * out.background_sum;

  out.intensity = foreground_sum - out.background_sum;
  // Leslie (11): G [ I + I_bg + (m/n) I_bg ]. The first two together are the
  // Poisson noise of what was actually counted in the foreground, so they are
  // the raw foreground sum; writing it that way keeps the variance positive
  // when the intensity is negative, which happens for weak reflections and is
  // not an error.
  out.variance = options.gain * foreground_sum + out.background_sum_variance;

  // The fitted background, left in the shoebox so a saved one carries what was
  // subtracted from it.
  box->background.assign(box->size(), static_cast<float>(background.mean));
  out.valid = true;
  return out;
}

}  // namespace mxi
