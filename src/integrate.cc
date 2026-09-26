#include "integrate.hh"
#include "refl.hh"

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

double resolution(const Crystal &crystal, int h, int k, int l) {
  // A rather than A_at: the static cell is what the column holds, even for a
  // scan-varying crystal.
  const Vec3 q = crystal.A * Vec3{static_cast<double>(h), static_cast<double>(k),
                                  static_cast<double>(l)};
  const double length = q.norm();
  if (!(length > 0.0)) return 0.0;
  return 1.0 / length;
}

double partiality(const Scan &scan, double phi, double zeta, double sigma_m,
                  std::int32_t z_first, std::int32_t z_last) {
  const double spread = Scan::radians(sigma_m) / std::fabs(zeta);
  if (!(spread > 0.0)) return 1.0;
  const double a = scan.phi_from_z(static_cast<double>(z_first));
  const double b = scan.phi_from_z(static_cast<double>(z_last));
  const double lo = std::fmin(a, b);
  const double hi = std::fmax(a, b);
  const double scale = std::sqrt(2.0) * spread;
  const double fraction =
      0.5 * (std::erf((hi - phi) / scale) - std::erf((lo - phi) / scale));
  return std::fmax(0.0, std::fmin(1.0, fraction));
}

double lorentz_polarization(const Beam &beam, const Goniometer &goniometer,
                            const Vec3 &s1) {
  const Vec3 s0 = beam.s0();
  const double s1_length = s1.norm();
  const double s0_length = s0.norm();
  if (!(s1_length > 0.0) || !(s0_length > 0.0)) return 0.0;
  const Vec3 m2 = goniometer.lab_axis();

  const double lorentz =
      std::fabs(s1.dot(m2.cross(s0))) / (s1_length * s0_length);

  const Vec3 u = s1 / s1_length;
  const double p = beam.polarization_fraction;
  const double to_normal = u.dot(beam.polarization_normal);
  const double to_beam = u.dot(s0 / s0_length);
  const double polarization = (1.0 - p) + (2.0 * p - 1.0) * to_normal * to_normal +
                              p * to_beam * to_beam;
  if (!(polarization > 0.0)) return 0.0;
  return lorentz / polarization;
}

std::int64_t summation_flags(const IntegratedReflection &r) {
  std::int64_t f = 0;
  if (r.n_foreground_bad > 0) {
    f |= flag::kForegroundIncludesBadPixels | flag::kFailedDuringSummation;
  } else if (r.valid) {
    f |= flag::kIntegratedSum;
  }
  if (r.n_background_bad > 0) f |= flag::kBackgroundIncludesBadPixels;
  return f;
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
    if ((m & shoebox_mask::kValid) == 0) {
      // Counted before it is skipped: a bad pixel keeps its region bit, so
      // this is where it is known which region lost a measurement.
      if (m & shoebox_mask::kForeground) ++out.n_foreground_bad;
      else if (m & shoebox_mask::kBackground) ++out.n_background_bad;
      continue;
    }
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

  // The observed centre of mass, over the foreground only and with the
  // background taken off. Over the whole box instead it is 0.15 to 0.19 pixels
  // away from DIALS', because the rim contributes its own noise with no signal
  // in it.
  //
  // Negative excesses are dropped rather than allowed to pull the centroid
  // about: a background-subtracted pixel below zero is noise, and on a weak
  // reflection there are enough of them to move the answer a long way.
  {
    double weight = 0.0, mf = 0.0, ms = 0.0, mz = 0.0;
    for (std::int32_t z = 0; z < box->nz(); ++z) {
      for (std::int32_t y = 0; y < box->ny(); ++y) {
        for (std::int32_t x = 0; x < box->nx(); ++x) {
          const std::size_t i = box->at(x, y, z);
          if ((box->mask[i] & shoebox_mask::kForeground) == 0) continue;
          if ((box->mask[i] & shoebox_mask::kValid) == 0) continue;
          const double w = static_cast<double>(box->data[i]) - background.mean;
          if (!(w > 0.0)) continue;
          weight += w;
          mf += w * (static_cast<double>(box->bbox[0] + x) + 0.5);
          ms += w * (static_cast<double>(box->bbox[2] + y) + 0.5);
          mz += w * (static_cast<double>(box->bbox[4] + z) + 0.5);
        }
      }
    }
    if (weight > 0.0) {
      out.centroid_valid = true;
      out.centroid_fast = mf / weight;
      out.centroid_slow = ms / weight;
      out.centroid_z = mz / weight;
      double sf = 0.0, ss = 0.0, sz = 0.0;
      for (std::int32_t z = 0; z < box->nz(); ++z) {
        for (std::int32_t y = 0; y < box->ny(); ++y) {
          for (std::int32_t x = 0; x < box->nx(); ++x) {
            const std::size_t i = box->at(x, y, z);
            if ((box->mask[i] & shoebox_mask::kForeground) == 0) continue;
            if ((box->mask[i] & shoebox_mask::kValid) == 0) continue;
            const double w = static_cast<double>(box->data[i]) - background.mean;
            if (!(w > 0.0)) continue;
            const double df =
                static_cast<double>(box->bbox[0] + x) + 0.5 - out.centroid_fast;
            const double ds =
                static_cast<double>(box->bbox[2] + y) + 0.5 - out.centroid_slow;
            const double dz =
                static_cast<double>(box->bbox[4] + z) + 0.5 - out.centroid_z;
            // The Poisson variance of this pixel's COUNT, not its excess over
            // the background. The centroid is sum(w x) / W with w = c - b, so
            // its variance is sum((x - xbar)^2 var(c)) / W^2, and var(c) is c:
            // the background photons are as noisy as the signal ones. Using w
            // here made the uncertainty of a weak reflection on a background
            // come out far too small, and independent of the background
            // altogether.
            const double count = static_cast<double>(box->data[i]);
            // And a count is not AT its pixel's centre but somewhere in the
            // pixel, uniformly, which is a variance of 1/12 of a pixel squared
            // on each axis -- and of an image squared along the scan. Without
            // it a spot whose counts share one coordinate, all on a single
            // image say, has every (x - xbar) zero and a variance of nothing,
            // and its pull came out at 1e12.
            constexpr double kQuantisation = 1.0 / 12.0;
            sf += count * (df * df + kQuantisation);
            ss += count * (ds * ds + kQuantisation);
            sz += count * (dz * dz + kQuantisation);
          }
        }
      }
      // The variance of the MEAN: the second moment of the distribution
      // divided by the weight in it, which is the standard quantity and goes
      // to zero as a spot gets stronger.
      //
      // Propagated from the Poisson noise on each count, so it grows with
      // the background and not only with the signal.
      out.centroid_variance_fast = sf / (weight * weight);
      out.centroid_variance_slow = ss / (weight * weight);
      out.centroid_variance_z = sz / (weight * weight);
    }
  }

  // The unclipped centre of mass. Every foreground pixel, with its excess
  // over the background whatever its sign, and the variance propagated from
  // the Poisson noise on each count plus the 1/12 of a pixel that a count's
  // position within its pixel adds. This is the estimator whose uncertainty
  // was checked against spots planted at known positions.
  {
    constexpr double kQuantisation = 1.0 / 12.0;
    double weight = 0.0, mf = 0.0, ms = 0.0, mz = 0.0;
    for (std::int32_t z = 0; z < box->nz(); ++z) {
      for (std::int32_t y = 0; y < box->ny(); ++y) {
        for (std::int32_t x = 0; x < box->nx(); ++x) {
          const std::size_t i = box->at(x, y, z);
          if ((box->mask[i] & shoebox_mask::kForeground) == 0) continue;
          if ((box->mask[i] & shoebox_mask::kValid) == 0) continue;
          const double w = static_cast<double>(box->data[i]) - background.mean;
          weight += w;
          mf += w * (static_cast<double>(box->bbox[0] + x) + 0.5);
          ms += w * (static_cast<double>(box->bbox[2] + y) + 0.5);
          mz += w * (static_cast<double>(box->bbox[4] + z) + 0.5);
        }
      }
    }
    // A centre exists only for a reflection that was DETECTED. The centre is
    // sum(w x) / W, and a W that is barely positive by noise throws it
    // anywhere: requiring only W > 0 put 117 of 17203 residuals more than
    // fifty pixels from their prediction, the worst at -31025 images, with a
    // median I/sigma of 0.07 among them. The variance was honest about it --
    // 1e21 square pixels -- but a min, a max, a plot or an unweighted mean is
    // wrecked by one such row.
    //
    // So W must stand clear of its own noise, whose variance is the sum of the
    // counts: the same quantity the instability comes from, rather than an
    // intensity cut chosen for the purpose.
    double counts = 0.0;
    for (std::size_t i = 0; i < box->size(); ++i) {
      if ((box->mask[i] & shoebox_mask::kForeground) == 0) continue;
      if ((box->mask[i] & shoebox_mask::kValid) == 0) continue;
      counts += static_cast<double>(box->data[i]);
    }
    const bool detected =
        weight > 0.0 && counts > 0.0 &&
        weight >= options.least_centroid_significance * std::sqrt(counts);
    if (detected) {
      out.unbiased_fast = mf / weight;
      out.unbiased_slow = ms / weight;
      out.unbiased_z = mz / weight;
      double vf = 0.0, vs = 0.0, vz = 0.0;
      for (std::int32_t z = 0; z < box->nz(); ++z) {
        for (std::int32_t y = 0; y < box->ny(); ++y) {
          for (std::int32_t x = 0; x < box->nx(); ++x) {
            const std::size_t i = box->at(x, y, z);
            if ((box->mask[i] & shoebox_mask::kForeground) == 0) continue;
            if ((box->mask[i] & shoebox_mask::kValid) == 0) continue;
            const double count = static_cast<double>(box->data[i]);
            const double df =
                static_cast<double>(box->bbox[0] + x) + 0.5 - out.unbiased_fast;
            const double ds =
                static_cast<double>(box->bbox[2] + y) + 0.5 - out.unbiased_slow;
            const double dz =
                static_cast<double>(box->bbox[4] + z) + 0.5 - out.unbiased_z;
            vf += count * (df * df + kQuantisation);
            vs += count * (ds * ds + kQuantisation);
            vz += count * (dz * dz + kQuantisation);
          }
        }
      }
      out.unbiased_variance_fast = vf / (weight * weight);
      out.unbiased_variance_slow = vs / (weight * weight);
      out.unbiased_variance_z = vz / (weight * weight);
      // And inside its own box. With every weight positive a centre of mass
      // cannot leave the region it was taken over; with signed weights it
      // can, and one that has is describing the noise and not the spot.
      const bool inside =
          out.unbiased_fast >= box->bbox[0] && out.unbiased_fast <= box->bbox[1] &&
          out.unbiased_slow >= box->bbox[2] && out.unbiased_slow <= box->bbox[3] &&
          out.unbiased_z >= box->bbox[4] && out.unbiased_z <= box->bbox[5];
      out.unbiased_valid = inside;
    }
  }

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
