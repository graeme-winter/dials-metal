#include "forward.hh"

#include <cmath>

namespace mxi {

std::vector<double> render_shoebox(const Experiment &e, const Shoebox &box,
                                   const Vec3 &s1, double phi_calculated,
                                   const ForwardOptions &options) {
  std::vector<double> out;
  if (box.panel < 0 ||
      static_cast<std::size_t>(box.panel) >= e.detector.size()) {
    return out;
  }
  const Panel &p = e.detector[static_cast<std::size_t>(box.panel)];
  const KabschFrame frame = kabsch_frame(e, s1);
  if (!frame.valid)
    return out;
  if (!(options.sigma_d > 0.0) || !(options.sigma_m > 0.0))
    return out;
  const int subdivisions = std::max(1, options.subdivisions);

  out.assign(box.size(), 0.0);
  const double width = Scan::radians(e.scan.osc_width);
  const double share = 1.0 / static_cast<double>(subdivisions * subdivisions);
  const double length = s1.norm();

  // Depths at which a photon is absorbed, and their weights: exponential, cut
  // off at the back of the sensor. Sampled at the midpoints of equal-
  // probability intervals, so a handful of samples carries the distribution
  // rather than its ends.
  std::vector<double> depth, depth_weight;
  if (options.sensor && p.mu > 0.0 && p.thickness > 0.0) {
    const double transmitted = std::exp(-p.mu * p.thickness);
    const double survive = 1.0 - transmitted;
    const int n = std::max(1, options.depth_samples);
    if (n == 1) {
      // One sample means every photon at the MEAN depth: the shift with none
      // of the spread, which is exactly what the parallax correction does.
      // Asking for one sample and silently getting no sensor at all, which is
      // what this did first, makes the two indistinguishable in the output.
      depth.push_back(1.0 / p.mu - p.thickness * transmitted / survive);
      depth_weight.push_back(1.0);
    } else {
      for (int i = 0; i < n; ++i) {
        const double q =
            (static_cast<double>(i) + 0.5) / static_cast<double>(n);
        depth.push_back(-std::log(1.0 - q * survive) / p.mu);
        depth_weight.push_back(1.0 / static_cast<double>(n));
      }
    }
  } else {
    depth.push_back(0.0);
    depth_weight.push_back(1.0);
  }

  for (std::int32_t z = 0; z < box.nz(); ++z) {
    // The angular range this image covers, as an eps3 range.
    const double image = static_cast<double>(box.bbox[4] + z);
    const double phi_low =
        Scan::radians(e.scan.osc_start) +
        (image - static_cast<double>(e.scan.z_offset)) * width;
    double e3_low = Scan::degrees(frame.zeta * (phi_low - phi_calculated));
    double e3_high =
        Scan::degrees(frame.zeta * (phi_low + width - phi_calculated));
    if (e3_low > e3_high)
      std::swap(e3_low, e3_high);
    // The Gaussian integrated over the image, which is what the image records.
    const double root_two = std::sqrt(2.0);
    const double along_rotation =
        0.5 * (std::erf(e3_high / (root_two * options.sigma_m)) -
               std::erf(e3_low / (root_two * options.sigma_m)));
    if (!(along_rotation > 0.0))
      continue;

    for (std::int32_t y = 0; y < box.ny(); ++y) {
      for (std::int32_t x = 0; x < box.nx(); ++x) {
        // Subdivisions of the pixel this photon would be RECORDED in.
        for (int sy = 0; sy < subdivisions; ++sy) {
          for (int sx = 0; sx < subdivisions; ++sx) {
            const double px = static_cast<double>(box.bbox[0] + x) +
                              (static_cast<double>(sx) + 0.5) /
                                  static_cast<double>(subdivisions);
            const double py = static_cast<double>(box.bbox[2] + y) +
                              (static_cast<double>(sy) + 0.5) /
                                  static_cast<double>(subdivisions);
            const double mm_fast = px * p.pixel_size[0];
            const double mm_slow = py * p.pixel_size[1];

            for (std::size_t d = 0; d < depth.size(); ++d) {
              // Where the ray recorded here ENTERED the sensor: back along its
              // own direction by the depth it travelled. The direction is
              // taken from the recorded point, which is accurate to well
              // inside a pixel over a depth of half a millimetre.
              const Vec3 lab = p.lab_coord_mm(mm_fast, mm_slow);
              const double lab_length = lab.norm();
              if (!(lab_length > 0.0))
                continue;
              const Vec3 unit = lab / lab_length;
              const double entry_fast = mm_fast - depth[d] * unit.dot(p.fast);
              const double entry_slow = mm_slow - depth[d] * unit.dot(p.slow);

              const Vec3 entry = p.lab_coord_mm(entry_fast, entry_slow);
              const double entry_length = entry.norm();
              if (!(entry_length > 0.0))
                continue;
              const Vec3 s_prime = entry * (length / entry_length);
              const Vec3 difference = s_prime - s1;
              const double eps1 =
                  Scan::degrees(frame.e1.dot(difference) / length);
              const double eps2 =
                  Scan::degrees(frame.e2.dot(difference) / length);

              const double u1 = eps1 / options.sigma_d;
              const double u2 = eps2 / options.sigma_d;
              out[box.at(x, y, z)] += share * depth_weight[d] * along_rotation *
                                      std::exp(-0.5 * (u1 * u1 + u2 * u2));
            }
          }
        }
      }
    }
  }

  double total = 0.0;
  for (std::size_t i = 0; i < out.size(); ++i) {
    if ((box.mask[i] & shoebox_mask::kValid) == 0)
      out[i] = 0.0;
    total += out[i];
  }
  if (total > 0.0) {
    for (double &v : out)
      v /= total;
  }
  return out;
}

Moments moments_of(const Shoebox &box, const std::vector<double> &counts,
                   bool subtract_background) {
  Moments out;
  if (counts.size() != box.size())
    return out;
  double sum = 0.0, mf = 0.0, ms = 0.0, mz = 0.0;
  double sf = 0.0, ss = 0.0, sz = 0.0;
  for (int pass = 0; pass < 2; ++pass) {
    for (std::int32_t z = 0; z < box.nz(); ++z) {
      for (std::int32_t y = 0; y < box.ny(); ++y) {
        for (std::int32_t x = 0; x < box.nx(); ++x) {
          const std::size_t at = box.at(x, y, z);
          if ((box.mask[at] & shoebox_mask::kValid) == 0)
            continue;
          double value = counts[at];
          if (subtract_background)
            value -= static_cast<double>(box.background[at]);
          if (!(value > 0.0))
            continue;
          // In pixels and images, counting from the corner of the box: the
          // units the detector records in, so neither side is transformed.
          const double fx = static_cast<double>(x) + 0.5;
          const double fy = static_cast<double>(y) + 0.5;
          const double fz = static_cast<double>(z) + 0.5;
          if (pass == 0) {
            sum += value;
            mf += value * fx;
            ms += value * fy;
            mz += value * fz;
          } else {
            sf += value * (fx - mf) * (fx - mf);
            ss += value * (fy - ms) * (fy - ms);
            sz += value * (fz - mz) * (fz - mz);
          }
        }
      }
    }
    if (pass == 0) {
      if (!(sum > 0.0))
        return out;
      mf /= sum;
      ms /= sum;
      mz /= sum;
    }
  }
  out.valid = true;
  out.total = sum;
  out.com_fast = mf;
  out.com_slow = ms;
  out.com_z = mz;
  out.width_fast = std::sqrt(sf / sum);
  out.width_slow = std::sqrt(ss / sum);
  out.width_z = std::sqrt(sz / sum);
  return out;
}

} // namespace mxi
