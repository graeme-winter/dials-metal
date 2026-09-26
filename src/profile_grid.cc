#include "profile_grid.hh"

#include <cmath>

namespace mxi {

void ProfileGrid::reset() {
  value.assign(size(), 0.0);
  n_spots = 0;
  counts_added = 0.0;
  counts_outside = 0.0;
}

void ProfileGrid::normalise() {
  double total = 0.0;
  for (double v : value)
    total += v;
  if (!(total > 0.0))
    return;
  for (double &v : value)
    v /= total;
}

ProfileGrid make_grid(int n, double sigma_d, double sigma_m,
                      double half_width) {
  ProfileGrid grid;
  grid.n = n;
  grid.sigma_d = sigma_d;
  grid.sigma_m = sigma_m;
  grid.half_width = half_width;
  grid.reset();
  return grid;
}

void add_to_grid(const Experiment &e, const Shoebox &box, const Vec3 &s1,
                 double phi_calculated, ProfileGrid *grid, int subdivisions,
                 bool recentre) {
  if (grid == nullptr || subdivisions < 1)
    return;
  if (box.panel < 0 ||
      static_cast<std::size_t>(box.panel) >= e.detector.size()) {
    return;
  }
  const Panel &p = e.detector[static_cast<std::size_t>(box.panel)];
  const KabschFrame frame = kabsch_frame(e, s1);
  if (!frame.valid)
    return;

  const int side = grid->side();
  const double span_d = grid->half_width * grid->sigma_d;
  const double span_m = grid->half_width * grid->sigma_m;
  if (!(span_d > 0.0) || !(span_m > 0.0))
    return;
  // Grid point i covers [-span + i * step, -span + (i + 1) * step).
  const double step_d = 2.0 * span_d / static_cast<double>(side);
  const double step_m = 2.0 * span_m / static_cast<double>(side);
  const double width = Scan::radians(e.scan.osc_width);
  const double share = 1.0 / static_cast<double>(subdivisions * subdivisions);

  // Where this spot's own centre of mass sits in the frame, if it is to be
  // shifted onto it.
  double centre1 = 0.0, centre2 = 0.0, centre3 = 0.0;
  if (recentre) {
    const SpotMoments moments = spot_moments(e, box, s1, phi_calculated);
    if (!moments.valid)
      return;
    centre1 = moments.centre1;
    centre2 = moments.centre2;
    centre3 = moments.centre3;
  }

  ++grid->n_spots;
  for (std::int32_t z = 0; z < box.nz(); ++z) {
    // The angular range this image covers, as an eps3 range. eps3 runs
    // backwards in phi when zeta is negative, so the ends are sorted rather
    // than assumed ordered -- a reflection on the other side of the rotation
    // axis would otherwise contribute nothing at all.
    const double image = static_cast<double>(box.bbox[4] + z);
    const double phi_low =
        Scan::radians(e.scan.osc_start) +
        (image - static_cast<double>(e.scan.z_offset)) * width;
    const double phi_high = phi_low + width;
    double e3_low =
        Scan::degrees(frame.zeta * (phi_low - phi_calculated)) - centre3;
    double e3_high =
        Scan::degrees(frame.zeta * (phi_high - phi_calculated)) - centre3;
    if (e3_low > e3_high)
      std::swap(e3_low, e3_high);
    const double e3_span = e3_high - e3_low;

    for (std::int32_t y = 0; y < box.ny(); ++y) {
      for (std::int32_t x = 0; x < box.nx(); ++x) {
        const std::size_t at = box.at(x, y, z);
        if ((box.mask[at] & shoebox_mask::kValid) == 0)
          continue;
        const double count = static_cast<double>(box.data[at]) -
                             static_cast<double>(box.background[at]);
        if (!(count > 0.0))
          continue;

        for (int sy = 0; sy < subdivisions; ++sy) {
          for (int sx = 0; sx < subdivisions; ++sx) {
            // The centre of this subdivision, in pixels.
            const double px = static_cast<double>(box.bbox[0] + x) +
                              (static_cast<double>(sx) + 0.5) /
                                  static_cast<double>(subdivisions);
            const double py = static_cast<double>(box.bbox[2] + y) +
                              (static_cast<double>(sy) + 0.5) /
                                  static_cast<double>(subdivisions);
            Epsilon eps =
                epsilon_of(e, frame, p, px, py, phi_low, phi_calculated);
            eps.e1 -= centre1;
            eps.e2 -= centre2;
            const int i1 =
                static_cast<int>(std::floor((eps.e1 + span_d) / step_d));
            const int i2 =
                static_cast<int>(std::floor((eps.e2 + span_d) / step_d));
            const double portion = count * share;
            grid->counts_added += portion;
            if (i1 < 0 || i1 >= side || i2 < 0 || i2 >= side) {
              grid->counts_outside += portion;
              continue;
            }

            // Share along e3 by how much of the image's angular range falls in
            // each plane. Geometric, with no model in it.
            bool landed = false;
            for (int i3 = 0; i3 < side; ++i3) {
              const double low = -span_m + static_cast<double>(i3) * step_m;
              const double high = low + step_m;
              const double overlap =
                  std::fmin(e3_high, high) - std::fmax(e3_low, low);
              if (!(overlap > 0.0))
                continue;
              const double fraction = e3_span > 0.0 ? overlap / e3_span : 1.0;
              grid->value[grid->at(i1, i2, i3)] += portion * fraction;
              landed = true;
            }
            if (!landed)
              grid->counts_outside += portion;
          }
        }
      }
    }
  }
}

SpotMoments spot_moments(const Experiment &e, const Shoebox &box,
                         const Vec3 &s1, double phi_calculated) {
  SpotMoments out;
  if (box.panel < 0 ||
      static_cast<std::size_t>(box.panel) >= e.detector.size()) {
    return out;
  }
  const Panel &p = e.detector[static_cast<std::size_t>(box.panel)];
  const KabschFrame frame = kabsch_frame(e, s1);
  if (!frame.valid)
    return out;
  const double width = Scan::radians(e.scan.osc_width);

  // Two passes: the centroid, then the spread about it.
  // Two axes fixed in the laboratory, perpendicular to the beam: the rotation
  // axis, and the direction perpendicular to both it and the beam.
  const Vec3 s0 = e.beam.s0();
  const Vec3 axis = e.goniometer.lab_axis();
  Vec3 along = axis - s0 * (axis.dot(s0) / s0.dot(s0));
  const double along_length = along.norm();
  if (!(along_length > 0.0))
    return out;
  along = along / along_length;
  const Vec3 across = s0.cross(along).normalized();

  double sum = 0.0, m1 = 0.0, m2 = 0.0, m3 = 0.0;
  double s11 = 0.0, s22 = 0.0, s33 = 0.0;
  double ma = 0.0, mc = 0.0, saa = 0.0, scc = 0.0;
  for (int pass = 0; pass < 2; ++pass) {
    for (std::int32_t z = 0; z < box.nz(); ++z) {
      const double image = static_cast<double>(box.bbox[4] + z) + 0.5;
      const double phi = Scan::radians(e.scan.osc_start) +
                         (image - static_cast<double>(e.scan.z_offset)) * width;
      for (std::int32_t y = 0; y < box.ny(); ++y) {
        for (std::int32_t x = 0; x < box.nx(); ++x) {
          const std::size_t at = box.at(x, y, z);
          if ((box.mask[at] & shoebox_mask::kValid) == 0)
            continue;
          const double count = static_cast<double>(box.data[at]) -
                               static_cast<double>(box.background[at]);
          if (!(count > 0.0))
            continue;
          const Epsilon eps = epsilon_of(
              e, frame, p, static_cast<double>(box.bbox[0] + x) + 0.5,
              static_cast<double>(box.bbox[2] + y) + 0.5, phi, phi_calculated);
          // The same offset, resolved in the laboratory instead.
          const auto mm =
              p.px_to_mm(static_cast<double>(box.bbox[0] + x) + 0.5,
                         static_cast<double>(box.bbox[2] + y) + 0.5);
          const Vec3 lab = p.lab_coord_mm(mm.first, mm.second);
          const double length = frame.s1.norm();
          const double lab_length = lab.norm();
          const Vec3 offset = lab_length > 0.0
                                  ? lab * (length / lab_length) - frame.s1
                                  : Vec3{0.0, 0.0, 0.0};
          const double ea = Scan::degrees(along.dot(offset) / length);
          const double ec = Scan::degrees(across.dot(offset) / length);
          if (pass == 0) {
            sum += count;
            m1 += count * eps.e1;
            m2 += count * eps.e2;
            m3 += count * eps.e3;
            ma += count * ea;
            mc += count * ec;
          } else {
            s11 += count * (eps.e1 - m1) * (eps.e1 - m1);
            s22 += count * (eps.e2 - m2) * (eps.e2 - m2);
            s33 += count * (eps.e3 - m3) * (eps.e3 - m3);
            saa += count * (ea - ma) * (ea - ma);
            scc += count * (ec - mc) * (ec - mc);
          }
        }
      }
    }
    if (pass == 0) {
      if (!(sum > 1.0))
        return out;
      m1 /= sum;
      m2 /= sum;
      m3 /= sum;
      ma /= sum;
      mc /= sum;
    }
  }
  out.valid = true;
  out.counts = sum;
  out.centre1 = m1;
  out.centre2 = m2;
  out.centre3 = m3;
  out.width1 = std::sqrt(s11 / sum);
  out.width2 = std::sqrt(s22 / sum);
  out.width3 = std::sqrt(s33 / sum);
  out.width_along_axis = std::sqrt(saa / sum);
  out.width_across_axis = std::sqrt(scc / sum);

  // Where the spot sits, from its predicted direction rather than its pixels.
  const Vec3 normal = p.fast.cross(p.slow).normalized();
  const Vec3 direction = s1 / s1.norm();
  out.obliquity = std::acos(std::fmin(1.0, std::fabs(direction.dot(normal))));
  const double to_panel = p.origin.dot(normal) / direction.dot(normal);
  const Vec3 hit = direction * to_panel;
  out.path_mm = hit.norm();
  // Distance from where the beam itself strikes.
  const Vec3 beam = e.beam.s0() * (-1.0 / e.beam.s0().norm());
  const double beam_along = p.origin.dot(normal) / beam.dot(normal);
  out.radius_mm = (hit - beam * beam_along).norm();
  return out;
}

double sensor_depth_width(double mu, double thickness, double obliquity,
                          double path) {
  if (!(mu > 0.0) || !(thickness > 0.0) || !(path > 0.0))
    return 0.0;
  // Depth of absorption: exponential, cut off at the back of the sensor.
  const double transmitted = std::exp(-mu * thickness);
  const double survive = 1.0 - transmitted;
  if (!(survive > 0.0))
    return 0.0;
  const double mean = 1.0 / mu - thickness * transmitted / survive;
  const double second =
      (2.0 / (mu * mu) * survive -
       (thickness * thickness + 2.0 * thickness / mu) * transmitted) /
      survive;
  const double variance = std::fmax(0.0, second - mean * mean);
  const double sigma_depth = std::sqrt(variance);
  // Projected onto the face it is tan(obliquity) * sigma_depth; seen from the
  // crystal that subtends sin(obliquity) * sigma_depth / path.
  return Scan::degrees(std::sin(obliquity) * sigma_depth / path);
}

double source_extent_width(double extent, double obliquity, double path) {
  if (!(extent > 0.0) || !(path > 0.0))
    return 0.0;
  // extent * tan(obliquity) on the face; seen from the crystal that subtends
  // extent * sin(obliquity) / path, the same form as the sensor's smear.
  return Scan::degrees(std::sin(obliquity) * extent / path);
}

} // namespace mxi
