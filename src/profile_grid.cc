#include "profile_grid.h"

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
  for (double v : value) total += v;
  if (!(total > 0.0)) return;
  for (double &v : value) v /= total;
}

ProfileGrid make_grid(int n, double sigma_d, double sigma_m, double half_width) {
  ProfileGrid grid;
  grid.n = n;
  grid.sigma_d = sigma_d;
  grid.sigma_m = sigma_m;
  grid.half_width = half_width;
  grid.reset();
  return grid;
}

void add_to_grid(const Experiment &e, const Shoebox &box, const Vec3 &s1,
                 double phi_calculated, ProfileGrid *grid, int subdivisions) {
  if (grid == nullptr || subdivisions < 1) return;
  if (box.panel < 0 || static_cast<std::size_t>(box.panel) >= e.detector.size()) {
    return;
  }
  const Panel &p = e.detector[static_cast<std::size_t>(box.panel)];
  const KabschFrame frame = kabsch_frame(e, s1);
  if (!frame.valid) return;

  const int side = grid->side();
  const double span_d = grid->half_width * grid->sigma_d;
  const double span_m = grid->half_width * grid->sigma_m;
  if (!(span_d > 0.0) || !(span_m > 0.0)) return;
  // Grid point i covers [-span + i * step, -span + (i + 1) * step).
  const double step_d = 2.0 * span_d / static_cast<double>(side);
  const double step_m = 2.0 * span_m / static_cast<double>(side);
  const double width = Scan::radians(e.scan.osc_width);
  const double share = 1.0 / static_cast<double>(subdivisions * subdivisions);

  ++grid->n_spots;
  for (std::int32_t z = 0; z < box.nz(); ++z) {
    // The angular range this image covers, as an eps3 range. eps3 runs
    // backwards in phi when zeta is negative, so the ends are sorted rather
    // than assumed ordered -- a reflection on the other side of the rotation
    // axis would otherwise contribute nothing at all.
    const double image = static_cast<double>(box.bbox[4] + z);
    const double phi_low = Scan::radians(e.scan.osc_start) +
                           (image - static_cast<double>(e.scan.z_offset)) * width;
    const double phi_high = phi_low + width;
    double e3_low = Scan::degrees(frame.zeta * (phi_low - phi_calculated));
    double e3_high = Scan::degrees(frame.zeta * (phi_high - phi_calculated));
    if (e3_low > e3_high) std::swap(e3_low, e3_high);
    const double e3_span = e3_high - e3_low;

    for (std::int32_t y = 0; y < box.ny(); ++y) {
      for (std::int32_t x = 0; x < box.nx(); ++x) {
        const std::size_t at = box.at(x, y, z);
        if (box.mask[at] == 0) continue;
        const double count = static_cast<double>(box.data[at]) -
                             static_cast<double>(box.background[at]);
        if (!(count > 0.0)) continue;

        for (int sy = 0; sy < subdivisions; ++sy) {
          for (int sx = 0; sx < subdivisions; ++sx) {
            // The centre of this subdivision, in pixels.
            const double px = static_cast<double>(box.bbox[0] + x) +
                              (static_cast<double>(sx) + 0.5) /
                                  static_cast<double>(subdivisions);
            const double py = static_cast<double>(box.bbox[2] + y) +
                              (static_cast<double>(sy) + 0.5) /
                                  static_cast<double>(subdivisions);
            const Epsilon eps =
                epsilon_of(e, frame, p, px, py, phi_low, phi_calculated);
            const int i1 = static_cast<int>(
                std::floor((eps.e1 + span_d) / step_d));
            const int i2 = static_cast<int>(
                std::floor((eps.e2 + span_d) / step_d));
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
              if (!(overlap > 0.0)) continue;
              const double fraction =
                  e3_span > 0.0 ? overlap / e3_span : 1.0;
              grid->value[grid->at(i1, i2, i3)] += portion * fraction;
              landed = true;
            }
            if (!landed) grid->counts_outside += portion;
          }
        }
      }
    }
  }
}

}  // namespace mxi
