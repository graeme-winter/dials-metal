#include "reference.h"

#include <algorithm>
#include <cmath>

#include "profile_model.h"

namespace mxi {

Transformed transform_shoebox(const Experiment &e, const Shoebox &box,
                              const Vec3 &s1, double phi_calculated,
                              const GridSpec &spec) {
  Transformed out;
  if (box.panel < 0 || static_cast<std::size_t>(box.panel) >= e.detector.size()) {
    return out;
  }
  if (spec.subdivisions < 1) return out;
  const Panel &p = e.detector[static_cast<std::size_t>(box.panel)];
  const KabschFrame frame = kabsch_frame(e, s1);
  if (!frame.valid) return out;

  const int side = spec.side();
  const double span_d = spec.half_width * spec.sigma_d;
  const double span_m = spec.half_width * spec.sigma_m;
  if (!(span_d > 0.0) || !(span_m > 0.0)) return out;
  const double step_d = 2.0 * span_d / static_cast<double>(side);
  const double step_m = 2.0 * span_m / static_cast<double>(side);
  const double osc = Scan::radians(e.scan.osc_width);
  const double share =
      1.0 / static_cast<double>(spec.subdivisions * spec.subdivisions);

  out.data.assign(spec.size(), 0.0);
  out.background.assign(spec.size(), 0.0);
  out.coverage.assign(spec.size(), 0.0);
  double inside = 0.0;
  double total = 0.0;

  // eps1 and eps2 depend only on where a subdivision sits on the detector, so
  // the subdivided face is computed once rather than once per image. With five
  // subdivisions and eight images that is 88000 evaluations of a mapping with
  // an exp() in it per reflection, against 11000 -- and it was 215 seconds of
  // a 234 second run before this.
  const std::int32_t fine_x = box.nx() * spec.subdivisions;
  const std::int32_t fine_y = box.ny() * spec.subdivisions;
  std::vector<int> face_i1(static_cast<std::size_t>(fine_x) *
                           static_cast<std::size_t>(fine_y));
  std::vector<int> face_i2(face_i1.size());
  {
    // phi does not matter for eps1 and eps2, so any value in the box will do;
    // the reflection's own is the natural one.
    for (std::int32_t fy = 0; fy < fine_y; ++fy) {
      const double py = static_cast<double>(box.bbox[2]) +
                        (static_cast<double>(fy) + 0.5) /
                            static_cast<double>(spec.subdivisions);
      for (std::int32_t fx = 0; fx < fine_x; ++fx) {
        const double px = static_cast<double>(box.bbox[0]) +
                          (static_cast<double>(fx) + 0.5) /
                              static_cast<double>(spec.subdivisions);
        const Epsilon eps =
            epsilon_of(e, frame, p, px, py, phi_calculated, phi_calculated);
        const std::size_t at = static_cast<std::size_t>(fy) *
                                   static_cast<std::size_t>(fine_x) +
                               static_cast<std::size_t>(fx);
        face_i1[at] = static_cast<int>(std::floor((eps.e1 + span_d) / step_d));
        face_i2[at] = static_cast<int>(std::floor((eps.e2 + span_d) / step_d));
      }
    }
  }

  for (std::int32_t z = 0; z < box.nz(); ++z) {
    const double image = static_cast<double>(box.bbox[4] + z);
    const double phi_low = Scan::radians(e.scan.osc_start) +
                           (image - static_cast<double>(e.scan.z_offset)) * osc;
    const double phi_high = phi_low + osc;
    // eps3 runs backwards in phi where zeta is negative, so the ends are
    // sorted rather than assumed ordered: a reflection on the other side of
    // the rotation axis would otherwise contribute nothing at all.
    double e3_low = Scan::degrees(frame.zeta * (phi_low - phi_calculated));
    double e3_high = Scan::degrees(frame.zeta * (phi_high - phi_calculated));
    if (e3_low > e3_high) std::swap(e3_low, e3_high);
    const double e3_span = e3_high - e3_low;
    if (!(e3_span > 0.0)) continue;

    for (std::int32_t y = 0; y < box.ny(); ++y) {
      for (std::int32_t x = 0; x < box.nx(); ++x) {
        const std::size_t at = box.at(x, y, z);
        // A masked voxel is no measurement at all, not a zero one: it must not
        // reach the grid as a count and must not claim coverage there either,
        // or a module gap becomes a hole in the profile.
        if (box.mask[at] == 0) continue;
        const double count = static_cast<double>(box.data[at]);
        const double back = static_cast<double>(box.background[at]);
        // What counts as "inside" is the SIGNAL in the foreground, not every
        // count in the box. The box is 1.9 times wider than the foreground on
        // the detector so that it has a background rim, and the grid spans the
        // foreground -- so the rim lies outside the grid by construction, and
        // measuring against every count made every reflection look as though
        // it mostly missed. That rejected all but 8 of 20472 as unfit to learn
        // from.
        if (box.mask[at] & shoebox_mask::kForeground) {
          total += std::max(count - back, 0.0);
        }

        for (int sy = 0; sy < spec.subdivisions; ++sy) {
          for (int sx = 0; sx < spec.subdivisions; ++sx) {
            const std::size_t on_face =
                static_cast<std::size_t>(y * spec.subdivisions + sy) *
                    static_cast<std::size_t>(fine_x) +
                static_cast<std::size_t>(x * spec.subdivisions + sx);
            const int i1 = face_i1[on_face];
            const int i2 = face_i2[on_face];
            if (i1 < 0 || i1 >= side || i2 < 0 || i2 >= side) continue;

            // One image covers a RANGE of eps3, which is shared between the
            // grid planes it overlaps in proportion to the overlap. A whole
            // image dropped into one plane is what makes a rocking curve look
            // like a step.
            const int j_low =
                static_cast<int>(std::floor((e3_low + span_m) / step_m));
            const int j_high =
                static_cast<int>(std::floor((e3_high + span_m) / step_m));
            for (int j = std::max(j_low, 0); j <= std::min(j_high, side - 1);
                 ++j) {
              const double plane_low = -span_m + static_cast<double>(j) * step_m;
              const double plane_high = plane_low + step_m;
              const double overlap = std::min(e3_high, plane_high) -
                                     std::max(e3_low, plane_low);
              if (!(overlap > 0.0)) continue;
              const double weight = share * (overlap / e3_span);
              const std::size_t into = spec.at(i1, i2, j);
              out.data[into] += weight * count;
              out.background[into] += weight * back;
              out.coverage[into] += weight;
              if (box.mask[at] & shoebox_mask::kForeground) {
                inside += weight * std::max(count - back, 0.0);
              }
            }
          }
        }
      }
    }
  }
  out.outside = total > 0.0 ? 1.0 - inside / total : 1.0;
  out.valid = true;
  return out;
}

std::size_t ReferenceProfiles::region_of(const Panel &panel,
                                         std::size_t which_panel,
                                         double px_fast, double px_slow) const {
  const double fast_size = static_cast<double>(std::max<std::int64_t>(panel.image_size[0], 1));
  const double slow_size = static_cast<double>(std::max<std::int64_t>(panel.image_size[1], 1));
  const double d = static_cast<double>(std::max(divisions, 1));
  int i = static_cast<int>(px_fast / fast_size * d);
  int j = static_cast<int>(px_slow / slow_size * d);
  i = std::clamp(i, 0, divisions - 1);
  j = std::clamp(j, 0, divisions - 1);
  const std::size_t within =
      static_cast<std::size_t>(j) * static_cast<std::size_t>(divisions) +
      static_cast<std::size_t>(i);
  return std::min(which_panel, panels - 1) * regions_per_panel() + within;
}

ReferenceProfiles make_reference(const GridSpec &spec, int divisions,
                                 std::size_t panels) {
  ReferenceProfiles out;
  out.spec = spec;
  out.divisions = std::max(divisions, 1);
  out.panels = std::max<std::size_t>(panels, 1);
  out.profile.assign(out.region_count(), std::vector<double>(spec.size(), 0.0));
  out.spots.assign(out.region_count(), 0);
  return out;
}

bool add_reference(ReferenceProfiles *reference, std::size_t region,
                   const Transformed &t) {
  if (reference == nullptr || !t.valid) return false;
  if (region >= reference->profile.size()) return false;
  if (t.data.size() != reference->spec.size()) return false;

  // Normalised to unit sum, so the profile is the average shape rather than
  // the average spot: without this a single strong reflection outvotes a
  // hundred weak ones and the profile is that reflection's.
  double sum = 0.0;
  for (std::size_t i = 0; i < t.data.size(); ++i) {
    sum += t.data[i] - t.background[i];
  }
  if (!(sum > 0.0)) return false;
  std::vector<double> &into = reference->profile[region];
  for (std::size_t i = 0; i < t.data.size(); ++i) {
    into[i] += (t.data[i] - t.background[i]) / sum;
  }
  ++reference->spots[region];
  return true;
}

void finalise_reference(ReferenceProfiles *reference, std::size_t least) {
  if (reference == nullptr) return;
  const std::size_t n = reference->spec.size();

  // The whole-detector average, for regions that saw too few spots to have a
  // profile of their own. A region at the corner of a detector may have almost
  // nothing in it, and an empty profile fits nothing at all.
  std::vector<double> everything(n, 0.0);
  std::size_t everything_spots = 0;
  for (std::size_t r = 0; r < reference->profile.size(); ++r) {
    if (reference->spots[r] == 0) continue;
    for (std::size_t i = 0; i < n; ++i) everything[i] += reference->profile[r][i];
    everything_spots += reference->spots[r];
  }
  double whole = 0.0;
  for (double v : everything) whole += v;
  if (whole > 0.0) {
    for (double &v : everything) v /= whole;
  }

  for (std::size_t r = 0; r < reference->profile.size(); ++r) {
    if (reference->spots[r] < least) {
      reference->profile[r] = everything;
      continue;
    }
    double sum = 0.0;
    for (double v : reference->profile[r]) sum += v;
    if (!(sum > 0.0)) {
      reference->profile[r] = everything;
      continue;
    }
    for (double &v : reference->profile[r]) v /= sum;
  }
  (void)everything_spots;
  reference->finalised = true;
}

ProfileFit fit_profile(const std::vector<double> &reference, const Transformed &t,
                       double gain, int iterations) {
  ProfileFit out;
  if (!t.valid || reference.size() != t.data.size()) return out;

  double profile_sum = 0.0;
  for (double v : reference) profile_sum += v;
  if (!(profile_sum > 0.0)) return out;

  // A first guess from the plain sum, to weight the first round with.
  double scale = 0.0;
  for (std::size_t i = 0; i < t.data.size(); ++i) {
    if (t.coverage[i] > 0.0) scale += t.data[i] - t.background[i];
  }
  scale /= profile_sum;

  double variance = 0.0;
  for (int round = 0; round < std::max(iterations, 1); ++round) {
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t i = 0; i < t.data.size(); ++i) {
      // A grid point no pixel reached is not a measurement and cannot be
      // fitted to. Counting it as an observed zero pulls every intensity down.
      if (!(t.coverage[i] > 0.0)) continue;
      const double model = t.background[i] + scale * reference[i];
      // The Poisson variance of what the model says should be there, floored
      // so that a point the model puts at zero still has a weight. Using the
      // OBSERVED count instead biases the fit: a point that happened to record
      // nothing would be given infinite weight.
      const double v = std::max(gain * model, 1e-6);
      numerator += reference[i] * (t.data[i] - t.background[i]) / v;
      denominator += reference[i] * reference[i] / v;
    }
    if (!(denominator > 0.0)) return out;
    scale = numerator / denominator;
    variance = 1.0 / denominator;
    out.iterations = round + 1;
  }

  out.intensity = scale * profile_sum;
  out.variance = variance * profile_sum * profile_sum;

  // How well the profile actually describes this reflection.
  double n = 0.0, sp = 0.0, sd = 0.0, spp = 0.0, sdd = 0.0, spd = 0.0;
  for (std::size_t i = 0; i < t.data.size(); ++i) {
    if (!(t.coverage[i] > 0.0)) continue;
    const double a = reference[i];
    const double b = t.data[i] - t.background[i];
    n += 1.0;
    sp += a;
    sd += b;
    spp += a * a;
    sdd += b * b;
    spd += a * b;
  }
  if (n > 1.0) {
    const double cov = spd / n - (sp / n) * (sd / n);
    const double va = spp / n - (sp / n) * (sp / n);
    const double vb = sdd / n - (sd / n) * (sd / n);
    if (va > 0.0 && vb > 0.0) out.correlation = cov / std::sqrt(va * vb);
  }
  out.valid = true;
  return out;
}

}  // namespace mxi
