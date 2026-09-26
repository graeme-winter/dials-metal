#include "reference.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

#include "profile_model.hh"

namespace mxi {

namespace {
//: The least variance a measurement can have, in counts. A point with no
//: signal still carries the background's noise; anything smaller is a weight
//: this fit has no business trusting.
constexpr double kLeastVariance = 0.05;
} // namespace

namespace {

//: Which grid cells each pixel of a shoebox reaches, and how many of its
//: subdivisions reach each: for pixel k, entries [start[k], start[k + 1]).
//:
//: The geometry both transforms share, in one place so that the two cannot
//: drift. eps1 and eps2 are computed at the pixel CORNERS and interpolated
//: within each pixel -- (nx + 1)(ny + 1) calls to epsilon_of rather than 25 nx
//: ny -- and a subdivision within a thousandth of a cell of a boundary is done
//: exactly. Interpolation alone is good to about 1e-7 degrees and put a
//: subdivision across a boundary in 147 of 3000 real boxes; with the guard the
//: cells are the direct computation's.
struct PixelCells {
  std::vector<std::uint32_t> start;
  std::vector<std::uint32_t> cell; // i2 * side + i1
  std::vector<std::uint16_t> hits;
};

PixelCells pixel_cells(const Experiment &e, const KabschFrame &frame,
                       const Panel &p, const Shoebox &box, int sub, int side,
                       double span_d, double step_d, double phi_calculated) {
  const std::int32_t nx = box.nx(), ny = box.ny();
  const std::size_t cw = static_cast<std::size_t>(nx) + 1;
  std::vector<double> c1(cw * (static_cast<std::size_t>(ny) + 1));
  std::vector<double> c2(c1.size());
  for (std::int32_t cy = 0; cy <= ny; ++cy) {
    const double py = static_cast<double>(box.bbox[2] + cy);
    for (std::int32_t cx = 0; cx <= nx; ++cx) {
      const double px = static_cast<double>(box.bbox[0] + cx);
      const Epsilon eps =
          epsilon_of(e, frame, p, px, py, phi_calculated, phi_calculated);
      const std::size_t at =
          static_cast<std::size_t>(cy) * cw + static_cast<std::size_t>(cx);
      c1[at] = eps.e1;
      c2[at] = eps.e2;
    }
  }
  const std::size_t pixels =
      static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
  PixelCells out;
  out.start.assign(pixels + 1, 0);
  out.cell.reserve(pixels * 4);
  out.hits.reserve(pixels * 4);
  const auto near_boundary = [](double f) {
    const double r = f - std::floor(f);
    return std::min(r, 1.0 - r) < 1e-3;
  };
  for (std::int32_t y = 0; y < ny; ++y) {
    for (std::int32_t x = 0; x < nx; ++x) {
      const std::size_t k =
          static_cast<std::size_t>(y) * static_cast<std::size_t>(nx) +
          static_cast<std::size_t>(x);
      out.start[k] = static_cast<std::uint32_t>(out.cell.size());
      const std::size_t a00 =
          static_cast<std::size_t>(y) * cw + static_cast<std::size_t>(x);
      const std::size_t a10 = a00 + 1, a01 = a00 + cw, a11 = a01 + 1;
      for (int sy = 0; sy < sub; ++sy) {
        const double v =
            (static_cast<double>(sy) + 0.5) / static_cast<double>(sub);
        for (int sx = 0; sx < sub; ++sx) {
          const double u =
              (static_cast<double>(sx) + 0.5) / static_cast<double>(sub);
          double e1 = (1.0 - u) * (1.0 - v) * c1[a00] +
                      u * (1.0 - v) * c1[a10] + (1.0 - u) * v * c1[a01] +
                      u * v * c1[a11];
          double e2 = (1.0 - u) * (1.0 - v) * c2[a00] +
                      u * (1.0 - v) * c2[a10] + (1.0 - u) * v * c2[a01] +
                      u * v * c2[a11];
          if (near_boundary((e1 + span_d) / step_d) ||
              near_boundary((e2 + span_d) / step_d)) {
            const Epsilon eps = epsilon_of(
                e, frame, p, static_cast<double>(box.bbox[0] + x) + u,
                static_cast<double>(box.bbox[2] + y) + v, phi_calculated,
                phi_calculated);
            e1 = eps.e1;
            e2 = eps.e2;
          }
          const int i1 = static_cast<int>(std::floor((e1 + span_d) / step_d));
          const int i2 = static_cast<int>(std::floor((e2 + span_d) / step_d));
          if (i1 < 0 || i1 >= side || i2 < 0 || i2 >= side)
            continue;
          const std::uint32_t id = static_cast<std::uint32_t>(i2 * side + i1);
          bool found = false;
          for (std::size_t q = out.start[k]; q < out.cell.size(); ++q) {
            if (out.cell[q] == id) {
              ++out.hits[q];
              found = true;
              break;
            }
          }
          if (!found) {
            out.cell.push_back(id);
            out.hits.push_back(1);
          }
        }
      }
    }
  }
  out.start[pixels] = static_cast<std::uint32_t>(out.cell.size());
  return out;
}

//: The eps3 planes one image of the box overlaps, each with the fraction of
//: the image's eps3 range that falls in it. `spans` is false when the image
//: covers no eps3 range at all, which the direct transforms skip entirely --
//: as distinct from an image that spans a range but meets no plane of the grid,
//: which still counts toward how much of the reflection lies outside it.
struct Planes {
  bool spans = false;
  std::vector<std::pair<int, double>> planes; // j, overlap / e3 range
};

Planes plane_weights(const Experiment &e, const KabschFrame &frame,
                     const Shoebox &box, std::int32_t z, int side,
                     double span_m, double step_m, double phi_calculated) {
  Planes out;
  const double osc = Scan::radians(e.scan.osc_width);
  const double image = static_cast<double>(box.bbox[4] + z);
  const double phi_low = Scan::radians(e.scan.osc_start) +
                         (image - static_cast<double>(e.scan.z_offset)) * osc;
  const double phi_high = phi_low + osc;
  double e3_low = Scan::degrees(frame.zeta * (phi_low - phi_calculated));
  double e3_high = Scan::degrees(frame.zeta * (phi_high - phi_calculated));
  if (e3_low > e3_high)
    std::swap(e3_low, e3_high);
  const double e3_span = e3_high - e3_low;
  if (!(e3_span > 0.0))
    return out;
  out.spans = true;
  const int j_low = static_cast<int>(std::floor((e3_low + span_m) / step_m));
  const int j_high = static_cast<int>(std::floor((e3_high + span_m) / step_m));
  for (int j = std::max(j_low, 0); j <= std::min(j_high, side - 1); ++j) {
    const double plane_low = -span_m + static_cast<double>(j) * step_m;
    const double plane_high = plane_low + step_m;
    const double overlap =
        std::min(e3_high, plane_high) - std::max(e3_low, plane_low);
    if (overlap > 0.0)
      out.planes.emplace_back(j, overlap / e3_span);
  }
  return out;
}

} // namespace

Transformed transform_shoebox_direct(const Experiment &e, const Shoebox &box,
                                     const Vec3 &s1, double phi_calculated,
                                     const GridSpec &spec) {
  Transformed out;
  if (box.panel < 0 ||
      static_cast<std::size_t>(box.panel) >= e.detector.size()) {
    return out;
  }
  if (spec.subdivisions < 1)
    return out;
  const Panel &p = e.detector[static_cast<std::size_t>(box.panel)];
  const KabschFrame frame = kabsch_frame(e, s1);
  if (!frame.valid)
    return out;

  const int side = spec.side();
  const double span_d = spec.half_width * spec.sigma_d;
  const double span_m = spec.half_width * spec.sigma_m;
  if (!(span_d > 0.0) || !(span_m > 0.0))
    return out;
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
        const std::size_t at =
            static_cast<std::size_t>(fy) * static_cast<std::size_t>(fine_x) +
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
    if (e3_low > e3_high)
      std::swap(e3_low, e3_high);
    const double e3_span = e3_high - e3_low;
    if (!(e3_span > 0.0))
      continue;

    for (std::int32_t y = 0; y < box.ny(); ++y) {
      for (std::int32_t x = 0; x < box.nx(); ++x) {
        const std::size_t at = box.at(x, y, z);
        // A masked voxel is no measurement at all, not a zero one: it must not
        // reach the grid as a count and must not claim coverage there either,
        // or a module gap becomes a hole in the profile.
        // Validity, not the whole mask. A bad pixel keeps its region flag so
        // that the fit knows part of the reflection is missing, so testing the
        // mask against zero stopped excluding them -- and they went onto the
        // grid as genuine zeroes, which is a hole in every profile learned
        // from a reflection that crosses a module gap.
        if ((box.mask[at] & shoebox_mask::kValid) == 0)
          continue;
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
            if (i1 < 0 || i1 >= side || i2 < 0 || i2 >= side)
              continue;

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
              const double plane_low =
                  -span_m + static_cast<double>(j) * step_m;
              const double plane_high = plane_low + step_m;
              const double overlap =
                  std::min(e3_high, plane_high) - std::max(e3_low, plane_low);
              if (!(overlap > 0.0))
                continue;
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

Transformed transform_shoebox(const Experiment &e, const Shoebox &box,
                              const Vec3 &s1, double phi_calculated,
                              const GridSpec &spec) {
  // transform_shoebox_direct's answer, tested against it, faster two ways.
  // The face through pixel_cells: 49 per cent of the direct version's time
  // was epsilon_of for every subdivision. And the voxels: every valid voxel
  // ran twenty-five subdivisions times its planes, writing three grid arrays
  // each time, where the subdivisions of a pixel reach one to four cells --
  // so each reached cell is written once, with the count of subdivisions that
  // reached it.
  Transformed out;
  if (box.panel < 0 ||
      static_cast<std::size_t>(box.panel) >= e.detector.size()) {
    return out;
  }
  if (spec.subdivisions < 1)
    return out;
  const Panel &p = e.detector[static_cast<std::size_t>(box.panel)];
  const KabschFrame frame = kabsch_frame(e, s1);
  if (!frame.valid)
    return out;
  const int side = spec.side();
  const double span_d = spec.half_width * spec.sigma_d;
  const double span_m = spec.half_width * spec.sigma_m;
  if (!(span_d > 0.0) || !(span_m > 0.0))
    return out;
  const double step_d = 2.0 * span_d / static_cast<double>(side);
  const double step_m = 2.0 * span_m / static_cast<double>(side);
  const double share =
      1.0 / static_cast<double>(spec.subdivisions * spec.subdivisions);
  const std::size_t slice =
      static_cast<std::size_t>(side) * static_cast<std::size_t>(side);

  out.data.assign(spec.size(), 0.0);
  out.background.assign(spec.size(), 0.0);
  out.coverage.assign(spec.size(), 0.0);
  double inside = 0.0;
  double total = 0.0;

  const PixelCells cells = pixel_cells(e, frame, p, box, spec.subdivisions,
                                       side, span_d, step_d, phi_calculated);
  const std::int32_t nx = box.nx(), ny = box.ny();
  for (std::int32_t z = 0; z < box.nz(); ++z) {
    const Planes planes =
        plane_weights(e, frame, box, z, side, span_m, step_m, phi_calculated);
    if (!planes.spans)
      continue;
    for (std::int32_t y = 0; y < ny; ++y) {
      for (std::int32_t x = 0; x < nx; ++x) {
        const std::size_t at = box.at(x, y, z);
        // Validity, not the whole mask: a bad pixel keeps its region bit.
        if ((box.mask[at] & shoebox_mask::kValid) == 0)
          continue;
        const double count = static_cast<double>(box.data[at]);
        const double back = static_cast<double>(box.background[at]);
        const bool foreground = (box.mask[at] & shoebox_mask::kForeground) != 0;
        const double signal = std::max(count - back, 0.0);
        if (foreground)
          total += signal;
        const std::size_t k =
            static_cast<std::size_t>(y) * static_cast<std::size_t>(nx) +
            static_cast<std::size_t>(x);
        for (std::uint32_t q = cells.start[k]; q < cells.start[k + 1]; ++q) {
          const double reach = share * static_cast<double>(cells.hits[q]);
          for (const auto &plane : planes.planes) {
            const double weight = reach * plane.second;
            const std::size_t into =
                static_cast<std::size_t>(plane.first) * slice + cells.cell[q];
            out.data[into] += weight * count;
            out.background[into] += weight * back;
            out.coverage[into] += weight;
            if (foreground)
              inside += weight * signal;
          }
        }
      }
    }
  }
  out.outside = total > 0.0 ? 1.0 - inside / total : 1.0;
  out.valid = true;
  return out;
}

namespace {

//: Where a coordinate sits in a division, in units of cells, measured from the
//: centre of the first cell. Cell centres are then at 0, 1, 2, ...
double cell_position(double value, double size, int divisions) {
  const double n = static_cast<double>(std::max(divisions, 1));
  const double extent = size > 0.0 ? size : 1.0;
  return std::clamp(value / extent * n - 0.5, 0.0, n - 1.0);
}

} // namespace

std::size_t ReferenceProfiles::index_of(std::size_t which_panel, int block,
                                        int j, int i) const {
  const std::size_t d = static_cast<std::size_t>(std::max(divisions, 1));
  const std::size_t b = static_cast<std::size_t>(std::max(blocks, 1));
  const std::size_t panel = std::min(which_panel, panels - 1);
  const std::size_t bl =
      static_cast<std::size_t>(std::clamp(block, 0, std::max(blocks, 1) - 1));
  const std::size_t jj =
      static_cast<std::size_t>(std::clamp(j, 0, std::max(divisions, 1) - 1));
  const std::size_t ii =
      static_cast<std::size_t>(std::clamp(i, 0, std::max(divisions, 1) - 1));
  return ((panel * b + bl) * d + jj) * d + ii;
}

std::vector<double>
profile_on_pixels_direct(const Experiment &e, const Shoebox &box,
                         const Vec3 &s1, double phi_calculated,
                         const GridSpec &spec,
                         const std::vector<double> &reference) {
  std::vector<double> out;
  if (box.panel < 0 ||
      static_cast<std::size_t>(box.panel) >= e.detector.size()) {
    return out;
  }
  if (reference.size() != spec.size() || spec.subdivisions < 1)
    return out;
  const Panel &p = e.detector[static_cast<std::size_t>(box.panel)];
  const KabschFrame frame = kabsch_frame(e, s1);
  if (!frame.valid)
    return out;

  const int side = spec.side();
  const double span_d = spec.half_width * spec.sigma_d;
  const double span_m = spec.half_width * spec.sigma_m;
  if (!(span_d > 0.0) || !(span_m > 0.0))
    return out;
  const double step_d = 2.0 * span_d / static_cast<double>(side);
  const double step_m = 2.0 * span_m / static_cast<double>(side);
  const double osc = Scan::radians(e.scan.osc_width);
  const double share =
      1.0 / static_cast<double>(spec.subdivisions * spec.subdivisions);

  out.assign(box.size(), 0.0);

  // The subdivided face once, as the forward transform does: eps1 and eps2
  // depend only on where a subdivision sits on the detector.
  const std::int32_t fine_x = box.nx() * spec.subdivisions;
  const std::int32_t fine_y = box.ny() * spec.subdivisions;
  std::vector<int> face_i1(static_cast<std::size_t>(fine_x) *
                           static_cast<std::size_t>(fine_y));
  std::vector<int> face_i2(face_i1.size());
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
      const std::size_t at =
          static_cast<std::size_t>(fy) * static_cast<std::size_t>(fine_x) +
          static_cast<std::size_t>(fx);
      face_i1[at] = static_cast<int>(std::floor((eps.e1 + span_d) / step_d));
      face_i2[at] = static_cast<int>(std::floor((eps.e2 + span_d) / step_d));
    }
  }

  for (std::int32_t z = 0; z < box.nz(); ++z) {
    const double image = static_cast<double>(box.bbox[4] + z);
    const double phi_low = Scan::radians(e.scan.osc_start) +
                           (image - static_cast<double>(e.scan.z_offset)) * osc;
    const double phi_high = phi_low + osc;
    double e3_low = Scan::degrees(frame.zeta * (phi_low - phi_calculated));
    double e3_high = Scan::degrees(frame.zeta * (phi_high - phi_calculated));
    if (e3_low > e3_high)
      std::swap(e3_low, e3_high);
    const double e3_span = e3_high - e3_low;
    if (!(e3_span > 0.0))
      continue;

    for (std::int32_t y = 0; y < box.ny(); ++y) {
      for (std::int32_t x = 0; x < box.nx(); ++x) {
        const std::size_t at = box.at(x, y, z);
        // Every voxel, including ones with no measurement in them: the profile
        // is geometry and does not depend on whether the detector recorded
        // anything there. What is missing is the count, not the shape.
        double value = 0.0;
        for (int sy = 0; sy < spec.subdivisions; ++sy) {
          for (int sx = 0; sx < spec.subdivisions; ++sx) {
            const std::size_t on_face =
                static_cast<std::size_t>(y * spec.subdivisions + sy) *
                    static_cast<std::size_t>(fine_x) +
                static_cast<std::size_t>(x * spec.subdivisions + sx);
            const int i1 = face_i1[on_face];
            const int i2 = face_i2[on_face];
            if (i1 < 0 || i1 >= side || i2 < 0 || i2 >= side)
              continue;
            const int j_low =
                static_cast<int>(std::floor((e3_low + span_m) / step_m));
            const int j_high =
                static_cast<int>(std::floor((e3_high + span_m) / step_m));
            for (int j = std::max(j_low, 0); j <= std::min(j_high, side - 1);
                 ++j) {
              const double plane_low =
                  -span_m + static_cast<double>(j) * step_m;
              const double plane_high = plane_low + step_m;
              const double overlap =
                  std::min(e3_high, plane_high) - std::max(e3_low, plane_low);
              if (!(overlap > 0.0))
                continue;
              value +=
                  share * (overlap / e3_span) * reference[spec.at(i1, i2, j)];
            }
          }
        }
        out[at] = value;
      }
    }
  }
  return out;
}

std::vector<double> profile_on_pixels(const Experiment &e, const Shoebox &box,
                                      const Vec3 &s1, double phi_calculated,
                                      const GridSpec &spec,
                                      const std::vector<double> &reference) {
  // profile_on_pixels_direct's answer, tested against it. The face through
  // pixel_cells, as transform_shoebox; and each image's eps3 planes blended
  // into one slice of the grid once, so that a voxel is its pixel's short list
  // of cells against that slice. It was 0.772 ms a box direct, 97 per cent of
  // profile fitting.
  std::vector<double> out;
  if (box.panel < 0 ||
      static_cast<std::size_t>(box.panel) >= e.detector.size()) {
    return out;
  }
  if (reference.size() != spec.size() || spec.subdivisions < 1)
    return out;
  const Panel &p = e.detector[static_cast<std::size_t>(box.panel)];
  const KabschFrame frame = kabsch_frame(e, s1);
  if (!frame.valid)
    return out;
  const int side = spec.side();
  const double span_d = spec.half_width * spec.sigma_d;
  const double span_m = spec.half_width * spec.sigma_m;
  if (!(span_d > 0.0) || !(span_m > 0.0))
    return out;
  const double step_d = 2.0 * span_d / static_cast<double>(side);
  const double step_m = 2.0 * span_m / static_cast<double>(side);
  const double share =
      1.0 / static_cast<double>(spec.subdivisions * spec.subdivisions);
  const std::int32_t nx = box.nx(), ny = box.ny();
  out.assign(box.size(), 0.0);

  const PixelCells cells = pixel_cells(e, frame, p, box, spec.subdivisions,
                                       side, span_d, step_d, phi_calculated);
  const std::size_t slice =
      static_cast<std::size_t>(side) * static_cast<std::size_t>(side);
  std::vector<double> blended(slice);
  for (std::int32_t z = 0; z < box.nz(); ++z) {
    const Planes planes =
        plane_weights(e, frame, box, z, side, span_m, step_m, phi_calculated);
    if (!planes.spans || planes.planes.empty())
      continue;
    std::fill(blended.begin(), blended.end(), 0.0);
    for (const auto &plane : planes.planes) {
      const double *grid =
          reference.data() + static_cast<std::size_t>(plane.first) * slice;
      for (std::size_t q = 0; q < slice; ++q)
        blended[q] += plane.second * grid[q];
    }
    for (std::int32_t y = 0; y < ny; ++y) {
      for (std::int32_t x = 0; x < nx; ++x) {
        const std::size_t k =
            static_cast<std::size_t>(y) * static_cast<std::size_t>(nx) +
            static_cast<std::size_t>(x);
        double value = 0.0;
        for (std::uint32_t q = cells.start[k]; q < cells.start[k + 1]; ++q) {
          value += static_cast<double>(cells.hits[q]) * blended[cells.cell[q]];
        }
        out[box.at(x, y, z)] = share * value;
      }
    }
  }
  return out;
}

ProfileFit fit_on_pixels(const Shoebox &box,
                         const std::vector<double> &pixel_profile, double gain,
                         int iterations) {
  ProfileFit out;
  if (pixel_profile.size() != box.size())
    return out;
  if (box.data.size() != box.size() || box.mask.size() != box.size())
    return out;
  if (box.background.size() != box.size())
    return out;

  // TWO SUMS, AND THE DIFFERENCE BETWEEN THEM IS THE POINT.
  //
  // `profile_sum` is over the WHOLE foreground, including voxels with no
  // measurement in them, because that is what the intensity is: the scale
  // times the whole profile. `seen` is over the ones that were measured, and
  // is what the fit is done against.
  //
  // Normalising by the seen part instead gives a reflection with a third of
  // its foreground in a module gap two thirds of its intensity -- which is the
  // one thing profile fitting exists to avoid. Leslie sections 6.7.2 and
  // 6.7.3: a fitted profile recovers a reflection whose pixels are missing or
  // saturated, and that is most of its value beyond the variance.
  double profile_sum = 0.0;
  double seen = 0.0;
  for (std::size_t i = 0; i < box.size(); ++i) {
    if ((box.mask[i] & shoebox_mask::kForeground) == 0)
      continue;
    profile_sum += pixel_profile[i];
    if (box.mask[i] & shoebox_mask::kValid)
      seen += pixel_profile[i];
  }
  if (!(profile_sum > 0.0) || !(seen > 0.0))
    return out;
  out.measured = seen / profile_sum;

  double scale = 0.0;
  for (std::size_t i = 0; i < box.size(); ++i) {
    if ((box.mask[i] & shoebox_mask::kForeground) == 0)
      continue;
    if ((box.mask[i] & shoebox_mask::kValid) == 0)
      continue;
    scale += static_cast<double>(box.data[i]) -
             static_cast<double>(box.background[i]);
  }
  scale /= profile_sum;

  double variance = 0.0;
  for (int round = 0; round < std::max(iterations, 1); ++round) {
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t i = 0; i < box.size(); ++i) {
      if ((box.mask[i] & shoebox_mask::kForeground) == 0)
        continue;
      if ((box.mask[i] & shoebox_mask::kValid) == 0)
        continue;
      const double b = static_cast<double>(box.background[i]);
      const double expected = b + std::max(scale, 0.0) * pixel_profile[i];
      const double v = std::max(gain * expected, gain * kLeastVariance);
      const double residual = static_cast<double>(box.data[i]) - b;
      numerator += pixel_profile[i] * residual / v;
      denominator += pixel_profile[i] * pixel_profile[i] / v;
    }
    if (!(denominator > 0.0))
      return out;
    scale = numerator / denominator;
    variance = 1.0 / denominator;
    out.iterations = round + 1;
  }
  out.intensity = scale * profile_sum;
  out.variance = variance * profile_sum * profile_sum;

  // Leslie equation 34: the fitted variance has TWO parts, the fit itself and
  // the background. The second is the same (m/n) I_bg that summation carries,
  // because the background was estimated from n pixels and subtracted from m
  // of them, and that uncertainty does not go away because the foreground was
  // weighted by a profile.
  //
  // Leaving it out made the fitted variance 0.32 of the summed one where DIALS
  // has 0.85 -- and 0.32 is below the floor Leslie section 6.6 derives, which
  // is about 0.5 for a typical profile. A ratio better than the theory allows
  // is not a better algorithm; it is a term that has been forgotten.
  {
    double foreground = 0.0, background_pixels = 0.0, background_sum = 0.0;
    for (std::size_t i = 0; i < box.size(); ++i) {
      if ((box.mask[i] & shoebox_mask::kValid) == 0)
        continue;
      if (box.mask[i] & shoebox_mask::kForeground) {
        foreground += 1.0;
        background_sum += static_cast<double>(box.background[i]);
      } else if (box.mask[i] & shoebox_mask::kBackground) {
        background_pixels += 1.0;
      }
    }
    if (background_pixels > 0.0) {
      out.variance += gain * (foreground / background_pixels) * background_sum;
    }
  }

  double n = 0.0, sp = 0.0, sd = 0.0, spp = 0.0, sdd = 0.0, spd = 0.0;
  for (std::size_t i = 0; i < box.size(); ++i) {
    if ((box.mask[i] & shoebox_mask::kForeground) == 0)
      continue;
    if ((box.mask[i] & shoebox_mask::kValid) == 0)
      continue;
    const double a = pixel_profile[i];
    const double b = static_cast<double>(box.data[i]) -
                     static_cast<double>(box.background[i]);
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
    if (va > 0.0 && vb > 0.0)
      out.correlation = cov / std::sqrt(va * vb);
  }
  out.valid = true;
  return out;
}

std::size_t ReferenceProfiles::region_of(const Panel &panel,
                                         std::size_t which_panel,
                                         double px_fast, double px_slow,
                                         double z) const {
  const double fast_size =
      static_cast<double>(std::max<std::int64_t>(panel.image_size[0], 1));
  const double slow_size =
      static_cast<double>(std::max<std::int64_t>(panel.image_size[1], 1));
  const int i = static_cast<int>(
      std::lround(cell_position(px_fast, fast_size, divisions)));
  const int j = static_cast<int>(
      std::lround(cell_position(px_slow, slow_size, divisions)));
  const int b = static_cast<int>(std::lround(cell_position(
      z - first_image, std::max(last_image - first_image, 1.0), blocks)));
  return index_of(which_panel, b, j, i);
}

std::vector<Neighbour> neighbours_of(const ReferenceProfiles &reference,
                                     const Panel &panel,
                                     std::size_t which_panel, double px_fast,
                                     double px_slow, double z) {
  const double fast_size =
      static_cast<double>(std::max<std::int64_t>(panel.image_size[0], 1));
  const double slow_size =
      static_cast<double>(std::max<std::int64_t>(panel.image_size[1], 1));
  const double span =
      std::max(reference.last_image - reference.first_image, 1.0);

  const double x = cell_position(px_fast, fast_size, reference.divisions);
  const double y = cell_position(px_slow, slow_size, reference.divisions);
  const double t =
      cell_position(z - reference.first_image, span, reference.blocks);

  // The cell below and the one above, in each division, with the weight
  // falling linearly between their centres. At an edge both land on the same
  // cell and the weights add back to one there, which is the fallback rather
  // than a special case.
  const auto pair = [](double v, int limit) {
    const int low = std::clamp(static_cast<int>(std::floor(v)), 0, limit - 1);
    const int high = std::clamp(low + 1, 0, limit - 1);
    const double up = (high == low) ? 0.0 : v - static_cast<double>(low);
    return std::array<std::pair<int, double>, 2>{
        std::pair<int, double>{low, 1.0 - up},
        std::pair<int, double>{high, up}};
  };

  const auto in_x = pair(x, std::max(reference.divisions, 1));
  const auto in_y = pair(y, std::max(reference.divisions, 1));
  const auto in_t = pair(t, std::max(reference.blocks, 1));

  std::vector<Neighbour> out;
  out.reserve(8);
  for (const auto &a : in_t) {
    for (const auto &b : in_y) {
      for (const auto &c : in_x) {
        const double weight = a.second * b.second * c.second;
        if (!(weight > 0.0))
          continue;
        const std::size_t region =
            reference.index_of(which_panel, a.first, b.first, c.first);
        // An edge can name the same cell twice; add rather than append, so the
        // weights still sum to one.
        bool merged = false;
        for (Neighbour &n : out) {
          if (n.region == region) {
            n.weight += weight;
            merged = true;
            break;
          }
        }
        if (!merged)
          out.push_back({region, weight});
      }
    }
  }
  return out;
}

std::vector<double> profile_at(const ReferenceProfiles &reference,
                               const Panel &panel, std::size_t which_panel,
                               double px_fast, double px_slow, double z) {
  std::vector<double> out(reference.spec.size(), 0.0);
  const std::vector<Neighbour> near =
      neighbours_of(reference, panel, which_panel, px_fast, px_slow, z);
  double total = 0.0;
  for (const Neighbour &n : near) {
    if (n.region >= reference.profile.size())
      continue;
    const std::vector<double> &p = reference.profile[n.region];
    if (p.size() != out.size())
      continue;
    for (std::size_t i = 0; i < out.size(); ++i)
      out[i] += n.weight * p[i];
    total += n.weight;
  }
  if (!(total > 0.0))
    return out;
  double sum = 0.0;
  for (double v : out)
    sum += v;
  if (sum > 0.0) {
    for (double &v : out)
      v /= sum;
  }
  return out;
}

ReferenceProfiles make_reference(const GridSpec &spec, int divisions,
                                 int blocks, std::size_t panels,
                                 double first_image, double last_image) {
  ReferenceProfiles out;
  out.spec = spec;
  out.divisions = std::max(divisions, 1);
  out.blocks = std::max(blocks, 1);
  out.panels = std::max<std::size_t>(panels, 1);
  out.first_image = first_image;
  out.last_image = std::max(last_image, first_image + 1.0);
  out.profile.assign(out.region_count(), std::vector<double>(spec.size(), 0.0));
  out.spots.assign(out.region_count(), 0);
  return out;
}

bool add_reference(ReferenceProfiles *reference, std::size_t region,
                   const Transformed &t) {
  if (reference == nullptr || !t.valid)
    return false;
  if (region >= reference->profile.size())
    return false;
  if (t.data.size() != reference->spec.size())
    return false;

  // Normalised to unit sum, so the profile is the average shape rather than
  // the average spot: without this a single strong reflection outvotes a
  // hundred weak ones and the profile is that reflection's.
  double sum = 0.0;
  for (std::size_t i = 0; i < t.data.size(); ++i) {
    sum += t.data[i] - t.background[i];
  }
  if (!(sum > 0.0))
    return false;
  std::vector<double> &into = reference->profile[region];
  for (std::size_t i = 0; i < t.data.size(); ++i) {
    into[i] += (t.data[i] - t.background[i]) / sum;
  }
  ++reference->spots[region];
  return true;
}

void finalise_reference(ReferenceProfiles *reference, std::size_t least) {
  if (reference == nullptr)
    return;
  const std::size_t n = reference->spec.size();

  // The whole-detector average, for regions that saw too few spots to have a
  // profile of their own. A region at the corner of a detector may have almost
  // nothing in it, and an empty profile fits nothing at all.
  std::vector<double> everything(n, 0.0);
  std::size_t everything_spots = 0;
  for (std::size_t r = 0; r < reference->profile.size(); ++r) {
    if (reference->spots[r] == 0)
      continue;
    for (std::size_t i = 0; i < n; ++i)
      everything[i] += reference->profile[r][i];
    everything_spots += reference->spots[r];
  }
  double whole = 0.0;
  for (double v : everything)
    whole += v;
  if (whole > 0.0) {
    for (double &v : everything)
      v /= whole;
  }

  for (std::size_t r = 0; r < reference->profile.size(); ++r) {
    if (reference->spots[r] < least) {
      reference->profile[r] = everything;
      continue;
    }
    double sum = 0.0;
    for (double v : reference->profile[r])
      sum += v;
    if (!(sum > 0.0)) {
      reference->profile[r] = everything;
      continue;
    }
    for (double &v : reference->profile[r])
      v /= sum;
  }
  (void)everything_spots;
  reference->finalised = true;
}

ProfileFit fit_profile(const std::vector<double> &reference,
                       const Transformed &t, double gain, int iterations) {
  ProfileFit out;
  if (!t.valid || reference.size() != t.data.size())
    return out;

  double profile_sum = 0.0;
  for (double v : reference)
    profile_sum += v;
  if (!(profile_sum > 0.0))
    return out;

  // A first guess from the plain sum, to weight the first round with.
  double scale = 0.0;
  for (std::size_t i = 0; i < t.data.size(); ++i) {
    if (t.coverage[i] > 0.0)
      scale += t.data[i] - t.background[i];
  }
  scale /= profile_sum;

  double variance = 0.0;
  for (int round = 0; round < std::max(iterations, 1); ++round) {
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t i = 0; i < t.data.size(); ++i) {
      // A grid point no pixel reached is not a measurement and cannot be
      // fitted to. Counting it as an observed zero pulls every intensity down.
      if (!(t.coverage[i] > 0.0))
        continue;
      // The Poisson variance of what the model says should be there. Using the
      // OBSERVED count instead biases the fit: a point that happened to record
      // nothing would be given infinite weight.
      //
      // THE FLOOR IS THE POINT. It was an epsilon, 1e-6, and a grid point whose
      // expected counts are near zero then carries a weight of a million and
      // the fit does what those few points say. That is not hypothetical: 79 of
      // 19771 reflections came back fitted with the opposite sign to their
      // summed intensity, their sums having a median of 14845 -- the strongest
      // reflections in the dataset, fitted large and negative. Raising the
      // floor took that to none, and the agreement with DIALS from 0.379 to
      // 0.968, on the same data at the same resolution.
      //
      // Such points are ordinary here rather than pathological: the background
      // at a grid point is shared out in proportion to how much of a pixel
      // reached it, so a point at the edge of a spot's footprint has a
      // background of a ten-thousandth of a count and an expected value to
      // match.
      //
      // 0.05 counts is the least variance a measurement of anything is allowed
      // to have. It is arbitrary in the way a floor must be, and it is four
      // orders of magnitude away from the value that failed.
      //
      // The clamp on the scale below is NOT what fixed this -- it was measured
      // separately and moved 79 flips to 78. It stays because a negative
      // expected count is not a Poisson mean whatever else is true.
      const double expected =
          t.background[i] + std::max(scale, 0.0) * reference[i];
      const double v = std::max(gain * expected, gain * kLeastVariance);
      numerator += reference[i] * (t.data[i] - t.background[i]) / v;
      denominator += reference[i] * reference[i] / v;
    }
    if (!(denominator > 0.0))
      return out;
    scale = numerator / denominator;
    variance = 1.0 / denominator;
    out.iterations = round + 1;
  }

  out.intensity = scale * profile_sum;
  out.variance = variance * profile_sum * profile_sum;

  // How well the profile actually describes this reflection.
  double n = 0.0, sp = 0.0, sd = 0.0, spp = 0.0, sdd = 0.0, spd = 0.0;
  for (std::size_t i = 0; i < t.data.size(); ++i) {
    if (!(t.coverage[i] > 0.0))
      continue;
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
    if (va > 0.0 && vb > 0.0)
      out.correlation = cov / std::sqrt(va * vb);
  }
  out.valid = true;
  return out;
}

} // namespace mxi
