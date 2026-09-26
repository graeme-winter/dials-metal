#include "mask.hh"

#include <algorithm>
#include <cmath>

namespace mxi {

namespace {

// Where a ray offset from s1 by (t1, t2) radians in the frame's own plane
// lands on the panel, in pixels. The offsets are small -- three sigma is a
// hundredth of a degree here -- so moving along the tangent plane and
// renormalising is accurate well past the pixel it is used to find.
bool offset_pixel(const Panel &p, const KabschFrame &frame, double t1,
                  double t2, double *px_fast, double *px_slow) {
  const double length = frame.s1.norm();
  if (!(length > 0.0))
    return false;
  Vec3 s = frame.s1 + (frame.e1 * t1 + frame.e2 * t2) * length;
  const double n = s.norm();
  if (!(n > 0.0))
    return false;
  s = s * (length / n);
  const auto hit = p.intersect(s);
  if (!hit.has_value())
    return false;
  *px_fast = hit->first;
  *px_slow = hit->second;
  return true;
}

} // namespace

const char *describe(BoxRejection why) {
  switch (why) {
  case BoxRejection::kNone:
    return "kept";
  case BoxRejection::kNoPanel:
    return "no such panel";
  case BoxRejection::kNoFrame:
    return "no Kabsch frame";
  case BoxRejection::kSmallZeta:
    return "zeta below the cut";
  case BoxRejection::kNoIntersection:
    return "a corner misses the detector";
  case BoxRejection::kTooManyImages:
    return "spans too many images";
  case BoxRejection::kOffDetector:
    return "off the detector or the scan";
  }
  return "unknown";
}

bool integration_bbox(const Experiment &e, const Prediction &p,
                      const MaskOptions &options, std::int32_t bbox[6],
                      BoxRejection *why) {
  const auto refuse = [&](BoxRejection reason) {
    if (why != nullptr)
      *why = reason;
    return false;
  };
  if (why != nullptr)
    *why = BoxRejection::kNone;
  if (p.panel >= e.detector.size())
    return refuse(BoxRejection::kNoPanel);
  const Panel &panel = e.detector[p.panel];
  const Vec3 s1 = p.s1;
  const KabschFrame frame = kabsch_frame(e, s1);
  if (!frame.valid)
    return refuse(BoxRejection::kNoFrame);
  if (std::fabs(frame.zeta) < options.min_zeta) {
    return refuse(BoxRejection::kSmallZeta);
  }

  // The BOX is wider than the foreground region, to hold background. The mask
  // below still marks only the n_sigma region as foreground.
  const double scale = options.box_scale > 1.0 ? options.box_scale : 1.0;
  const double t = Scan::radians(options.n_sigma * options.sigma_d) * scale;
  if (!(t > 0.0))
    return refuse(BoxRejection::kNoFrame);

  // The four corners of the region in the tangent plane. The mapping is
  // smooth and monotonic over an offset this small, so the corners bound it.
  double low_fast = 0.0, high_fast = 0.0, low_slow = 0.0, high_slow = 0.0;
  bool first = true;
  for (int i = 0; i < 4; ++i) {
    const double t1 = (i & 1) ? t : -t;
    const double t2 = (i & 2) ? t : -t;
    double px_fast = 0.0, px_slow = 0.0;
    if (!offset_pixel(panel, frame, t1, t2, &px_fast, &px_slow)) {
      return refuse(BoxRejection::kNoIntersection);
    }
    if (first) {
      low_fast = high_fast = px_fast;
      low_slow = high_slow = px_slow;
      first = false;
    } else {
      low_fast = std::fmin(low_fast, px_fast);
      high_fast = std::fmax(high_fast, px_fast);
      low_slow = std::fmin(low_slow, px_slow);
      high_slow = std::fmax(high_slow, px_slow);
    }
  }

  // The rotation half-width, which is sigma_M divided by zeta: a reflection
  // crossing the sphere obliquely takes longer to do it.
  // In the rotation direction the box follows the foreground: DIALS' image
  // extents agree with ours already, 99.3 per cent of the time, so widening
  // here would break what is right to fix what is not.
  const double half =
      Scan::radians(options.n_sigma * options.sigma_m) / std::fabs(frame.zeta);
  const double z_low = e.scan.z_from_phi(p.phi - half);
  const double z_high = e.scan.z_from_phi(p.phi + half);

  bbox[0] = static_cast<std::int32_t>(std::floor(low_fast));
  bbox[1] = static_cast<std::int32_t>(std::ceil(high_fast));
  bbox[2] = static_cast<std::int32_t>(std::floor(low_slow));
  bbox[3] = static_cast<std::int32_t>(std::ceil(high_slow));
  bbox[4] = static_cast<std::int32_t>(std::floor(std::fmin(z_low, z_high)));
  bbox[5] = static_cast<std::int32_t>(std::ceil(std::fmax(z_low, z_high)));

  if (bbox[5] - bbox[4] > options.max_images) {
    return refuse(BoxRejection::kTooManyImages);
  }

  // Clipped to what exists. A box hanging off the edge of the detector or the
  // end of the scan is kept, cut down to the part that is real.
  bbox[0] = std::max<std::int32_t>(bbox[0], 0);
  bbox[2] = std::max<std::int32_t>(bbox[2], 0);
  bbox[1] = std::min<std::int32_t>(
      bbox[1], static_cast<std::int32_t>(panel.image_size[0]));
  bbox[3] = std::min<std::int32_t>(
      bbox[3], static_cast<std::int32_t>(panel.image_size[1]));
  bbox[4] = std::max<std::int32_t>(bbox[4], 0);
  bbox[5] = std::min<std::int32_t>(
      bbox[5], static_cast<std::int32_t>(e.scan.num_images()));
  if (!(bbox[1] > bbox[0] && bbox[3] > bbox[2] && bbox[5] > bbox[4])) {
    return refuse(BoxRejection::kOffDetector);
  }
  return true;
}

bool build_shoebox(const Experiment &e, const Prediction &p,
                   const MaskOptions &options, Shoebox *box,
                   BoxRejection *why) {
  if (box == nullptr)
    return false;
  std::int32_t bbox[6];
  if (!integration_bbox(e, p, options, bbox, why))
    return false;
  const Panel &panel = e.detector[p.panel];
  const KabschFrame frame = kabsch_frame(e, p.s1);

  box->panel = static_cast<std::int32_t>(p.panel);
  for (int i = 0; i < 6; ++i)
    box->bbox[i] = bbox[i];
  box->flag = 2;
  const std::size_t n = box->size();
  box->data.assign(n, 0.0f);
  box->background.assign(n, 0.0f);
  box->mask.assign(n, 0);

  const double width = Scan::radians(e.scan.osc_width);

  // eps1 and eps2 depend only on where a pixel is on the detector, and eps3
  // only on which image it is. Computing all three per voxel evaluates the
  // detector mapping nz times over for every pixel -- and that mapping has a
  // parallax correction with an exp() in it. On a thirty degree sweep this
  // loop was 63 per cent of the whole integration, more than reading and
  // decompressing the images put together.
  //
  // So the face is computed once and the images once: nx*ny + nz evaluations
  // instead of nx*ny*nz. The arithmetic is unchanged and so is every mask.
  const std::int32_t face = box->nx() * box->ny();
  std::vector<double> eps1(static_cast<std::size_t>(face));
  std::vector<double> eps2(static_cast<std::size_t>(face));
  for (std::int32_t y = 0; y < box->ny(); ++y) {
    for (std::int32_t x = 0; x < box->nx(); ++x) {
      const Epsilon eps =
          epsilon_of(e, frame, panel, static_cast<double>(bbox[0] + x) + 0.5,
                     static_cast<double>(bbox[2] + y) + 0.5, p.phi, p.phi);
      const std::size_t at =
          static_cast<std::size_t>(y) * static_cast<std::size_t>(box->nx()) +
          static_cast<std::size_t>(x);
      eps1[at] = eps.e1;
      eps2[at] = eps.e2;
    }
  }
  std::vector<double> eps3(static_cast<std::size_t>(box->nz()));
  for (std::int32_t z = 0; z < box->nz(); ++z) {
    const double image = static_cast<double>(bbox[4] + z) + 0.5;
    const double phi = Scan::radians(e.scan.osc_start) +
                       (image - static_cast<double>(e.scan.z_offset)) * width;
    eps3[static_cast<std::size_t>(z)] =
        Scan::degrees(frame.zeta * (phi - p.phi));
  }

  for (std::int32_t z = 0; z < box->nz(); ++z) {
    for (std::int32_t y = 0; y < box->ny(); ++y) {
      for (std::int32_t x = 0; x < box->nx(); ++x) {
        Epsilon eps;
        const std::size_t on_face =
            static_cast<std::size_t>(y) * static_cast<std::size_t>(box->nx()) +
            static_cast<std::size_t>(x);
        eps.e1 = eps1[on_face];
        eps.e2 = eps2[on_face];
        eps.e3 = eps3[static_cast<std::size_t>(z)];
        // In sigmas, which is the only scale on which the three directions
        // are comparable: eps1 and eps2 are measured against sigma_D and eps3
        // against sigma_M, and here they differ by a factor of four.
        const double u1 = eps.e1 / options.sigma_d;
        const double u2 = eps.e2 / options.sigma_d;
        const double u3 = eps.e3 / options.sigma_m;
        const bool inside = options.shape == RegionShape::kEllipsoid
                                ? (u1 * u1 + u2 * u2 + u3 * u3) <=
                                      options.n_sigma * options.n_sigma
                                : (std::fabs(u1) <= options.n_sigma &&
                                   std::fabs(u2) <= options.n_sigma &&
                                   std::fabs(u3) <= options.n_sigma);
        // Every voxel is valid -- nothing here has read an image, so nothing
        // is known to be bad -- and each is either in the region or around it.
        box->mask[box->at(x, y, z)] = static_cast<std::uint8_t>(
            shoebox_mask::kValid |
            (inside ? shoebox_mask::kForeground : shoebox_mask::kBackground));
      }
    }
  }
  return true;
}

} // namespace mxi
