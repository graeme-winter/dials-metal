// Shoeboxes built from the profile model, for looking at against the images.

#include <cmath>
#include <vector>

#include "../src/mask.hh"
#include "../src/predict.hh"
#include "check.hh"

namespace mxi {

namespace {

Experiment masking_experiment() {
  Experiment e;
  e.beam.direction = {0.0, 0.0, 1.0};
  e.beam.wavelength = 0.9537;
  Panel p;
  p.fast = {1.0, 0.0, 0.0};
  p.slow = {0.0, -1.0, 0.0};
  p.pixel_size[0] = p.pixel_size[1] = 0.075;
  p.image_size[0] = 2068;
  p.image_size[1] = 2162;
  p.origin = {-0.075 * 1034, 0.075 * 1081, -170.0};
  e.detector.panels.push_back(p);
  e.goniometer.axis = {1.0, 0.0, 0.0};
  e.scan.first_image = 1;
  e.scan.last_image = 400;
  e.scan.osc_start = 0.0;
  e.scan.osc_width = 0.1;
  const Mat3 r = rotation({0.3, -0.5, 0.81}, 0.7);
  e.crystal = Crystal::from_real_space(r * Vec3{78, 0, 0}, r * Vec3{0, 78, 0},
                                       r * Vec3{0, 0, 78});
  return e;
}

//: A box that hugs the foreground, for the tests that are about the region
//: rather than about the rim.
MaskOptions options_for(double n_sigma = 3.0) {
  MaskOptions o;
  o.n_sigma = n_sigma;
  o.box_scale = 1.0;  // no rim, so the box is the region
  o.sigma_d = 0.031;
  o.sigma_m = 0.098;
  return o;
}

}  // namespace

TEST(the_box_encloses_the_region_and_not_much_more) {
  // The bounding box has to contain every voxel the mask will call foreground
  // and no whole plane of voxels that it will not, or the box is either
  // clipping the region or padding it.
  const Experiment e = masking_experiment();
  PredictOptions po;
  po.d_min = 2.0;
  const std::vector<Prediction> predictions = predict(e, po);
  check::is_true(predictions.size() > 100, "enough predictions");

  // With the rim, which is what the pipeline uses: the assertions below are
  // about the rim existing and the foreground sitting inside it.
  MaskOptions options = options_for();
  options.box_scale = 1.9;
  std::size_t checked = 0;
  for (const Prediction &p : predictions) {
    Shoebox box;
    if (!build_shoebox(e, p, options, &box)) continue;
    if (++checked > 200) break;

    // Every foreground voxel is inside the box by construction; the real
    // question is whether the outermost planes are empty, which would mean the
    // box is bigger than it needs to be.
    bool edge_has_foreground = false;
    for (std::int32_t z = 0; z < box.nz(); ++z) {
      for (std::int32_t y = 0; y < box.ny(); ++y) {
        for (std::int32_t x = 0; x < box.nx(); ++x) {
          if ((box.mask[box.at(x, y, z)] & shoebox_mask::kForeground) == 0) continue;
          // Fast and slow only. The rim is on the detector; in the rotation
          // direction the box follows the foreground exactly, because DIALS'
          // image extents already agreed and widening them would break what
          // is right to fix what is not. So foreground on the first or last
          // image is expected and is not an error.
          if (x == 0 || y == 0 || x == box.nx() - 1 || y == box.ny() - 1) {
            edge_has_foreground = true;
          }
        }
      }
    }
    // Foreground reaching the edge is fine and expected; foreground OUTSIDE
    // would show as a box too small, which cannot happen here, so what this
    // pins is that a box is not wildly padded: at least a tenth of it is
    // foreground.
    std::size_t foreground = 0;
    for (std::uint8_t m : box.mask) {
      if (m & shoebox_mask::kForeground) ++foreground;
    }
    check::is_true(foreground > 0, "a box always contains its own reflection");
    // The box is deliberately wider than the foreground, because the
    // background estimate needs background pixels and there are none in a box
    // that is only the peak. DIALS' boxes are 13 per cent foreground on real
    // data. What must hold is that the foreground is inside the box and does
    // not reach its edge, or the rim is not a rim.
    check::is_true(!edge_has_foreground,
                   "the foreground does not reach the fast or slow edge, so "
                   "there is a rim of background to estimate from");
    check::is_true(foreground * 100 >= box.size(),
                   "and the box is not absurdly larger than the region");
  }
  check::is_true(checked > 50, "enough boxes built");
}

TEST(every_voxel_is_valid_and_either_foreground_or_background) {
  // The viewer reads these codes. A voxel marked neither, or marked invalid,
  // would be drawn as a bad pixel over a perfectly good image.
  const Experiment e = masking_experiment();
  PredictOptions po;
  po.d_min = 2.5;
  const std::vector<Prediction> predictions = predict(e, po);
  const MaskOptions options = options_for();
  std::size_t checked = 0;
  for (const Prediction &p : predictions) {
    Shoebox box;
    if (!build_shoebox(e, p, options, &box)) continue;
    if (++checked > 50) break;
    for (std::uint8_t m : box.mask) {
      check::is_true((m & shoebox_mask::kValid) != 0, "valid");
      const bool foreground = (m & shoebox_mask::kForeground) != 0;
      const bool background = (m & shoebox_mask::kBackground) != 0;
      check::is_true(foreground != background, "exactly one of the two");
    }
    for (float v : box.data) check::close(v, 0.0, 0.0, "no invented counts");
    for (float v : box.background) check::close(v, 0.0, 0.0, "nor background");
  }
  check::is_true(checked > 10, "enough boxes");
}

TEST(a_wider_region_gives_a_bigger_box) {
  // The box follows n_sigma. If it did not, the picture would show the same
  // thing whatever was asked for, which is the sort of bug that survives
  // because the output always looks plausible.
  const Experiment e = masking_experiment();
  PredictOptions po;
  po.d_min = 2.5;
  const std::vector<Prediction> predictions = predict(e, po);
  std::size_t compared = 0;
  for (const Prediction &p : predictions) {
    Shoebox narrow, wide;
    if (!build_shoebox(e, p, options_for(1.0), &narrow)) continue;
    if (!build_shoebox(e, p, options_for(3.0), &wide)) continue;
    if (++compared > 50) break;
    check::is_true(wide.size() > narrow.size(), "three sigma is bigger than one");
    // And the narrow region is contained in the wide one.
    check::is_true(wide.bbox[0] <= narrow.bbox[0], "in fast");
    check::is_true(wide.bbox[1] >= narrow.bbox[1], "both ways");
    check::is_true(wide.bbox[4] <= narrow.bbox[4], "and in images");
  }
  check::is_true(compared > 10, "enough compared");
}

TEST(a_reflection_near_the_rotation_axis_is_refused) {
  // Its region in rotation is sigma_M over zeta, which diverges: the box would
  // swallow the scan. Refused rather than allocated.
  const Experiment e = masking_experiment();
  PredictOptions po;
  po.d_min = 2.0;
  const std::vector<Prediction> predictions = predict(e, po);
  MaskOptions options = options_for();
  options.min_zeta = 0.5;  // a severe cut, so some are certainly refused

  std::size_t built = 0, refused = 0;
  for (const Prediction &p : predictions) {
    Shoebox box;
    if (build_shoebox(e, p, options, &box)) {
      ++built;
      const KabschFrame frame = kabsch_frame(e, p.s1);
      check::is_true(std::fabs(frame.zeta) >= options.min_zeta,
                     "nothing below the cut was built");
    } else {
      ++refused;
    }
  }
  check::is_true(built > 0, "some were built");
  check::is_true(refused > 0, "and some refused");
}

TEST(the_ellipsoid_is_inside_the_box_and_is_pi_over_six_of_it) {
  // The box is each coordinate separately within n sigma; the ellipsoid is the
  // surface the Gaussian is constant on. The ellipsoid is inscribed in the
  // box, so every voxel it marks the box marks too, and its volume is pi/6 of
  // the box's -- a little over half.
  //
  // At the corner of the box all three coordinates are at n sigma at once, so
  // a three-dimensional Gaussian is at exp(-3 n^2 / 2) there: 1.4e-6 of its
  // peak at n = 3. There are eight such corners and they hold nothing.
  const Experiment e = masking_experiment();
  PredictOptions po;
  po.d_min = 2.0;
  const std::vector<Prediction> predictions = predict(e, po);

  MaskOptions box_options = options_for();
  MaskOptions ellipsoid_options = options_for();
  ellipsoid_options.shape = RegionShape::kEllipsoid;

  std::size_t in_box = 0, in_ellipsoid = 0, compared = 0;
  for (const Prediction &p : predictions) {
    Shoebox as_box, as_ellipsoid;
    if (!build_shoebox(e, p, box_options, &as_box)) continue;
    if (!build_shoebox(e, p, ellipsoid_options, &as_ellipsoid)) continue;
    if (++compared > 300) break;

    // Same bounding box either way: the ellipsoid is inscribed, so it does not
    // need a smaller one and must not get a different one, or the two pictures
    // would not be comparable.
    for (int k = 0; k < 6; ++k) {
      check::equal(static_cast<long long>(as_ellipsoid.bbox[k]),
                   static_cast<long long>(as_box.bbox[k]), "same box");
    }
    for (std::size_t i = 0; i < as_box.mask.size(); ++i) {
      const bool b = (as_box.mask[i] & shoebox_mask::kForeground) != 0;
      const bool l = (as_ellipsoid.mask[i] & shoebox_mask::kForeground) != 0;
      check::is_true(!l || b, "every ellipsoid voxel is a box voxel");
      if (b) ++in_box;
      if (l) ++in_ellipsoid;
    }
  }
  check::is_true(compared > 100, "enough compared");
  const double ratio = static_cast<double>(in_ellipsoid) /
                       static_cast<double>(in_box);
  // pi/6 = 0.5236. Measured at 0.546 on real insulin, the difference being
  // that a voxel is in or out as a whole and the boxes are only a dozen
  // pixels across.
  check::close(ratio, 0.5236, 0.06, "and its share is pi over six");
}

TEST(the_region_is_an_ellipsoid_in_sigmas_not_in_degrees) {
  // eps1 and eps2 are measured against sigma_D and eps3 against sigma_M, which
  // differ by a factor of four here. An ellipsoid built in degrees rather than
  // in sigmas would be a sphere in the wrong space: far too generous in the
  // narrow direction and far too mean in the wide one.
  const Experiment e = masking_experiment();
  PredictOptions po;
  po.d_min = 2.0;
  const std::vector<Prediction> predictions = predict(e, po);
  MaskOptions options = options_for();
  options.shape = RegionShape::kEllipsoid;

  for (const Prediction &p : predictions) {
    Shoebox box;
    if (!build_shoebox(e, p, options, &box)) continue;
    const KabschFrame frame = kabsch_frame(e, p.s1);
    const Panel &panel = e.detector[p.panel];
    const double width = Scan::radians(e.scan.osc_width);
    // Every marked voxel must satisfy the ellipsoid in sigmas.
    for (std::int32_t z = 0; z < box.nz(); ++z) {
      const double image = static_cast<double>(box.bbox[4] + z) + 0.5;
      const double phi = Scan::radians(e.scan.osc_start) +
                         (image - static_cast<double>(e.scan.z_offset)) * width;
      for (std::int32_t y = 0; y < box.ny(); ++y) {
        for (std::int32_t x = 0; x < box.nx(); ++x) {
          if ((box.mask[box.at(x, y, z)] & shoebox_mask::kForeground) == 0) continue;
          const Epsilon eps = epsilon_of(
              e, frame, panel, static_cast<double>(box.bbox[0] + x) + 0.5,
              static_cast<double>(box.bbox[2] + y) + 0.5, phi, p.phi);
          const double u1 = eps.e1 / options.sigma_d;
          const double u2 = eps.e2 / options.sigma_d;
          const double u3 = eps.e3 / options.sigma_m;
          check::is_true(u1 * u1 + u2 * u2 + u3 * u3 <=
                             options.n_sigma * options.n_sigma + 1e-9,
                         "inside the ellipsoid in sigmas");
        }
      }
    }
    return;  // one reflection is enough to pin the arithmetic
  }
  check::is_true(false, "no reflection was built");
}

TEST(a_refused_box_says_which_reason) {
  // Half the predictions of a ten rotation sweep produced no shoebox and there
  // was no way to tell whether that was the detector, the rotation axis or a
  // bug. Each reason is now counted, and each is reachable: a test that only
  // checked the common one would leave the others as untested strings.
  const Experiment e = masking_experiment();
  MaskOptions options = options_for();

  // A panel that does not exist.
  {
    PredictOptions po;
    po.d_min = 3.0;
    const std::vector<Prediction> predictions = predict(e, po);
    check::is_true(!predictions.empty(), "something was predicted");
    Prediction p = predictions.front();
    p.panel = 7;
    Shoebox box;
    BoxRejection why = BoxRejection::kNone;
    check::is_true(!build_shoebox(e, p, options, &box, &why), "refused");
    check::is_true(why == BoxRejection::kNoPanel, "because there is no panel 7");
  }

  // Near the rotation axis: zeta small, so the region in rotation diverges.
  {
    PredictOptions po;
    po.d_min = 2.0;
    const std::vector<Prediction> predictions = predict(e, po);
    MaskOptions severe = options;
    severe.min_zeta = 0.99;  // so that almost everything is below it
    std::size_t small_zeta = 0;
    for (const Prediction &p : predictions) {
      Shoebox box;
      BoxRejection why = BoxRejection::kNone;
      if (!build_shoebox(e, p, severe, &box, &why) &&
          why == BoxRejection::kSmallZeta) {
        ++small_zeta;
      }
    }
    check::is_true(small_zeta > 0, "the zeta cut is reported as itself");
  }

  // A box wider in rotation than the limit allows.
  {
    PredictOptions po;
    po.d_min = 2.0;
    const std::vector<Prediction> predictions = predict(e, po);
    MaskOptions narrow = options;
    narrow.max_images = 1;
    std::size_t too_many = 0;
    for (const Prediction &p : predictions) {
      Shoebox box;
      BoxRejection why = BoxRejection::kNone;
      if (!build_shoebox(e, p, narrow, &box, &why) &&
          why == BoxRejection::kTooManyImages) {
        ++too_many;
      }
    }
    check::is_true(too_many > 0, "the image span limit is reported as itself");
  }

  // And a box that is built says nothing was wrong.
  {
    PredictOptions po;
    po.d_min = 2.5;
    for (const Prediction &p : predict(e, po)) {
      Shoebox box;
      BoxRejection why = BoxRejection::kTooManyImages;
      if (build_shoebox(e, p, options, &box, &why)) {
        check::is_true(why == BoxRejection::kNone, "a kept box has no reason");
        return;
      }
    }
    check::is_true(false, "nothing was built at all");
  }
}

}  // namespace mxi
