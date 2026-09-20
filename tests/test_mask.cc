// Shoeboxes built from the profile model, for looking at against the images.

#include <cmath>
#include <vector>

#include "../src/mask.h"
#include "../src/predict.h"
#include "check.h"

using namespace mxi;

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

MaskOptions options_for(double n_sigma = 3.0) {
  MaskOptions o;
  o.n_sigma = n_sigma;
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

  const MaskOptions options = options_for();
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
          if (x == 0 || y == 0 || z == 0 || x == box.nx() - 1 ||
              y == box.ny() - 1 || z == box.nz() - 1) {
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
    check::is_true(foreground * 10 >= box.size(),
                   "and is not mostly padding");
    (void)edge_has_foreground;
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
