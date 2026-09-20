// Shoeboxes built from the profile model, for looking at.
//
// Given sigma_D and sigma_M, the integration region of a reflection is
//
//     |eps1| <= n sigma_D,  |eps2| <= n sigma_D,  |eps3| <= n sigma_M
//
// in the frame of section 2.3. That is a region in Kabsch space; a shoebox is
// a box in pixels and images. This turns one into the other: the bounding box
// that encloses the region, and a mask saying which voxels inside the box are
// actually inside the region.
//
// The pixel values are left at zero. Nothing here has read an image, and a
// shoebox with invented counts in it would be worse than an empty one -- the
// point is to see WHERE the model says the signal is, against the real image
// in a viewer that draws both.

#pragma once

#include <cstddef>
#include <vector>

#include "predict.h"
#include "profile_model.h"
#include "shoebox.h"

namespace mxi {

struct MaskOptions {
  //: How many sigmas the region spans, in each direction.
  double n_sigma = 3.0;
  double sigma_d = 0.0;  //: degrees
  double sigma_m = 0.0;  //: degrees
  //: Reflections whose zeta is smaller than this are skipped: their region in
  //: rotation is sigma_M / |zeta|, which diverges, and the box would swallow
  //: the whole scan.
  double min_zeta = 0.05;
  //: A box wider than this many images is refused rather than allocated. A
  //: mistake in zeta or sigma_M turns into gigabytes of empty shoebox
  //: otherwise, and the first sign of it should not be the machine swapping.
  std::int32_t max_images = 200;
};

//: The bounding box of a reflection's integration region, half open, clipped
//: to the panel and the scan. Returns false if it does not intersect them.
bool integration_bbox(const Experiment &e, const Prediction &p,
                      const MaskOptions &options, std::int32_t bbox[6]);

//: A shoebox for that box, with the mask marking which voxels are inside the
//: region: Valid | Foreground inside, Valid | Background outside. Values and
//: background are zero.
bool build_shoebox(const Experiment &e, const Prediction &p,
                   const MaskOptions &options, Shoebox *box);

}  // namespace mxi
