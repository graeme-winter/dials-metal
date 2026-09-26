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

#include "predict.hh"
#include "profile_model.hh"
#include "shoebox.hh"

namespace mxi {

//: What "within n sigma" means for a three-dimensional Gaussian.
enum class RegionShape {
  //: |eps1| <= n sigma_D and |eps2| <= n sigma_D and |eps3| <= n sigma_M.
  //: What Kabsch section 3.1 writes for the reflection mask, and what DIALS
  //: uses, so it is the default.
  kBox,
  //: (eps1/sigma_D)^2 + (eps2/sigma_D)^2 + (eps3/sigma_M)^2 <= n^2. The
  //: surface on which the model's density is actually constant.
  kEllipsoid,
};

struct MaskOptions {
  //: How many sigmas the FOREGROUND region spans, in each direction.
  double n_sigma = 3.0;
  //: How much wider the box is than the foreground region, on the detector.
  //:
  //: A box that is only the foreground has no background in it, and the
  //: background is where the background estimate comes from. DIALS' boxes are
  //: measurement boxes in Leslie's sense: foreground plus a rim, and the rim is
  //: most of the volume. Measured on 4018 reflections of insulin against a
  //: DIALS `integrated.refl`: foreground 462 voxels, background 3094, and
  //: their sum is the whole box for every row.
  //:
  //: The value here reproduces DIALS' box widths on that data. It is
  //: empirical: DIALS' own box is not simply n_sigma times anything this can
  //: see, and the factor is recorded as measured rather than derived.
  double box_scale = 1.9;
  //: The box is a region in which each coordinate is separately within n
  //: sigma; the ellipsoid is the one on which the Gaussian is constant. They
  //: are not the same set, and the difference is not small: the ellipsoid is
  //: pi/6 of the box, a little over half.
  //:
  //: At the corner of the box all three coordinates are at n sigma at once,
  //: so a three-dimensional Gaussian is at exp(-3 n^2 / 2) there -- 1.4e-6 of
  //: its peak at n = 3. That corner is pure background, and there are eight of
  //: them. The box is the larger region and captures slightly more of the
  //: density, 0.9919 against 0.9707, but it buys that with volume that holds
  //: essentially nothing.
  RegionShape shape = RegionShape::kBox;
  double sigma_d = 0.0; //: degrees
  double sigma_m = 0.0; //: degrees
  //: Reflections whose zeta is smaller than this are skipped: their region in
  //: rotation is sigma_M / |zeta|, which diverges, and the box would swallow
  //: the whole scan.
  double min_zeta = 0.05;
  //: A box wider than this many images is refused rather than allocated. A
  //: mistake in zeta or sigma_M turns into gigabytes of empty shoebox
  //: otherwise, and the first sign of it should not be the machine swapping.
  std::int32_t max_images = 200;
};

//: Why a reflection got no box. Counted rather than guessed: on a ten rotation
//: sweep half the predictions produced no shoebox and there was no way to tell
//: which of these it was.
enum class BoxRejection {
  kNone,
  kNoPanel,        //: the prediction names a panel that is not there
  kNoFrame,        //: s1 and s0 are parallel, so there is no Kabsch frame
  kSmallZeta,      //: too near the rotation axis; the box would span the scan
  kNoIntersection, //: a corner of the region misses the detector plane
  kTooManyImages,  //: wider in rotation than any reflection should be
  kOffDetector,    //: clipped away entirely by the panel or the scan
};

const char *describe(BoxRejection why);

//: The bounding box of a reflection's integration region, half open, clipped
//: to the panel and the scan. Returns false if it does not intersect them.
bool integration_bbox(const Experiment &e, const Prediction &p,
                      const MaskOptions &options, std::int32_t bbox[6],
                      BoxRejection *why = nullptr);

//: A shoebox for that box, with the mask marking which voxels are inside the
//: region: Valid | Foreground inside, Valid | Background outside. Values and
//: background are zero.
bool build_shoebox(const Experiment &e, const Prediction &p,
                   const MaskOptions &options, Shoebox *box,
                   BoxRejection *why = nullptr);

} // namespace mxi
