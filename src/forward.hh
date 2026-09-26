// The spot the model predicts, rendered onto the pixels it would be recorded
// on, so that it can be compared with the spot that was.
//
// Every measurement before this one compared a number taken from the data with
// a number taken from the model, and the two were computed differently: the
// data through a pixel grid that truncates and quantises it, the model in
// closed form. At widths below a pixel that difference is most of what is
// being measured.
//
// Rendering removes it. The model is integrated over the same pixels, over the
// same images, inside the same mask, and its centre of mass and moments are
// taken the same way. Whatever the grid does to the data it does to the model,
// and what is left over is the model being wrong.
//
// The physics put in, beyond the Gaussian:
//
//   * a photon is absorbed at a depth in the sensor, exponentially distributed
//     and cut off at the back, and is recorded where it stopped rather than
//     where it entered. That displaces it along the ray's projection on the
//     face -- radially -- by an amount that varies, so it both shifts and
//     smears. The shift is what the parallax correction removes; the smear is
//     what it cannot.

#pragma once

#include <cstddef>
#include <vector>

#include "profile_model.hh"
#include "shoebox.hh"

namespace mxi {

struct ForwardOptions {
  double sigma_d = 0.0;  //: degrees
  double sigma_m = 0.0;  //: degrees
  //: Subdivisions of each pixel per axis when integrating the model over it.
  int subdivisions = 5;
  //: Samples through the sensor depth. One puts every photon at the mean
  //: depth, which is the parallax correction and no smear.
  int depth_samples = 8;
  //: Whether to model the sensor at all.
  bool sensor = true;
};

//: Counts the model predicts in each voxel of `box`, normalised to sum to one
//: over the voxels the mask calls valid. Empty if the geometry fails.
std::vector<double> render_shoebox(const Experiment &e, const Shoebox &box,
                                   const Vec3 &s1, double phi_calculated,
                                   const ForwardOptions &options);

//: Centre of mass and second moments of a set of counts over a shoebox, in
//: PIXELS and IMAGES -- the units the detector records in, so that a model and
//: an observation can be compared without either being transformed.
struct Moments {
  bool valid = false;
  double com_fast = 0.0, com_slow = 0.0, com_z = 0.0;
  double width_fast = 0.0, width_slow = 0.0, width_z = 0.0;
  double total = 0.0;
};
Moments moments_of(const Shoebox &box, const std::vector<double> &counts,
                   bool subtract_background);

}  // namespace mxi
