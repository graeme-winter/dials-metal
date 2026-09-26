// Spots in Kabsch space: the density of a reflection on a grid in
// (eps1, eps2, eps3), where every reflection is supposed to look the same.
//
// This is a diagnostic, not the profile fitting. What it is for is looking at
// the thing the Gaussian is supposed to describe, next to an average of the
// spots around it, and seeing whether either resembles the other.
//
// Two departures from Kabsch section 3.3, both because this measures the data
// rather than applies a model:
//
//   * a pixel is subdivided in the detector plane and its counts shared
//     between the subdivisions, as Kabsch does, because the data are badly
//     undersampled -- a shoebox is a handful of pixels across and assigning
//     each whole pixel to one grid point aliases the result into a staircase;
//
//   * along e3, an image's counts are shared between grid planes by the plain
//     geometric overlap of the image's angular range with the plane's, NOT by
//     the Gaussian weights of Kabsch's f_3j. Using the model to place the
//     counts would beg the question this is asked to answer.

#pragma once

#include <cstddef>
#include <vector>

#include "profile_model.hh"
#include "shoebox.hh"

namespace mxi {

//: A cube of (2n+1) points a side, spanning plus and minus `half_width` sigmas.
struct ProfileGrid {
  int n = 4; //: so the side is 2n + 1
  double sigma_d = 0.0;
  double sigma_m = 0.0;
  double half_width = 3.0; //: in sigmas
  std::vector<double> value;
  //: How many spots were added, and how much of their signal landed inside.
  std::size_t n_spots = 0;
  double counts_added = 0.0;
  double counts_outside = 0.0;

  int side() const { return 2 * n + 1; }
  std::size_t size() const {
    return static_cast<std::size_t>(side()) * static_cast<std::size_t>(side()) *
           static_cast<std::size_t>(side());
  }
  std::size_t at(int i1, int i2, int i3) const {
    return (static_cast<std::size_t>(i3) * static_cast<std::size_t>(side()) +
            static_cast<std::size_t>(i2)) *
               static_cast<std::size_t>(side()) +
           static_cast<std::size_t>(i1);
  }
  void reset();
  //: Divide through by the total, so grids from different numbers of spots can
  //: be looked at side by side.
  void normalise();
};

ProfileGrid make_grid(int n, double sigma_d, double sigma_m, double half_width);

//: Add one shoebox, transformed into the frame of its own reflection.
//:
//: `subdivisions` splits each pixel that many ways along each detector axis;
//: five is what Kabsch uses. One means no subdivision, which is what makes the
//: undersampling visible.
//: `recentre` shifts each spot onto its own centroid in the frame before it is
//: added, instead of onto its predicted position.
//:
//: This matters more than it sounds. The prediction is off by a mean of 0.0135
//: degrees on the detector, which is 0.44 sigma_D and is the same size as the
//: measured width of a spot. Added about the prediction, an aggregate carries
//: that error convolved into it and is blurred by roughly a factor of the
//: square root of two -- so the average looks wider, rounder and less
//: structured than the spots actually are.
//:
//: Which is right depends on the question. About the prediction is what
//: integration sees, so it is the right centre for asking whether a mask
//: captures the signal. About the spot is the right centre for asking what
//: shape a spot is, which is what an aggregate profile is for.
void add_to_grid(const Experiment &e, const Shoebox &box, const Vec3 &s1,
                 double phi_calculated, ProfileGrid *grid, int subdivisions = 5,
                 bool recentre = false);

//: Counts-weighted second moments of one spot's density in the Kabsch frame,
//: about its own centroid there.
//:
//: About the centroid, not about zero: a spot whose predicted position is a
//: little off would otherwise report that error as width, and the question
//: here is the SHAPE.
struct SpotMoments {
  bool valid = false;
  //: Where the spot's centre of mass sits in the frame, relative to where the
  //: model predicts it, in degrees. Not small: a mean of 0.0135 on the
  //: detector, which is 0.44 sigma_D.
  double centre1 = 0.0, centre2 = 0.0, centre3 = 0.0;
  double width1 = 0.0, width2 = 0.0, width3 = 0.0; //: degrees
  //: The same spread resolved along two axes FIXED IN THE LABORATORY rather
  //: than in the reflection's own frame: along the rotation axis, and
  //: perpendicular to both it and the beam.
  //:
  //: The distinction decides what the anisotropy is. A property of the
  //: detector or of the scattering geometry is radial and belongs in
  //: width1 and width2; a property of the beam -- and a synchrotron beam is
  //: routinely wider horizontally than vertically -- is fixed in the
  //: laboratory and belongs here. One averages away in the other.
  double width_along_axis = 0.0, width_across_axis = 0.0;
  double counts = 0.0;
  //: Where it sits, for binning: distance from the beam centre in millimetres,
  //: and the angle between the incident ray and the detector normal.
  double radius_mm = 0.0;
  double obliquity = 0.0; //: radians
  double path_mm = 0.0;   //: crystal to pixel
};
SpotMoments spot_moments(const Experiment &e, const Shoebox &box,
                         const Vec3 &s1, double phi_calculated);

//: The width a sensor of this thickness and absorption adds to a spot, in
//: degrees, for a ray striking at `obliquity` and recorded `path` away.
//:
//: A photon entering at an angle is absorbed at a random depth, exponentially
//: distributed and cut off at the back of the sensor, and the charge is
//: recorded where it lands rather than where the ray entered. Projected onto
//: the face that is a smear of `tan(obliquity) * sigma_depth` along the plane
//: containing the ray and the normal -- the RADIAL direction -- and none at
//: all across it.
//:
//: So it belongs entirely to eps2 and not at all to eps1, which is what makes
//: it testable: the two should differ by this much and by nothing else.
double sensor_depth_width(double mu, double thickness, double obliquity,
                          double path);

//: The width a crystal of this extent along the beam adds to a spot, in
//: degrees, for a ray leaving at `obliquity` and recorded `path` away.
//:
//: Diffraction from the front and the back of the crystal starts from points
//: separated along the beam, and at a scattering angle those two rays land
//: `extent * tan(2 theta)` apart on the detector: a hundred micrometres at
//: thirty degrees is fifty-eight. Purely radial, like the sensor.
//:
//: Which is the difficulty. Seen from the crystal both smears go as
//: sin(2 theta) cos(2 theta) / distance, so they have the SAME dependence on
//: angle and this data cannot tell them apart -- only their sum can be
//: measured, and only if it follows that shape at all.
//:
//: A uniform slab of thickness T illuminated throughout has a standard
//: deviation of T / sqrt(12), which is what `extent` means here.
double source_extent_width(double extent, double obliquity, double path);

} // namespace mxi
