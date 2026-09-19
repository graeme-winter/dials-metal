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

#include "profile_model.h"
#include "shoebox.h"

namespace mxi {

//: A cube of (2n+1) points a side, spanning plus and minus `half_width` sigmas.
struct ProfileGrid {
  int n = 4;  //: so the side is 2n + 1
  double sigma_d = 0.0;
  double sigma_m = 0.0;
  double half_width = 3.0;  //: in sigmas
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
void add_to_grid(const Experiment &e, const Shoebox &box, const Vec3 &s1,
                 double phi_calculated, ProfileGrid *grid, int subdivisions = 5);

}  // namespace mxi
