// The Gaussian profile model: sigma_D, the beam divergence, and sigma_M, the
// reflecting range. Kabsch (2010a), Acta Cryst. D66, 133-144, section 3.1.
//
// These two numbers fix the integration region. The reflection mask is
//
//     |eps1| <= delta_D / 2,  |eps2| <= delta_D / 2,  |eps3| <= delta_M / 2
//
// in the per-reflection frame of section 2.3, with delta_D and delta_M taken
// as six to ten times sigma_D and sigma_M. Get them wrong and every shoebox is
// the wrong size, which costs intensity if too small and background if too
// large.
//
// STATUS
// ------
// `beam_divergence` implements Kabsch's estimator and is verified against
// synthetic spots whose angular spread is known in advance. On real insulin it
// gives 0.0304 degrees where dials.integrate records 0.0317, a difference of
// 4.2 per cent that is NOT yet explained -- see docs/integration.md. Do not
// treat the two as interchangeable until it is.

#pragma once

#include <cstddef>
#include <vector>

#include "geometry.h"
#include "shoebox.h"

namespace mxi {

struct ProfileModel {
  //: Standard deviation of the beam divergence, in degrees.
  double sigma_d = 0.0;
  //: Standard deviation of the reflecting range, in degrees.
  double sigma_m = 0.0;
  //: How many standard deviations the integration region spans.
  double n_sigma = 3.0;
  //: Reflections that contributed.
  std::size_t n_used = 0;
};

//: One reflection's contribution: the variance of the angles between its
//: pixels' diffracted-beam directions and its own diffracted beam, weighted by
//: the counts.
//:
//: Returns false if the box carries too little signal to have a variance --
//: one count cannot have a spread, and dividing by zero to say so would put a
//: NaN into the mean and lose every other reflection with it.
bool spot_angular_variance(const Experiment &e, const Shoebox &box, const Vec3 &s1,
                           double *variance);

//: One contribution to the reflecting-range estimate: a reflection seen on one
//: image, the angular gap between its Bragg maximum and that image's centre,
//: and its zeta.
//:
//: One per IMAGE, not one per reflection. With one per reflection the gap is
//: bounded by half an oscillation width by construction, the likelihood is
//: maximised by driving sigma to zero, and the estimate is meaningless. A
//: mosaic crystal puts a spot on images well away from its Bragg angle, and
//: that spread is the entire signal.
struct RangeSample {
  //: Radians between the Bragg maximum and the centre of the image.
  double delta = 0.0;
  //: |m2 . e1|, which corrects for the path length through the Ewald sphere.
  double zeta = 0.0;
};

//: The fraction of a reflection's intensity recorded on an image whose centre
//: is `delta` away, for a reflecting range of `sigma` and the given zeta.
//: Kabsch section 2.4, as a density in delta.
double recorded_fraction(double delta, double zeta, double sigma,
                         double oscillation);

//: Gather one sample per image a reflection was seen on, from the frames of
//: its shoebox that carry foreground counts.
std::vector<RangeSample> range_samples(const Experiment &e, const Shoebox &box,
                                       double phi_calculated, double zeta);

//: zeta = m2 . e1, with e1 the axis a reciprocal lattice point would cross the
//: Ewald sphere about by the shortest route.
//:
//: e1 is the CROSS product of the diffracted and incident beams, s1 x s0,
//: normalised -- perpendicular to both, as Kabsch section 2.3 defines it after
//: Schutt & Winkler. Using the difference s1 - s0 instead, which is the
//: reciprocal lattice vector and points somewhere else entirely, put sigma_M
//: out by fifty per cent.
double compute_zeta(const Experiment &e, const Vec3 &s1);

//: sigma_M by maximum likelihood over those samples. Kabsch section 3.1.
//:
//: Samples with |zeta| below `min_zeta` are dropped: their reflecting range is
//: sigma_M / |zeta|, which diverges, and they carry no information about
//: sigma_M while dominating the likelihood.
double reflecting_range(const std::vector<RangeSample> &samples, double oscillation,
                        double min_zeta = 0.05);

//: sigma_D as the root mean of those variances, over the reflections given.
//:
//: Kabsch section 3.1: "determine the centroid and variance s^2 of the
//: intensity-weighted diffracted-beam directions", then
//: sigma_D^2 = (1/n) sum s_j^2.
double beam_divergence(const Experiment &e, const std::vector<Shoebox> &boxes,
                       const std::vector<Vec3> &s1, std::size_t *n_used);

}  // namespace mxi
