#pragma once

// The scaling model: g = C(r) . exp(B(t) / 2d^2) . S(s0, s1), the physical
// model of Beilsten-Edmands et al. (2020), parameterised with cubic B-splines
// for the scale C and relative B factor B, where DIALS averages the nearest
// parameters under a Gaussian, and the paper's spherical harmonics (eqn 7) for
// the absorption surface S.

#include <cstddef>
#include <vector>

#include "linalg.hh"

namespace mxi {

//: Real spherical harmonics of degree 1 to lmax at a unit vector, orthonormal
//: on the sphere, in the order l = 1..lmax and m = -l..l: lmax (lmax + 2) of
//: them. Without the Condon-Shortley sign, which would only flip the sign of a
//: fitted coefficient: degree 1 is sqrt(3 / 4 pi) (y, z, x) for m = -1, 0, 1.
void real_spherical_harmonics(int lmax, const Vec3 &unit,
                              std::vector<double> *out);
std::size_t harmonic_count(int lmax);

//: What the model needs of one observation.
struct ScaleObservation {
  double rotation = 0.0; //: position in the scan's rotation, 0 to 1, for C
  double time = 0.0;     //: position in the scan's time, 0 to 1, for B
  double inv_2d2 = 0.0;  //: 1 / (2 d^2)
  //: [Y(s1) + Y(s0)] / 2 in the crystal frame, harmonic_count(lmax) of them;
  //: empty without absorption.
  std::vector<double> absorption;
};

struct ScaleModelShape {
  std::size_t scale_points = 1;
  std::size_t decay_points = 0; //: none: no decay term
  int lmax = 0;                 //: 0: no absorption term
};

//: The number of control points DIALS' physical model uses for a sweep, in
//: effect: 6 and 5 for scale and decay on 30 degrees; for a sweep of 90 degrees
//: or more, one per 15 and 20 degrees and two more. Absorption from 60
//: degrees, as dials.scale's automatic choice.
ScaleModelShape default_shape(double degrees);

class ScaleModel {
public:
  explicit ScaleModel(const ScaleModelShape &shape);

  const ScaleModelShape &shape() const { return shape_; }
  std::size_t size() const { return parameters.size(); }
  //: C control points, then B, then P_lm; C starts at 1, the rest at 0.
  std::vector<double> parameters;
  std::size_t first_decay() const { return shape_.scale_points; }
  std::size_t first_absorption() const {
    return shape_.scale_points + shape_.decay_points;
  }

  //: g for an observation, and if asked its derivatives: pairs of parameter
  //: index and dg/dp, repeated indices to be summed.
  double inverse_scale(
      const ScaleObservation &o,
      std::vector<std::pair<std::size_t, double>> *gradient = nullptr) const;

  //: The scale's control points divided by their mean. g and the merged
  //: intensities trade a common factor freely, so this changes no fit; it keeps
  //: the degenerate direction from wandering.
  void normalise();

private:
  ScaleModelShape shape_;
};

} // namespace mxi
