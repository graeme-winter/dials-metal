// Analytical first derivatives of the prediction formula.
//
// Following Waterman et al. (2016), Acta Cryst. D72, 558-575, Appendix A. The
// equation numbers in the comments are that paper's.
//
// WHY, given that the numerical version works: on a device without double
// precision a finite difference is a difference of two nearly equal residuals,
// and however carefully the step is chosen some of the significance is spent
// on the cancellation. An analytical derivative spends none. Measured earlier,
// float32 leaves about four digits in a numerical derivative at the step this
// code uses; an analytical one would leave seven.
//
// And the cost: a wrong analytical derivative does not crash. It converges
// smoothly to the wrong answer and reports a small residual doing it. So the
// numerical version stays, it is the oracle, and `tests/test_derivatives.cc`
// compares the two element by element on real geometry. Nothing here should be
// trusted further than that test.
//
// The chain, for a parameter p of the crystal:
//
//   dphi/dp  = -(R_phi dr0/dp . s1) / ((e x r_phi) . s0)        eqn (40)
//   dr_phi/dp = (e x r_phi) dphi/dp + R_phi dr0/dp              eqn (46)
//   ds1/dp   = dr_phi/dp                    (s0 does not depend on p)
//   dv/dp    = D ds1/dp                     (detector fixed)     eqn (45)
//   dX/dp    = (w du/dp - u dw/dp) / w^2                         eqn (43)
//
// with v = (u, v, w) = D s1 and X = u/w, Y = v/w in millimetres on the panel.
//
// The denominator of eqn (40) is the volume of the parallelepiped formed by
// the rotation axis, the reciprocal lattice vector and the beam. It vanishes
// for reflections near the rotation axis, where phi is genuinely ill
// determined -- the same reflections that carry large Lorentz factors. DIALS
// discards any below 0.05 by default. `volume` is returned so a caller can do
// the same rather than dividing by something close to zero.

#pragma once

#include <array>
#include <vector>

#include "geometry.hh"

namespace mxi {

//: Derivative of the predicted centroid with respect to one parameter, in
//: millimetres on the panel and radians of rotation.
struct CentroidDerivative {
  double dX = 0.0;
  double dY = 0.0;
  double dphi = 0.0;
};

//: Everything about one predicted reflection that the derivatives need, so it
//: is computed once rather than nine times.
struct PredictionState {
  bool valid = false;
  Vec3 r0;    // reciprocal lattice point in the crystal frame
  Vec3 r_phi; // rotated into the laboratory frame
  Vec3 s1;    // diffracted beam
  Vec3 v;     // D s1, homogeneous panel coordinates
  Mat3 D;     // inverse of (fast | slow | origin) as columns
  Mat3 R_phi; // the goniometer rotation at the diffracting angle
  Vec3 axis;  // the rotation axis in the laboratory frame
  double phi = 0.0;
  //: (e x r_phi) . s0, the denominator of eqn (40). Small means phi is ill
  //: determined; DIALS discards below 0.05.
  double volume = 0.0;
};

// Assemble the state for a reflection observed at scan position `z`, choosing
// the Ewald root nearer the observation as the refinement target does.
PredictionState prediction_state(const Experiment &e, std::size_t panel, int h,
                                 int k, int l, double z);

// Derivatives with respect to the nine elements of the setting matrix A, in
// row-major order. For a scan-varying crystal these are the derivatives with
// respect to A at this reflection's own scan position; multiplying by the
// B-spline weight of a control point gives the derivative with respect to that
// control point, which is where the banding comes from.
std::array<CentroidDerivative, 9> crystal_derivatives(const PredictionState &s,
                                                      int h, int k, int l);

// The B-spline weights of the control points at scan position `z`, and the
// index of the first of the four. Everything outside those four is exactly
// zero, which is the fact a device implementation should exploit.
struct SplineWeights {
  //: Control points touched, and their weights. Indices repeat near the ends
  //: of the scan, where the clamping duplicates the outermost control point;
  //: a caller accumulates, so repeats are correct rather than a special case.
  std::size_t index[4] = {0, 0, 0, 0};
  double weight[4] = {0.0, 0.0, 0.0, 0.0};
  std::size_t count = 0;
};
SplineWeights spline_weights(const Experiment &e, double z);

// --------------------------------------------------------------------------
// Detector and beam
// --------------------------------------------------------------------------
//
// The perturbations live here, next to their derivatives, and refinement uses
// these same functions. A derivative and the thing it differentiates drifting
// apart is a silent failure: the refinement would take confident steps in a
// direction that does not correspond to how it then moves the model.

//: Six parameters: three translations along the laboratory axes, then three
//: rotations about the laboratory axes THROUGH THE PANEL CENTRE. Rotating
//: about the laboratory origin instead would be mostly a translation for a
//: panel two hundred millimetres away, and the two groups would correlate
//: badly enough to make the normal matrix near singular.
Panel perturb_panel(const Panel &p, const double shift[6]);

//: Two parameters: tilts about the two directions perpendicular to the beam.
//: A third would be a rotation about the beam itself and would do nothing.
Beam perturb_beam(const Beam &b, const double shift[2]);

//: Derivatives with respect to the six detector parameters. The rotation angle
//: does not depend on them at all -- neither r0 nor s0 does -- so `dphi` is
//: exactly zero, which is worth knowing when building a normal matrix.
std::array<CentroidDerivative, 6> detector_derivatives(const PredictionState &s,
                                                       const Panel &p);

//: Derivatives with respect to the two beam parameters.
std::array<CentroidDerivative, 2> beam_derivatives(const PredictionState &s,
                                                   const Beam &b);

} // namespace mxi
