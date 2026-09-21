// Summation integration: Leslie (1999), sections 5.2 and 5.3.
//
// The intensity is the background-subtracted sum over the foreground,
//
//     I = sum_{i in foreground} (c_i - b)
//
// and the variance is the part that matters, because it is what says a weak
// reflection is weak:
//
//     var(I) = G [ I + I_bg + (m/n) I_bg ]
//
// with m foreground pixels, n background pixels, I_bg the background summed
// over the foreground, and G the gain. The first term is the Poisson noise of
// the signal, the second the Poisson noise of the background under it, and the
// third the uncertainty in the background estimate itself -- which is why the
// background region wants to be large. Leslie's point is that for a weak
// reflection the second and third dominate entirely: the error is set by the
// background, not by the signal.
//
// WHAT IS NOT APPLIED HERE
// ------------------------
// DIALS stores `intensity.sum.value` uncorrected and keeps `lp` and `qe`
// beside it for the scaler to apply later. This does the same, so that the two
// can be compared directly, and computes `qe` for the table.
//
// The quantum efficiency is the fraction of photons the sensor stops rather
// than passes through:
//
//     qe = 1 - exp(-mu t / cos(theta))
//
// with theta the angle between the diffracted beam and the detector normal.
// Checked against the `qe` column of a DIALS integrated.refl: identical for
// 100 per cent of reflections, to 2e-16.

#pragma once

#include <cstddef>
#include <vector>

#include "background.h"
#include "geometry.h"
#include "shoebox.h"

namespace mxi {

struct IntegrateOptions {
  //: Detector gain, counts per photon. One for a photon counting detector,
  //: which is what this is for.
  double gain = 1.0;
  BackgroundOptions background;
  //: Reject a shoebox with fewer than this many background pixels; the
  //: background estimate and its contribution to the variance both need them.
  std::size_t min_background = 10;
};

struct IntegratedReflection {
  bool valid = false;
  //: Background-subtracted sum over the foreground, uncorrected.
  double intensity = 0.0;
  double variance = 0.0;
  //: The fitted background, per pixel, and summed over the foreground.
  double background_mean = 0.0;
  double background_sum = 0.0;
  double background_sum_variance = 0.0;
  std::size_t n_foreground = 0;
  std::size_t n_background = 0;
  std::size_t n_valid = 0;
  //: Why it failed, when it did.
  bool background_failed = false;
  bool too_few_background = false;
};

//: The fraction of photons the sensor absorbs for a beam along `s1`.
double quantum_efficiency(const Panel &panel, const Vec3 &s1);

//: The Lorentz-polarization factor, as DIALS writes it into the `lp` column.
//:
//:     L  = |s1 . (m2 x s0)| / (|s1| |s0|)
//:     P  = (1 - p) + (2p - 1) (u . n)^2 + p (u . s0hat)^2
//:     lp = L / P
//:
//: with `u` the unit diffracted beam, `n` the polarization normal and `p` the
//: polarization fraction. Recovered from a DIALS `integrated.refl` rather than
//: recalled: the Lorentz part was identifiable by having the least scatter
//: against the column, and the three coefficients of P were then solved for by
//: least squares and came back as 0.001000, 0.998000 and 0.999000 with a
//: residual of 2e-16 -- which is (1-p), (2p-1) and p for the p = 0.999 in the
//: file, not a curve fit.
//:
//: At p = 0.5 it collapses to the unpolarized (1 + cos^2 2theta)/2, which is
//: the check that the form is right rather than merely fitted.
//:
//: Stored rather than applied, as DIALS does, so the scaler applies it later.
double lorentz_polarization(const Beam &beam, const Goniometer &goniometer,
                            const Vec3 &s1);

//: Integrate one shoebox whose `data` holds counts and whose `mask` marks
//: foreground and background. Fills the shoebox's `background` array with the
//: fitted value, so that a saved shoebox carries what was subtracted.
IntegratedReflection integrate_shoebox(Shoebox *box,
                                       const IntegrateOptions &options = {});

}  // namespace mxi
