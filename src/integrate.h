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
  //: The observed centre of mass, in pixels and images, over the foreground
  //: with the background taken off. dials.scale wants this; it is also the
  //: only thing in the output that says where the spot actually was rather
  //: than where the model put it.
  bool centroid_valid = false;
  double centroid_fast = 0.0, centroid_slow = 0.0, centroid_z = 0.0;
  //: Variance of that centre of mass: the second moment of the distribution
  //: over the square of the weight, which is the variance of a weighted mean
  //: and falls as a spot gets stronger.
  //:
  //: It does not reproduce DIALS' `xyzobs.px.variance`, and the comment in
  //: integrate.cc says what was measured. The centroid VALUE agrees to 0.07
  //: pixels, which is what `dials.scale` uses.
  double centroid_variance_fast = 0.0;
  double centroid_variance_slow = 0.0;
  double centroid_variance_z = 0.0;
  std::size_t n_foreground = 0;
  std::size_t n_background = 0;
  std::size_t n_valid = 0;
  //: Why it failed, when it did.
  bool background_failed = false;
  bool too_few_background = false;
};

//: The fraction of photons the sensor absorbs for a beam along `s1`.
double quantum_efficiency(const Panel &panel, const Vec3 &s1);

//: The resolution of a reflection, in Angstroms, as DIALS writes it into the
//: `d` column.
//:
//: From the STATIC cell and the Miller index, `1 / |A h|`, not from `s1` and
//: not from the scan-varying `A`. All three are nearly the same number and
//: only one of them is the column: against a DIALS `integrated.refl`, the
//: static cell agrees to 1.2e-15 for every reflection, `1/|s1 - s0|` to 2.6e-4
//: and the scan-varying `A` to 2.6e-4. Redundant information, and expected
//: anyway -- dials.scale stops with an unknown column error without it.
double resolution(const Crystal &crystal, int h, int k, int l);

//: The fraction of a reflection's rocking curve that the shoebox recorded.
//:
//: The box, not the scan. What was summed is what is in the box, so a
//: reflection wholly inside the scan is still not wholly measured: the box
//: spans n_sigma of the rocking curve and misses the tails beyond it. Against
//: DIALS that is the difference between 1.0000 and 0.99914, which is small and
//: is exactly the kind of small that a scale factor multiplies through a whole
//: dataset.
//:
//: Kabsch section 2.4: the Gaussian in phi has width sigma_m / |zeta|, and
//: this is the part of it between the two ends of the box, which are already
//: clipped to the scan.
double partiality(const Scan &scan, double phi, double zeta, double sigma_m,
                  std::int32_t z_first, std::int32_t z_last);

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
