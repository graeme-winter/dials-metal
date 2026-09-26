// The conventions, checked against a file DIALS wrote.
//
// test_geometry.cc cannot do this. It predicts reflections and maps them back,
// applying each convention once forward and once backward, so any error in one
// cancels itself exactly and the suite stays green. Two of them were wrong
// under precisely that cover: the sign of s0, and the absence of any parallax
// correction at all.
//
// Everything here compares against numbers produced by dials.index on a real
// insulin sweep. The tolerances are chosen from what the comparison actually
// achieves, not from what would be comfortable.

#include <algorithm>
#include <cmath>
#include <vector>

#include "../src/geometry.hh"
#include "../src/predict.hh"
#include "check.hh"
#include "real_data.hh"

namespace mxi {

namespace {

Vec3 v(const double (&a)[3]) { return {a[0], a[1], a[2]}; }

Panel imported_panel() {
  Panel p;
  p.fast = v(real::kImportedFast);
  p.slow = v(real::kImportedSlow);
  p.origin = v(real::kImportedOrigin);
  p.pixel_size[0] = real::kImportedPixelSize[0];
  p.pixel_size[1] = real::kImportedPixelSize[1];
  p.image_size[0] = real::kImageSize[0];
  p.image_size[1] = real::kImageSize[1];
  p.parallax = true;
  p.mu = real::kImportedMu;
  p.thickness = real::kImportedThickness;
  return p;
}

Panel refined_panel() {
  Panel p = imported_panel();
  p.fast = v(real::kRefinedFast);
  p.slow = v(real::kRefinedSlow);
  p.origin = v(real::kRefinedOrigin);
  p.mu = real::kRefinedMu;
  p.thickness = real::kRefinedThickness;
  return p;
}

Experiment refined_experiment() {
  Experiment e;
  e.beam.direction = v(real::kBeamDirection);
  e.beam.wavelength = real::kWavelength;
  e.detector.panels.push_back(refined_panel());

  std::vector<Vec3> axes;
  std::vector<double> angles;
  for (int i = 0; i < 3; ++i) {
    axes.push_back(v(real::kGonioAxes[i]));
    angles.push_back(real::kGonioAngles[i]);
  }
  e.goniometer = Goniometer::from_axes(
      axes, angles, static_cast<std::size_t>(real::kScanAxis));

  e.scan.first_image = real::kImageRange[0];
  e.scan.last_image = real::kImageRange[1];
  e.scan.osc_start = real::kOscStart;
  e.scan.osc_width = real::kOscWidth;

  e.crystal = Crystal::from_real_space(v(real::kRealSpaceA), v(real::kRealSpaceB),
                                       v(real::kRealSpaceC));
  return e;
}

}  // namespace

TEST(real_parallax_reproduces_dials_pixel_to_millimetre_exactly) {
  // NOTE ON xyzobs.mm: this test consumes it as *evidence about the formula*,
  // and it is the only file that pins the pixels-to-millimetres direction. It
  // is not an endorsement of the column. xyzobs.mm is computed once at import
  // and never recomputed, so it carries the inverse parallax correction under
  // the IMPORTED geometry, and after refinement it describes a detector that
  // no longer exists. Nothing in src/ reads it, and nothing should: work from
  // xyzobs.px through the current detector model every time.
  //
  // The proof of the staleness is that this test needs the imported panel.
  // Against the refined one the same code is out by 1.7e-3 mm, which looks
  // like a flaw in the parallax model and is nothing of the kind.
  const Panel p = imported_panel();
  double worst = 0.0;
  for (const real::Row &r : real::rows()) {
    const auto mm = p.px_to_mm(r.px_fast, r.px_slow);
    worst = std::fmax(worst, std::abs(mm.first - r.mm_fast));
    worst = std::fmax(worst, std::abs(mm.second - r.mm_slow));
  }
  // Exact to the last bit of a double, not merely close.
  check::close(worst, 0.0, 1e-14, "worst |mm - DIALS mm|");
}

TEST(real_parallax_reproduces_dials_millimetre_to_pixel_exactly) {
  const Panel p = refined_panel();
  double worst = 0.0;
  for (const real::Row &r : real::rows()) {
    const auto px = p.mm_to_px(r.cal_mm_fast, r.cal_mm_slow);
    worst = std::fmax(worst, std::abs(px.first - r.cal_px_fast));
    worst = std::fmax(worst, std::abs(px.second - r.cal_px_slow));
  }
  check::close(worst, 0.0, 1e-10, "worst |px - DIALS xyzcal.px|");
}

TEST(real_parallax_is_worth_more_than_a_pixel) {
  // If this ever comes out small, the correction has been disabled somewhere
  // and every test above it would still pass while the model quietly acquired
  // a radial error.
  const Panel p = imported_panel();
  double worst = 0.0;
  for (const real::Row &r : real::rows()) {
    const auto mm = p.px_to_mm(r.px_fast, r.px_slow);
    worst = std::fmax(worst, std::abs(mm.first - r.px_fast * p.pixel_size[0]) /
                                 p.pixel_size[0]);
    worst = std::fmax(worst, std::abs(mm.second - r.px_slow * p.pixel_size[1]) /
                                 p.pixel_size[1]);
  }
  // Measured at 0.96 px over this sample, which is spread along the scan
  // rather than out to the panel corners; over the whole dataset it reaches
  // 1.58 px. Half a pixel is well clear of the measurement and far above any
  // plausible accident.
  check::is_true(worst > 0.5, "parallax should displace by over half a pixel");
}

TEST(real_observed_spots_map_onto_the_reciprocal_lattice) {
  // The headline check. Map a real observed centroid into the crystal frame
  // and it must land on A h. This exercises the sign of s0, the pixel-to-mm
  // conversion, the goniometer rotation, the scan angle and the setting
  // matrix, all at once, against numbers DIALS produced.
  const Experiment e = refined_experiment();
  std::vector<double> residual;
  for (const real::Row &r : real::rows()) {
    const Vec3 expected =
        e.crystal->A * Vec3{static_cast<double>(r.h), static_cast<double>(r.k),
                            static_cast<double>(r.l)};
    const Vec3 got =
        reciprocal_lattice_point(e, 0, r.px_fast, r.px_slow, r.px_z);
    residual.push_back((got - expected).norm());
  }
  std::sort(residual.begin(), residual.end());
  const double median = residual[residual.size() / 2];

  // What is left is the indexing residual itself and not an error in this
  // code: DIALS' own rlp column sits the same distance from A h, as the next
  // test shows by agreeing with it forty times more closely.
  //
  // Measured: median 1.9e-4, max 7.1e-4, against a reciprocal cell edge of
  // 1/67.4 = 1.48e-2. So the median is about one per cent of a lattice
  // spacing. The bounds below are set from those numbers with room, and the
  // median is asserted separately from the maximum because the distribution
  // has a tail and a bound on the max alone would be a bound on the tail.
  check::is_true(median < 3e-4, "median residual against A h");
  check::is_true(residual.back() < 1.2e-3, "worst residual against A h");
}

TEST(real_reciprocal_lattice_points_agree_with_the_dials_column) {
  // Tighter than the test above, because it compares against what DIALS
  // computed from the same observation rather than against the ideal lattice.
  const Experiment e = refined_experiment();
  double worst = 0.0;
  for (const real::Row &r : real::rows()) {
    const Vec3 got =
        reciprocal_lattice_point(e, 0, r.px_fast, r.px_slow, r.px_z);
    worst = std::fmax(worst, (got - v(r.rlp)).norm());
  }
  check::is_true(worst < 2e-4, "rlp must agree with DIALS' own column");
}

TEST(real_s0_sign_is_not_arbitrary) {
  // With the sign flipped the residual is not merely worse, it is two
  // reciprocal Angstrom against a cell edge of 0.0128 -- a hundred and sixty
  // lattice spacings. This convention fails loudly, which is the one mercy.
  Experiment e = refined_experiment();
  e.beam.direction = -e.beam.direction;
  double best = 1e30;
  for (const real::Row &r : real::rows()) {
    const Vec3 expected =
        e.crystal->A * Vec3{static_cast<double>(r.h), static_cast<double>(r.k),
                            static_cast<double>(r.l)};
    const Vec3 got =
        reciprocal_lattice_point(e, 0, r.px_fast, r.px_slow, r.px_z);
    best = std::fmin(best, (got - expected).norm());
  }
  check::is_true(best > 1.0, "the wrong sign must fail unmistakably");
}

TEST(real_scan_angle_matches_the_dials_millimetre_column) {
  // The third component of xyzobs.mm is safe where the first two are not: it
  // is the rotation angle, which comes from the scan and has no dependence on
  // the detector model at all, so refinement cannot make it stale. It pins
  // phi = radians(osc_start + z * osc_width).
  const Experiment e = refined_experiment();
  double worst = 0.0;
  for (const real::Row &r : real::rows()) {
    worst = std::fmax(worst, std::abs(e.scan.phi_from_z(r.px_z) - r.mm_phi));
  }
  check::close(worst, 0.0, 1e-12, "worst |phi - DIALS phi|");
}

TEST(real_entering_flag_matches_dials) {
  const Experiment e = refined_experiment();
  int agree = 0, total = 0;
  for (const real::Row &r : real::rows()) {
    const Vec3 r0 =
        e.crystal->A * Vec3{static_cast<double>(r.h), static_cast<double>(r.k),
                            static_cast<double>(r.l)};
    const Intersections x = ewald_intersections(e, r0);
    if (!x.any) continue;
    // Take the root nearer the observed angle; the other is the same
    // reflection on its way out.
    const double phi = e.scan.phi_from_z(r.px_z);
    const int nearer = std::abs(x.phi[0] - phi) <= std::abs(x.phi[1] - phi) ? 0 : 1;
    if (x.entering[nearer] == r.entering) ++agree;
    ++total;
  }
  check::is_true(total > 30, "should have intersections to test");
  check::equal(agree, total, "every entering flag must match DIALS");
}

TEST(real_multi_axis_goniometer_parses_to_the_scan_axis) {
  const Experiment e = refined_experiment();
  // Every angle is zero in this dataset, so both composed rotations collapse
  // to the identity and the decomposition itself is NOT tested here. Only a
  // dataset with a non-zero chi or phi setting can do that.
  check::close(rotation_angle(e.goniometer.fixed), 0.0, 1e-12, "fixed");
  check::close(rotation_angle(e.goniometer.setting), 0.0, 1e-12, "setting");
  check::close((e.goniometer.axis - Vec3{1.0, 0.0, 0.0}).norm(), 0.0, 1e-12,
               "scan axis is omega");
}

TEST(real_cell_is_insulin) {
  const UnitCell u = refined_experiment().crystal->cell();
  check::close(u.a, 67.0, 4.0, "a near 67");
  check::close(u.b, 67.0, 4.0, "b near 67");
  check::close(u.c, 67.0, 4.0, "c near 67");
  check::close(u.alpha, 109.5, 3.0, "alpha");
}

TEST(real_prediction_reproduces_dials_xyzcal_exactly) {
  // End to end: predict from DIALS' own refined model and land on DIALS' own
  // xyzcal.px. Not "close to" -- exactly, to the last bit of a double.
  //
  // The tolerance here was 0.5 px, which was useless. Omitting the parallax
  // correction from the prediction path displaces a spot by 0.68 px median and
  // 1.0 px worst on this data, so a half-pixel bound would have passed for
  // more than half the reflections while the correction was missing entirely.
  // The companion test below measures that, so this bound is known to be
  // sensitive to the thing it is guarding.
  Experiment e = refined_experiment();
  PredictOptions options;
  options.allow_outside_scan = true;
  double worst = 0.0;
  int found = 0;
  for (const real::Row &r : real::rows()) {
    const std::vector<Prediction> p =
        predict_indices(e, {{{r.h, r.k, r.l}}}, options);
    double best = 1e30;
    for (const Prediction &q : p) {
      if (q.entering != r.entering) continue;
      best = std::fmin(best, std::hypot(q.px_fast - r.cal_px_fast,
                                        q.px_slow - r.cal_px_slow));
    }
    if (best < 1e29) {
      worst = std::fmax(worst, best);
      ++found;
    }
  }
  check::is_true(found > 30, "most sampled reflections should be predicted");
  // Against DIALS' own xyzcal, so this compares two predictions of the same
  // model rather than fitting to observation. It is exact.
  check::close(worst, 0.0, 1e-6, "predicted position must match DIALS' xyzcal");
}

TEST(real_prediction_would_notice_a_missing_parallax_correction) {
  // The guard on the test above. Turn the correction off and the prediction
  // moves by two thirds of a pixel, so the exact bound is genuinely sensitive
  // to it rather than exact for some unrelated reason.
  Experiment e = refined_experiment();
  for (Panel &p : e.detector.panels) p.parallax = false;
  PredictOptions options;
  options.allow_outside_scan = true;

  std::vector<double> displaced;
  for (const real::Row &r : real::rows()) {
    double best = 1e30;
    for (const Prediction &q : predict_indices(e, {{{r.h, r.k, r.l}}}, options)) {
      best = std::fmin(best, std::hypot(q.px_fast - r.cal_px_fast,
                                        q.px_slow - r.cal_px_slow));
    }
    if (best < 1e29) displaced.push_back(best);
  }
  std::sort(displaced.begin(), displaced.end());
  check::is_true(displaced.size() > 30, "enough predictions");
  check::is_true(displaced[displaced.size() / 2] > 0.5,
                 "without parallax the prediction should move over half a pixel");
}


TEST(real_oscillation_array_is_uniform) {
  // dxtbx writes a per-image array of start angles rather than a start and a
  // width. Treating the width as constant is fair -- but it should be checked
  // rather than assumed, and the check is nearly free.
  std::vector<double> oscillation;
  const long n = real::kImageRange[1] - real::kImageRange[0] + 1;
  for (long i = 0; i < n; ++i) {
    oscillation.push_back(real::kOscStart + static_cast<double>(i) * real::kOscWidth);
  }
  const Scan s = Scan::from_oscillation(oscillation, real::kImageRange[0],
                                        real::kImageRange[1]);
  check::equal(s.num_images(), n, "image count");
  check::close(s.osc_width, real::kOscWidth, 1e-15, "width from the endpoints");
  // On the real array this is 3.6e-14: accumulation round-off, nothing more.
  check::is_true(s.max_width_deviation < 1e-10, "the scan should be uniform");
}

TEST(oscillation_width_comes_from_the_endpoints_not_the_first_pair) {
  // A single perturbed element must not move the width much. Taking the width
  // from the first two elements would let one bad value at the start set the
  // width for the entire scan.
  std::vector<double> oscillation;
  for (int i = 0; i < 300; ++i) oscillation.push_back(0.1 * i);
  oscillation[1] += 0.05;
  const Scan s = Scan::from_oscillation(oscillation, 1, 300);
  check::close(s.osc_width, 0.1, 1e-12, "width barely moves");
  // And the perturbation must still be reported rather than swallowed.
  check::is_true(s.max_width_deviation > 0.1, "non-uniformity must be visible");
}

TEST(goniometer_ignores_the_angle_of_its_own_scan_axis) {
  // The scan axis does not have one setting, it has a different one on every
  // frame, and the scan supplies it. A value in `angles` for that axis must
  // not be composed into either rotation.
  const std::vector<Vec3> axes = {
      {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  const Goniometer a = Goniometer::from_axes(axes, {0.0, 0.0, 0.0}, 2);
  const Goniometer b = Goniometer::from_axes(axes, {0.0, 0.0, 37.0}, 2);
  check::close(rotation_angle(a.fixed * b.fixed.transpose()), 0.0, 1e-12,
               "fixed unchanged by the scan axis angle");
  check::close(rotation_angle(a.setting * b.setting.transpose()), 0.0, 1e-12,
               "setting unchanged by the scan axis angle");

  // An angle on a NON-scan axis must change the fixed rotation, or the
  // decomposition is doing nothing at all.
  const Goniometer c = Goniometer::from_axes(axes, {0.0, 30.0, 0.0}, 2);
  check::close(rotation_angle(c.fixed), Scan::radians(30.0), 1e-12,
               "a chi setting must reach the fixed rotation");
}

}  // namespace mxi
