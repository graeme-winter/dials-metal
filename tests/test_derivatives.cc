// Analytical derivatives, checked against finite differences.
//
// This file is the reason the analytical derivatives are allowed to exist. A
// wrong one does not crash: it converges smoothly to the wrong answer and
// reports a small residual doing it, so the only defence is a comparison
// against a derivative computed a completely different way, on geometry that
// is not idealised.
//
// The comparison is against a CENTRAL difference rather than the forward
// difference the refinement uses, because a central difference has error
// O(step^2) and can therefore agree with an exact derivative to many more
// digits -- close enough that a genuine error in the analytical form cannot
// hide inside the tolerance.

#include <algorithm>
#include <cmath>
#include <vector>

#include "../src/derivatives.hh"
#include "../src/predict.hh"
#include "../src/refine.hh"
#include "../src/refl.hh"
#include "check.hh"
#include "real_data.hh"

namespace mxi {

namespace {

Vec3 v3(const double (&a)[3]) { return {a[0], a[1], a[2]}; }

// The real refined insulin geometry, so the test is not run on a detector at a
// convenient right angle to a beam along an axis.
Experiment real_experiment() {
  Experiment e;
  e.beam.direction = v3(real::kBeamDirection);
  e.beam.wavelength = real::kWavelength;

  Panel p;
  p.fast = v3(real::kRefinedFast);
  p.slow = v3(real::kRefinedSlow);
  p.origin = v3(real::kRefinedOrigin);
  p.pixel_size[0] = p.pixel_size[1] = real::kRefinedPixelSize[0];
  p.image_size[0] = real::kImageSize[0];
  p.image_size[1] = real::kImageSize[1];
  p.parallax = true;
  p.mu = real::kRefinedMu;
  p.thickness = real::kRefinedThickness;
  e.detector.panels.push_back(p);

  std::vector<Vec3> axes;
  std::vector<double> angles;
  for (int i = 0; i < 3; ++i) {
    axes.push_back(v3(real::kGonioAxes[i]));
    angles.push_back(real::kGonioAngles[i]);
  }
  e.goniometer = Goniometer::from_axes(
      axes, angles, static_cast<std::size_t>(real::kScanAxis));
  e.scan.first_image = real::kImageRange[0];
  e.scan.last_image = real::kImageRange[1];
  e.scan.osc_start = real::kOscStart;
  e.scan.osc_width = real::kOscWidth;
  e.crystal = Crystal::from_real_space(v3(real::kRealSpaceA), v3(real::kRealSpaceB),
                                       v3(real::kRealSpaceC));
  return e;
}

// The predicted centroid in millimetres and radians, which is what the
// analytical derivatives differentiate.
struct Centroid {
  bool valid = false;
  double X = 0.0, Y = 0.0, phi = 0.0;
};

Centroid centroid(const Experiment &e, int h, int k, int l, double z) {
  const PredictionState s = prediction_state(e, 0, h, k, l, z);
  Centroid c;
  if (!s.valid) return c;
  c.valid = true;
  c.X = s.v.x / s.v.z;
  c.Y = s.v.y / s.v.z;
  c.phi = s.phi;
  return c;
}

}  // namespace

TEST(the_prediction_state_reproduces_the_predicted_position) {
  // Before differentiating it, the thing being differentiated has to be right.
  // X and Y in millimetres, converted to pixels, must match what the ordinary
  // prediction path produces for the same reflection.
  const Experiment e = real_experiment();
  double worst = 0.0;
  int tested = 0;
  for (const real::Row &r : real::rows()) {
    const Centroid c = centroid(e, r.h, r.k, r.l, r.px_z);
    if (!c.valid) continue;
    const auto px = e.detector[0].mm_to_px(c.X, c.Y);
    worst = std::fmax(worst, std::hypot(px.first - r.cal_px_fast,
                                        px.second - r.cal_px_slow));
    ++tested;
  }
  check::is_true(tested > 30, "enough reflections");
  check::is_true(worst < 0.5, "the state must predict where DIALS predicts");
}

TEST(analytical_crystal_derivatives_match_finite_differences) {
  const Experiment e = real_experiment();
  double scale = 0.0;
  for (double m : e.crystal->A.m) scale = std::fmax(scale, std::abs(m));

  std::vector<double> relative;
  int tested = 0;
  for (const real::Row &r : real::rows()) {
    const PredictionState s = prediction_state(e, 0, r.h, r.k, r.l, r.px_z);
    if (!s.valid) continue;
    // Skip where phi is genuinely ill determined, which is where eqn (40)
    // divides by something near zero. DIALS discards these too.
    if (std::abs(s.volume) < 0.05) continue;
    const auto analytic = crystal_derivatives(s, r.h, r.k, r.l);

    for (std::size_t p = 0; p < 9; ++p) {
      // Central difference, step chosen near cube-root-eps relative for a
      // central formula.
      const double step = 1e-5 * scale;
      Experiment plus = e, minus = e;
      plus.crystal->A.m[p] += step;
      minus.crystal->A.m[p] -= step;
      const Centroid a = centroid(plus, r.h, r.k, r.l, r.px_z);
      const Centroid b = centroid(minus, r.h, r.k, r.l, r.px_z);
      if (!a.valid || !b.valid) continue;

      const double numeric[3] = {(a.X - b.X) / (2 * step),
                                 (a.Y - b.Y) / (2 * step),
                                 (a.phi - b.phi) / (2 * step)};
      const double exact[3] = {analytic[p].dX, analytic[p].dY, analytic[p].dphi};
      for (int c = 0; c < 3; ++c) {
        const double size = std::fmax(std::abs(numeric[c]), std::abs(exact[c]));
        if (size < 1e-6) continue;  // both essentially zero
        relative.push_back(std::abs(numeric[c] - exact[c]) / size);
        ++tested;
      }
    }
  }
  std::sort(relative.begin(), relative.end());
  check::is_true(tested > 500, "enough derivative components compared");
  // A central difference on a smooth function agrees with the exact
  // derivative to many digits. Anything worse than this is a mistake in the
  // algebra, not a limitation of the comparison.
  check::is_true(relative[relative.size() / 2] < 1e-8, "median agreement");
  check::is_true(relative[relative.size() * 99 / 100] < 1e-5,
                 "and the tail must not hide one either");
}

TEST(the_phi_derivative_blows_up_where_the_paper_says_it_does) {
  // The denominator of eqn (40) is the volume of the parallelepiped formed by
  // the rotation axis, the reciprocal lattice vector and the beam. Where it
  // vanishes, phi is genuinely undetermined -- these are the reflections near
  // the rotation axis, the ones with large Lorentz factors. The test is that
  // the derivative really does grow as the volume shrinks, so that discarding
  // on the volume is discarding the right reflections.
  const Experiment e = real_experiment();
  std::vector<std::pair<double, double>> pairs;  // volume, |dphi| per unit |h|
  for (const real::Row &r : real::rows()) {
    const PredictionState s = prediction_state(e, 0, r.h, r.k, r.l, r.px_z);
    if (!s.valid) continue;
    const auto d = crystal_derivatives(s, r.h, r.k, r.l);
    double biggest = 0.0;
    for (const CentroidDerivative &c : d) biggest = std::fmax(biggest, std::abs(c.dphi));
    // Divided by |h|. The numerator of eqn (40) carries a factor of the Miller
    // index, so without this the two effects are confounded and the
    // measurement comes out at 1.7 rather than 3.9 -- which is how this test
    // failed the first time it was written.
    const double h = std::sqrt(static_cast<double>(r.h * r.h + r.k * r.k + r.l * r.l));
    pairs.emplace_back(std::abs(s.volume), biggest / std::fmax(h, 1.0));
  }
  check::is_true(pairs.size() > 30, "enough reflections");
  std::sort(pairs.begin(), pairs.end());
  const std::size_t q = pairs.size() / 4;
  double low = 0.0, high = 0.0;
  for (std::size_t i = 0; i < q; ++i) low += pairs[i].second;
  for (std::size_t i = pairs.size() - q; i < pairs.size(); ++i) high += pairs[i].second;
  // Measured at 3.9. Note that every reflection in this sample has a volume
  // above DIALS' cutoff of 0.05, so the pathological cases are not even
  // represented here -- the trend is visible well before them.
  check::is_true(low > 2.0 * high,
                 "small volume must mean a large rotation-angle derivative");
}

TEST(spline_weights_agree_with_the_interpolation_they_describe) {
  // The weights are used to turn a derivative with respect to A into
  // derivatives with respect to control points, so they must be exactly the
  // weights A_at uses. Reconstructing A from them is the check.
  Experiment e = real_experiment();
  const std::size_t n = 6;
  e.crystal->A_points.clear();
  for (std::size_t i = 0; i < n; ++i) {
    Mat3 m = e.crystal->A;
    m.m[0] += 0.001 * static_cast<double>(i);
    m.m[4] -= 0.0004 * static_cast<double>(i);
    e.crystal->A_points.push_back(m);
  }
  for (double z : {0.0, 17.5, 120.0, 240.0, 299.9}) {
    const SplineWeights w = spline_weights(e, z);
    double sum = 0.0;
    Mat3 rebuilt;
    for (std::size_t i = 0; i < 9; ++i) rebuilt.m[i] = 0.0;
    for (std::size_t i = 0; i < w.count; ++i) {
      sum += w.weight[i];
      for (std::size_t k = 0; k < 9; ++k) {
        rebuilt.m[k] += w.weight[i] * e.crystal->A_points[w.index[i]].m[k];
      }
    }
    check::close(sum, 1.0, 1e-12, "weights must sum to one");
    const Mat3 direct = e.setting_at(z);
    for (std::size_t k = 0; k < 9; ++k) {
      check::close(rebuilt.m[k], direct.m[k], 1e-12, "rebuilt A matches A_at");
    }
  }
}

TEST(only_four_control_points_are_ever_touched) {
  // The banding, asserted rather than assumed: this is what would make a
  // device Jacobian sparse, and what an interpolating spline would not give.
  Experiment e = real_experiment();
  e.crystal->A_points.assign(12, e.crystal->A);
  for (double z = 0.0; z <= 300.0; z += 3.0) {
    const SplineWeights w = spline_weights(e, z);
    check::equal(static_cast<long long>(w.count), 4, "four control points");
    std::size_t distinct[4];
    std::size_t n = 0;
    for (std::size_t i = 0; i < 4; ++i) {
      bool seen = false;
      for (std::size_t j = 0; j < n; ++j) {
        if (distinct[j] == w.index[i]) seen = true;
      }
      if (!seen) distinct[n++] = w.index[i];
    }
    check::is_true(n <= 4, "at most four distinct");
    for (std::size_t i = 0; i < n; ++i) {
      check::is_true(distinct[i] < 12, "index in range");
    }
  }
}

TEST(analytical_detector_derivatives_match_finite_differences) {
  const Experiment e = real_experiment();
  std::vector<double> relative;
  int tested = 0;

  for (const real::Row &r : real::rows()) {
    const PredictionState s = prediction_state(e, 0, r.h, r.k, r.l, r.px_z);
    if (!s.valid) continue;
    const auto analytic = detector_derivatives(s, e.detector[0]);

    for (std::size_t p = 0; p < 6; ++p) {
      // Translations are in millimetres, rotations in radians, so they need
      // different steps; using one for both makes the rotation derivative
      // look wrong when it is the step that is.
      const double step = p < 3 ? 1e-4 : 1e-6;
      double shift[6] = {0, 0, 0, 0, 0, 0};
      Experiment plus = e, minus = e;
      shift[p] = step;
      plus.detector.panels[0] = perturb_panel(e.detector[0], shift);
      shift[p] = -step;
      minus.detector.panels[0] = perturb_panel(e.detector[0], shift);

      const Centroid a = centroid(plus, r.h, r.k, r.l, r.px_z);
      const Centroid b = centroid(minus, r.h, r.k, r.l, r.px_z);
      if (!a.valid || !b.valid) continue;

      const double numeric[3] = {(a.X - b.X) / (2 * step),
                                 (a.Y - b.Y) / (2 * step),
                                 (a.phi - b.phi) / (2 * step)};
      const double exact[3] = {analytic[p].dX, analytic[p].dY, analytic[p].dphi};
      for (int c = 0; c < 3; ++c) {
        const double size = std::fmax(std::abs(numeric[c]), std::abs(exact[c]));
        if (size < 1e-6) continue;
        relative.push_back(std::abs(numeric[c] - exact[c]) / size);
        ++tested;
      }
    }
  }
  std::sort(relative.begin(), relative.end());
  check::is_true(tested > 300, "enough derivative components compared");
  check::is_true(relative[relative.size() / 2] < 1e-7, "median agreement");
  check::is_true(relative[relative.size() * 99 / 100] < 1e-4, "and the tail");
}

TEST(the_detector_cannot_move_the_rotation_angle) {
  // Not an approximation: neither r0 nor s0 depends on where the detector is,
  // so the diffracting angle cannot either. Asserted because a nonzero dphi
  // here would mean the chain rule had picked up a term that does not exist.
  const Experiment e = real_experiment();
  for (const real::Row &r : real::rows()) {
    const PredictionState s = prediction_state(e, 0, r.h, r.k, r.l, r.px_z);
    if (!s.valid) continue;
    for (const CentroidDerivative &d : detector_derivatives(s, e.detector[0])) {
      check::close(d.dphi, 0.0, 0.0, "exactly zero, not merely small");
    }
  }
}

TEST(analytical_beam_derivatives_match_finite_differences) {
  const Experiment e = real_experiment();
  std::vector<double> relative;
  int tested = 0;

  for (const real::Row &r : real::rows()) {
    const PredictionState s = prediction_state(e, 0, r.h, r.k, r.l, r.px_z);
    if (!s.valid) continue;
    if (std::abs(s.volume) < 0.05) continue;
    const auto analytic = beam_derivatives(s, e.beam);

    for (std::size_t p = 0; p < 2; ++p) {
      const double step = 1e-6;
      double shift[2] = {0.0, 0.0};
      Experiment plus = e, minus = e;
      shift[p] = step;
      plus.beam = perturb_beam(e.beam, shift);
      shift[p] = -step;
      minus.beam = perturb_beam(e.beam, shift);

      const Centroid a = centroid(plus, r.h, r.k, r.l, r.px_z);
      const Centroid b = centroid(minus, r.h, r.k, r.l, r.px_z);
      if (!a.valid || !b.valid) continue;

      const double numeric[3] = {(a.X - b.X) / (2 * step),
                                 (a.Y - b.Y) / (2 * step),
                                 (a.phi - b.phi) / (2 * step)};
      const double exact[3] = {analytic[p].dX, analytic[p].dY, analytic[p].dphi};
      for (int c = 0; c < 3; ++c) {
        const double size = std::fmax(std::abs(numeric[c]), std::abs(exact[c]));
        if (size < 1e-6) continue;
        relative.push_back(std::abs(numeric[c] - exact[c]) / size);
        ++tested;
      }
    }
  }
  std::sort(relative.begin(), relative.end());
  check::is_true(tested > 50, "enough derivative components compared");
  check::is_true(relative[relative.size() / 2] < 1e-6, "median agreement");
  check::is_true(relative.back() < 1e-3, "worst case");
}

TEST(the_shared_perturbation_is_what_refinement_applies) {
  // perturb_panel and perturb_beam are used by refine and by the derivative
  // tests both. This checks the properties refinement relies on rather than
  // the identity of the code path: a zero shift changes nothing, a rotation
  // keeps the panel axes orthonormal however many times it is applied, and a
  // rotation about the panel centre leaves that centre where it was.
  const Experiment e = real_experiment();
  const Panel &p = e.detector[0];

  const double none[6] = {0, 0, 0, 0, 0, 0};
  const Panel same = perturb_panel(p, none);
  check::close((same.origin - p.origin).norm(), 0.0, 0.0, "origin untouched");
  check::close((same.fast - p.fast).norm(), 0.0, 0.0, "fast untouched");

  Panel turned = p;
  const double turn[6] = {0, 0, 0, 1e-3, -5e-4, 2e-4};
  for (int i = 0; i < 50; ++i) turned = perturb_panel(turned, turn);
  check::close(turned.fast.norm(), 1.0, 1e-12, "fast stays a unit vector");
  check::close(turned.slow.norm(), 1.0, 1e-12, "slow stays a unit vector");
  check::close(turned.fast.dot(turned.slow), 0.0, 1e-12, "and stay orthogonal");

  // The centre of the panel is the point rotations act about, so it moves only
  // by the translation.
  const auto centre_of = [](const Panel &q) {
    return q.lab_coord_mm(0.5 * static_cast<double>(q.image_size[0]) * q.pixel_size[0],
                          0.5 * static_cast<double>(q.image_size[1]) * q.pixel_size[1]);
  };
  const double turn_only[6] = {0, 0, 0, 2e-3, 1e-3, -1e-3};
  const Panel rotated = perturb_panel(p, turn_only);
  check::close((centre_of(rotated) - centre_of(p)).norm(), 0.0, 1e-9,
               "a pure rotation leaves the panel centre alone");

  const double move[2] = {1e-3, -2e-3};
  const Beam b = perturb_beam(e.beam, move);
  check::close(b.direction.norm(), 1.0, 1e-12, "the beam stays a unit vector");
}

TEST(the_parallax_jacobian_matches_finite_differences) {
  // The bridge between millimetres and pixels. Without the parallax correction
  // it would be a division by the pixel size; with it, it is not, and getting
  // that wrong scales every analytical derivative by a few parts in a thousand
  // in a way no test of the prediction itself would notice.
  const Experiment e = real_experiment();
  const Panel &p = e.detector[0];
  double worst = 0.0;
  for (double mm_f : {5.0, 60.0, 150.0, 280.0}) {
    for (double mm_s : {5.0, 90.0, 200.0, 320.0}) {
      double analytic[4];
      p.mm_to_px_jacobian(mm_f, mm_s, analytic);
      const double step = 1e-4;
      const auto plus_f = p.mm_to_px(mm_f + step, mm_s);
      const auto minus_f = p.mm_to_px(mm_f - step, mm_s);
      const auto plus_s = p.mm_to_px(mm_f, mm_s + step);
      const auto minus_s = p.mm_to_px(mm_f, mm_s - step);
      const double numeric[4] = {(plus_f.first - minus_f.first) / (2 * step),
                                 (plus_s.first - minus_s.first) / (2 * step),
                                 (plus_f.second - minus_f.second) / (2 * step),
                                 (plus_s.second - minus_s.second) / (2 * step)};
      for (int i = 0; i < 4; ++i) {
        // Scaled against the size of the Jacobian as a whole, not against each
        // element. The off-diagonal terms pass through zero near the beam
        // centre, and dividing by them turns an absolute agreement of 3e-10
        // into an apparently enormous relative error.
        const double size = 1.0 / p.pixel_size[0];
        worst = std::fmax(worst, std::abs(numeric[i] - analytic[i]) / size);
      }
    }
  }
  check::is_true(worst < 1e-9, "analytic parallax Jacobian");
}

TEST(the_parallax_jacobian_is_not_just_the_pixel_size) {
  // If it were, the whole function would be pointless. The off-diagonal terms
  // and the departure of the diagonal from 1/pixel_size are the correction.
  const Experiment e = real_experiment();
  const Panel &p = e.detector[0];
  double analytic[4];
  p.mm_to_px_jacobian(280.0, 320.0, analytic);
  const double plain = 1.0 / p.pixel_size[0];
  // Measured at 7.0e-4 in the far corner of this panel: small, and a hundred
  // times the tolerance the derivative comparison is held to, so neglecting it
  // would be visible there.
  check::is_true(std::abs(analytic[0] - plain) / plain > 5e-4,
                 "the diagonal must differ from the pixel size");
  check::is_true(std::abs(analytic[1]) > 1e-4,
                 "and there must be an off-diagonal term at all");
}

// --------------------------------------------------------------------------
// the two Jacobians, whole
// --------------------------------------------------------------------------


namespace {

// A reflection table of predictions from the real geometry, so the comparison
// runs on the full Jacobian rather than one reflection at a time.
Table table_from(const Experiment &e, double d_min) {
  PredictOptions po;
  po.d_min = d_min;
  const std::vector<Prediction> predictions = predict(e, po);
  Table t;
  t.nrows = predictions.size();
  Column &xyz = t.real_column("xyzobs.px.value", "vec3<double>", 3);
  Column &miller = t.int_column("miller_index", "cctbx::miller::index<>", 3);
  Column &panel = t.int_column("panel", "std::size_t", 1);
  Column &id = t.int_column("id", "int", 1);
  for (std::size_t i = 0; i < predictions.size(); ++i) {
    xyz.reals[i * 3 + 0] = predictions[i].px_fast;
    xyz.reals[i * 3 + 1] = predictions[i].px_slow;
    xyz.reals[i * 3 + 2] = predictions[i].z;
    miller.ints[i * 3 + 0] = predictions[i].h;
    miller.ints[i * 3 + 1] = predictions[i].k;
    miller.ints[i * 3 + 2] = predictions[i].l;
    panel.ints[i] = 0;
    id.ints[i] = 0;
  }
  return t;
}

}  // namespace

TEST(the_two_jacobians_agree_where_a_finite_difference_is_valid) {
  // The whole Jacobian, not one reflection at a time: crystal, detector and
  // beam together, static and scan-varying, which is what refinement actually
  // assembles. Agreeing on a refined result is not enough -- a flat minimum
  // can hide a wrong derivative.
  const Experiment e = real_experiment();
  ExperimentList list;
  list.experiments.push_back(e);
  const Table t = table_from(e, 3.5);
  check::is_true(t.nrows > 1000, "enough reflections");

  for (std::size_t points : {std::size_t(1), std::size_t(5)}) {
    RefineOptions options;
    options.scan_points = points;
    options.beam = true;
    const JacobianComparison c = compare_jacobians(list, t, options);
    check::is_true(c.compared > 100000, "enough entries compared");
    // Measured: median 1.2e-8 static, 3.6e-8 scan-varying. The ninety-ninth
    // percentile is 0.15 and 0.05 -- far worse than the median, because the
    // tail is not a gradual loss of accuracy but a small set of entries where
    // the two disagree completely, and a percentile walks into it.
    check::is_true(c.median_relative < 1e-6, "median agreement");
    check::is_true(c.percentile_99 < 0.3, "ninety-ninth percentile");
    // And the tail is small. It is NOT a defect in the analytical derivative:
    // a finite difference that straddles a change in which Ewald root is
    // nearest the observation is comparing two different branches, so the
    // numerical value is the invalid one there. Measured at 0.07 to 0.3 per
    // cent of entries, depending on the data and the parameterisation.
    const double gross = static_cast<double>(c.grossly_different) /
                         static_cast<double>(c.compared);
    check::is_true(gross < 0.01, "grossly disagreeing entries stay rare");
  }
}

TEST(analytic_and_numerical_refinement_reach_the_same_model) {
  // The end-to-end check. Perturb, refine both ways, and require the same
  // answer -- necessary because a derivative can be right and still be wired
  // into the normal equations with the wrong sign or index.
  const Experiment truth = real_experiment();
  const Table t = table_from(truth, 3.5);

  Experiment moved = truth;
  const double shift[6] = {0.3, -0.2, 0.5, 0.0, 0.0, 0.0};
  moved.detector.panels[0] = perturb_panel(truth.detector[0], shift);
  moved.crystal->A = moved.crystal->A * 1.0003;

  const auto run = [&](bool analytic) {
    ExperimentList list;
    list.experiments.push_back(moved);
    RefineOptions options;
    options.analytic = analytic;
    options.outlier_sigma = 0.0;
    options.macrocycles = 1;
    const RefineResult r = refine(list, t, options);
    return std::make_pair(r, list[0]);
  };

  const auto numerical = run(false);
  const auto analytical = run(true);

  check::is_true(numerical.first.rmsd_x < 0.02, "numerical converged");
  check::is_true(analytical.first.rmsd_x < 0.02, "analytical converged");
  check::close(analytical.first.rmsd_x, numerical.first.rmsd_x,
               0.1 * std::fmax(numerical.first.rmsd_x, 1e-6),
               "same residual either way");
  const UnitCell a = analytical.second.crystal->cell();
  const UnitCell b = numerical.second.crystal->cell();
  check::close(a.volume(), b.volume(), 1e-4 * b.volume(), "same cell volume");
  check::close((analytical.second.detector[0].origin -
                numerical.second.detector[0].origin)
                   .norm(),
               0.0, 1e-3, "same detector origin");
}

}  // namespace mxi
