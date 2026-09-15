// Refinement, tested by perturbing a model and requiring it back.
//
// The pattern throughout: take a known experiment, predict from it, perturb
// the model, and refine against the predictions. The answer is the model that
// went in, exactly, because the data was generated from it. A refinement that
// reduces its residual without recovering the truth has found a different
// minimum, and only ground truth distinguishes the two.

#include <algorithm>
#include <cmath>

#include "../src/expt.h"
#include "../src/predict.h"
#include "../src/derivatives.h"
#include "../src/refine.h"
#include "../src/refl.h"
#include "check.h"

using namespace mxi;

namespace {

Experiment base_experiment() {
  Experiment e;
  e.beam.direction = {0.0, 0.0, 1.0};
  e.beam.wavelength = 0.9537;

  Panel p;
  p.fast = {1.0, 0.0, 0.0};
  p.slow = {0.0, -1.0, 0.0};
  p.pixel_size[0] = p.pixel_size[1] = 0.075;
  p.image_size[0] = 4148;
  p.image_size[1] = 4362;
  p.origin = {-0.075 * 2074.0, 0.075 * 2181.0, -200.0};
  p.parallax = true;
  p.mu = 3.663;
  p.thickness = 0.45;
  e.detector.panels.push_back(p);

  e.goniometer.axis = {1.0, 0.0, 0.0};
  e.scan.first_image = 1;
  e.scan.last_image = 900;
  e.scan.osc_start = 0.0;
  e.scan.osc_width = 0.1;

  const Mat3 r = rotation({0.3, -0.5, 0.81}, 0.7);
  e.crystal = Crystal::from_real_space(r * Vec3{78.0, 0.0, 0.0},
                                       r * Vec3{0.0, 78.0, 0.0},
                                       r * Vec3{0.0, 0.0, 78.0});
  return e;
}

// Predictions from `truth`, presented as observations with their indices.
Table observations_from(const Experiment &truth, double d_min) {
  PredictOptions po;
  po.d_min = d_min;
  const std::vector<Prediction> predictions = predict(truth, po);

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
    panel.ints[i] = static_cast<std::int64_t>(predictions[i].panel);
    id.ints[i] = 0;
  }
  return t;
}

double cell_difference(const UnitCell &a, const UnitCell &b) {
  return std::fmax(std::fmax(std::abs(a.a - b.a), std::abs(a.b - b.b)),
                   std::abs(a.c - b.c));
}

}  // namespace

TEST(residual_is_zero_for_an_exact_model) {
  const Experiment truth = base_experiment();
  const Table t = observations_from(truth, 3.0);
  const Column &xyz = t.at("xyzobs.px.value");
  const Column &miller = t.at("miller_index");
  check::is_true(t.nrows > 1000, "enough predictions");

  double worst = 0.0;
  for (std::size_t i = 0; i < t.nrows; ++i) {
    const Residual r = centroid_residual(
        truth, 0, static_cast<int>(miller.integer(i, 0)),
        static_cast<int>(miller.integer(i, 1)),
        static_cast<int>(miller.integer(i, 2)), xyz.real(i, 0), xyz.real(i, 1),
        xyz.real(i, 2));
    check::is_true(r.valid, "every prediction must be reproducible");
    worst = std::fmax(worst, std::abs(r.dx));
    worst = std::fmax(worst, std::abs(r.dy));
    worst = std::fmax(worst, std::abs(r.dz));
  }
  check::close(worst, 0.0, 1e-8, "residual against the model that made it");
}

TEST(residual_picks_the_root_nearer_the_observation) {
  // Both Ewald roots are valid predictions of the same reflection at different
  // angles. Choosing by proximity rather than by the entering flag is what
  // makes refinement work from a model that is still far out, when the flag
  // cannot be trusted.
  const Experiment truth = base_experiment();
  const Table t = observations_from(truth, 3.5);
  const Column &xyz = t.at("xyzobs.px.value");
  const Column &miller = t.at("miller_index");
  std::size_t large = 0;
  for (std::size_t i = 0; i < t.nrows; ++i) {
    const Residual r = centroid_residual(
        truth, 0, static_cast<int>(miller.integer(i, 0)),
        static_cast<int>(miller.integer(i, 1)),
        static_cast<int>(miller.integer(i, 2)), xyz.real(i, 0), xyz.real(i, 1),
        xyz.real(i, 2));
    if (r.valid && std::abs(r.dz) > 1.0) ++large;
  }
  // Picking the wrong root puts a reflection tens or hundreds of images away,
  // so any such failure would show here immediately.
  check::equal(static_cast<long long>(large), 0, "no reflection on the wrong root");
}

TEST(refinement_recovers_a_perturbed_detector) {
  const Experiment truth = base_experiment();
  const Table t = observations_from(truth, 3.0);

  ExperimentList list;
  Experiment moved = truth;
  moved.detector.panels[0].origin += Vec3{0.4, -0.3, 1.2};
  list.experiments.push_back(moved);

  RefineOptions options;
  options.crystal = false;  // only the detector was moved
  options.outlier_sigma = 0.0;
  const RefineResult result = refine(list, t, options);

  check::is_true(result.n_used > 1000, "most reflections used");
  check::is_true(result.rmsd_x < 0.01 && result.rmsd_y < 0.01,
                 "residuals must come back to nothing");
  const Vec3 recovered = list[0].detector[0].origin;
  check::close((recovered - truth.detector[0].origin).norm(), 0.0, 0.01,
               "detector origin recovered");
}

TEST(refinement_recovers_a_perturbed_cell) {
  const Experiment truth = base_experiment();
  const Table t = observations_from(truth, 3.0);

  ExperimentList list;
  Experiment stretched = truth;
  // A tenth of a per cent on the cell, which is far more than refinement
  // should ever be asked to absorb and still an easy target.
  stretched.crystal->A = stretched.crystal->A * 1.001;
  list.experiments.push_back(stretched);

  RefineOptions options;
  options.detector = false;  // only the cell was changed
  options.outlier_sigma = 0.0;
  const RefineResult result = refine(list, t, options);

  check::is_true(result.rmsd_x < 0.01, "residual in x");
  check::is_true(cell_difference(list[0].crystal->cell(), truth.crystal->cell()) <
                     0.01,
                 "cell recovered to better than 0.01 Angstrom");
}

TEST(refinement_improves_a_model_wrong_in_both) {
  const Experiment truth = base_experiment();
  const Table t = observations_from(truth, 3.0);

  ExperimentList list;
  Experiment wrong = truth;
  wrong.detector.panels[0].origin += Vec3{0.0, 0.0, 0.8};
  wrong.crystal->A = wrong.crystal->A * 1.0004;
  list.experiments.push_back(wrong);

  ExperimentList before = list;
  RefineOptions options;
  options.outlier_sigma = 0.0;
  options.macrocycles = 1;
  const RefineResult result = refine(list, t, options);

  // Detector distance and cell scale are nearly degenerate -- moving the
  // detector further away and enlarging the cell predict almost the same
  // positions -- so what is asserted is that the residual is small, not that
  // each parameter came back individually. On real data that degeneracy is
  // the reason a refined cell can sit a part in a thousand from another
  // program's while predicting the same spots to a fifth of a pixel.
  check::is_true(result.rmsd_x < 0.05 && result.rmsd_y < 0.05,
                 "the fit must be good even if the parameters trade off");
  check::is_true(result.rmsd_z < 0.05, "and in the rotation angle");
}

TEST(outliers_are_rejected_rather_than_fitted) {
  const Experiment truth = base_experiment();
  Table t = observations_from(truth, 3.0);
  // Move fifty reflections a long way, as a second lattice or a misindexed
  // spot would be.
  Column &xyz = t.real_column("xyzobs.px.value", "vec3<double>", 3);
  {
    const Table fresh = observations_from(truth, 3.0);
    xyz.reals = fresh.at("xyzobs.px.value").reals;
  }
  for (std::size_t i = 0; i < 50; ++i) {
    xyz.reals[i * 3 + 0] += 25.0;
    xyz.reals[i * 3 + 1] -= 18.0;
  }

  ExperimentList list;
  list.experiments.push_back(truth);
  RefineOptions options;
  options.macrocycles = 3;
  options.outlier_sigma = 4.0;
  const RefineResult result = refine(list, t, options);

  check::is_true(result.n_rejected >= 50, "the planted outliers must go");
  // And with them gone the fit is exact again, which it would not be if they
  // had been absorbed.
  check::is_true(result.rmsd_x < 0.05 && result.rmsd_y < 0.05,
                 "residual after rejection");
}

TEST(separate_crystals_fit_sweeps_that_one_crystal_cannot) {
  // Two sweeps whose goniometer settings are not quite what the file says --
  // which is the situation in every real multi-sweep experiment, and the
  // reason refinement's first job is to break the single-crystal constraint.
  const Experiment truth = base_experiment();
  Experiment second = truth;
  second.goniometer.fixed = rotation({0.0, 1.0, 0.0}, Scan::radians(30.0));

  Table a = observations_from(truth, 3.5);
  // The second sweep really has a slightly different orientation, as a
  // goniometer that does not return to the same place would produce.
  Experiment drifted = second;
  drifted.crystal->A = rotation({0.1, 0.2, 0.97}, Scan::radians(0.15)) *
                       drifted.crystal->A;
  Table b = observations_from(drifted, 3.5);

  Table both;
  both.nrows = a.nrows + b.nrows;
  Column &xyz = both.real_column("xyzobs.px.value", "vec3<double>", 3);
  Column &miller = both.int_column("miller_index", "cctbx::miller::index<>", 3);
  Column &panel = both.int_column("panel", "std::size_t", 1);
  Column &id = both.int_column("id", "int", 1);
  const auto copy = [&](const Table &from, std::size_t offset, std::int64_t which) {
    for (std::size_t i = 0; i < from.nrows; ++i) {
      for (std::size_t k = 0; k < 3; ++k) {
        xyz.reals[(offset + i) * 3 + k] = from.at("xyzobs.px.value").real(i, k);
        miller.ints[(offset + i) * 3 + k] = from.at("miller_index").integer(i, k);
      }
      panel.ints[offset + i] = 0;
      id.ints[offset + i] = which;
    }
  };
  copy(a, 0, 0);
  copy(b, a.nrows, 1);

  ExperimentList list;
  list.experiments.push_back(truth);
  list.experiments.push_back(second);

  RefineOptions shared;
  shared.detector = false;
  shared.outlier_sigma = 0.0;
  shared.shared_crystal = true;
  ExperimentList one = list;
  const RefineResult with_one = refine(one, both, shared);

  RefineOptions separate = shared;
  separate.shared_crystal = false;
  ExperimentList many = list;
  const RefineResult with_many = refine(many, both, separate);

  check::is_true(with_many.rmsd_x < with_one.rmsd_x,
                 "separate crystals must fit better than one");
  // And the improvement should be substantial, not marginal: the sweeps really
  // do have different orientations.
  check::is_true(with_many.rmsd_x * 2.0 < with_one.rmsd_x,
                 "and by a wide margin");
}

TEST(update_predictions_writes_a_usable_xyzcal) {
  const Experiment truth = base_experiment();
  Table t = observations_from(truth, 3.5);
  ExperimentList list;
  list.experiments.push_back(truth);
  update_predictions(list, t);

  check::is_true(t.has("xyzcal.px"), "xyzcal.px written");
  const Column &cal = t.at("xyzcal.px");
  const Column &obs = t.at("xyzobs.px.value");
  double worst = 0.0;
  for (std::size_t i = 0; i < t.nrows; ++i) {
    worst = std::fmax(worst, std::abs(cal.real(i, 0) - obs.real(i, 0)));
    worst = std::fmax(worst, std::abs(cal.real(i, 1) - obs.real(i, 1)));
  }
  check::close(worst, 0.0, 1e-8, "calculated equals observed for exact data");
}

// --------------------------------------------------------------------------
// scan-varying crystals
// --------------------------------------------------------------------------

TEST(a_static_crystal_stays_static_under_scan_varying_refinement) {
  // The control points are free to drift and the truth does not. If they
  // wander anyway they are absorbing noise, and on real data they would be
  // absorbing the detector.
  const Experiment truth = base_experiment();
  const Table t = observations_from(truth, 3.0);
  ExperimentList list;
  list.experiments.push_back(truth);

  RefineOptions options;
  options.outlier_sigma = 0.0;
  options.scan_points = 5;
  options.macrocycles = 1;
  // No rejection needed: with a static truth there are no near-tangential
  // failures to reject, which is itself the point -- they appear only once the
  // model varies over the scan.
  const RefineResult result = refine(list, t, options);

  check::is_true(result.rmsd_x < 0.02, "still fits");
  check::is_true(list[0].crystal->scan_varying(), "the model is scan-varying");
  double worst = 0.0;
  for (const Mat3 &A : list[0].crystal->A_points) {
    for (std::size_t k = 0; k < 9; ++k) {
      worst = std::fmax(worst, std::abs(A.m[k] - truth.crystal->A.m[k]));
    }
  }
  // The elements of A are around 1/78, so this is a part in ten thousand.
  check::is_true(worst < 1e-6, "control points must not wander off the truth");
}

TEST(scan_varying_refinement_recovers_a_drifting_crystal) {
  // A crystal that rotates by a tenth of a degree across the scan, which is
  // the kind of drift a real goniometer and a real sample produce. A static
  // model cannot express it; the whole point of the scan-varying one is that
  // it can.
  Experiment truth = base_experiment();
  const std::size_t points = 5;
  truth.crystal->A_points.clear();
  for (std::size_t i = 0; i < points; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(points - 1);
    truth.crystal->A_points.push_back(
        rotation({0.2, 0.9, -0.39}, Scan::radians(0.10 * t)) * truth.crystal->A);
  }
  const Table observations = observations_from(truth, 3.0);
  check::is_true(observations.nrows > 1000, "enough predictions");

  ExperimentList list;
  Experiment start = truth;
  start.crystal->A_points.clear();  // begin from the static matrix
  list.experiments.push_back(start);

  // Outlier rejection is REQUIRED here, not incidental. About four per cent of
  // reflections cross the Ewald sphere near-tangentially, where the angle at
  // which they diffract is enormously sensitive to the crystal model; once the
  // model varies over the scan, the forward and reverse maps choose different
  // roots for those and their residuals are tens of images. They are the same
  // reflections that carry huge Lorentz factors. With them in, the fit is
  // dragged and the drift is not recovered at all; with them out, the planted
  // drift comes back to a per cent.
  RefineOptions statically;
  statically.outlier_sigma = 4.0;
  statically.macrocycles = 3;
  ExperimentList fixed = list;
  const RefineResult without = refine(fixed, observations, statically);

  RefineOptions varying = statically;
  // Five, matching the truth's own parameterisation. Under linear
  // interpolation this had to be dropped to three, because the outermost
  // control points of a finer model were constrained only by whatever lay at
  // the very ends of the scan and wandered. The B-spline does not have that
  // problem: its end control points are clamped to the curve, so the data near
  // the ends pins them directly, and five recovers the planted drift exactly
  // where three can only approximate it.
  varying.scan_points = points;
  const RefineResult with = refine(list, observations, varying);

  check::is_true(with.rmsd_x < 0.02 * without.rmsd_x,
                 "scan-varying must beat static on a drifting crystal");

  // And recover the drift itself, not merely fit it: the rotation between the
  // first and last control point should be the tenth of a degree put in.
  const Mat3 first = list[0].crystal->A_points.front();
  const Mat3 last = list[0].crystal->A_points.back();
  const double turn = Scan::degrees(rotation_angle(last * first.inverse()));
  // Exact, not approximate: the model being fitted is the model that made the
  // data, so anything short of exact recovery is a defect in the refinement.
  check::close(turn, 0.10, 0.001, "recovered drift across the scan");
}

TEST(scan_varying_models_survive_a_file_round_trip) {
  Experiment e = base_experiment();
  e.crystal->A_points.assign(4, e.crystal->A);
  e.crystal->A_points[3] = rotation({0, 0, 1}, Scan::radians(0.2)) * e.crystal->A;
  ExperimentList list;
  list.experiments.push_back(e);

  const ExperimentList back = experiments_from_json(experiments_to_json(list));
  check::is_true(back[0].crystal->scan_varying(), "still scan-varying");
  // Written per image, as DIALS does, so the control points are expanded --
  // what must survive is the model, not the parameterisation.
  double worst = 0.0;
  for (double t : {0.0, 0.25, 0.5, 0.75, 0.999}) {
    const Mat3 a = list[0].crystal->A_at(t);
    const Mat3 b = back[0].crystal->A_at(t);
    for (std::size_t k = 0; k < 9; ++k) {
      worst = std::fmax(worst, std::abs(a.m[k] - b.m[k]));
    }
  }
  check::is_true(worst < 1e-6, "A(t) must survive the round trip");
}

TEST(prediction_and_the_refinement_target_agree_exactly) {
  // The forward map (predict) and the reverse map (centroid_residual) must be
  // the same model seen from two directions. They were not for a scan-varying
  // crystal: the iteration that finds the setting matrix was seeded from the
  // start of the scan for both Ewald roots, so a reflection late in the sweep
  // converged on the wrong one. The existing round-trip test could not see it,
  // because it goes through reciprocal_lattice_point rather than through the
  // target function refinement actually minimises.
  //
  // Static first, where agreement is exact for every reflection.
  const Experiment truth = base_experiment();
  const Table t = observations_from(truth, 4.0);
  const Column &xyz = t.at("xyzobs.px.value");
  const Column &miller = t.at("miller_index");
  double worst = 0.0;
  for (std::size_t i = 0; i < t.nrows; ++i) {
    const Residual r = centroid_residual(
        truth, 0, static_cast<int>(miller.integer(i, 0)),
        static_cast<int>(miller.integer(i, 1)),
        static_cast<int>(miller.integer(i, 2)), xyz.real(i, 0), xyz.real(i, 1),
        xyz.real(i, 2));
    if (!r.valid) continue;
    worst = std::fmax(worst, std::hypot(r.dx, r.dy));
  }
  check::close(worst, 0.0, 1e-9, "static: every reflection must agree exactly");
}

TEST(a_scan_varying_model_agrees_except_where_the_angle_is_ill_conditioned) {
  // With a drifting crystal the two maps agree exactly for the great majority
  // and disagree wildly for a few per cent. Those few are reflections crossing
  // the Ewald sphere near-tangentially, where the diffracting angle is hugely
  // sensitive to the model -- the same population that carries large Lorentz
  // factors. This test pins both halves of that: the bulk must be exact, and
  // the tail must be small.
  Experiment truth = base_experiment();
  truth.crystal->A_points.clear();
  for (std::size_t i = 0; i < 5; ++i) {
    const double t = static_cast<double>(i) / 4.0;
    truth.crystal->A_points.push_back(
        rotation({0.2, 0.9, -0.39}, Scan::radians(0.10 * t)) * truth.crystal->A);
  }
  const Table t = observations_from(truth, 4.0);
  const Column &xyz = t.at("xyzobs.px.value");
  const Column &miller = t.at("miller_index");

  std::vector<double> offsets;
  for (std::size_t i = 0; i < t.nrows; ++i) {
    const Residual r = centroid_residual(
        truth, 0, static_cast<int>(miller.integer(i, 0)),
        static_cast<int>(miller.integer(i, 1)),
        static_cast<int>(miller.integer(i, 2)), xyz.real(i, 0), xyz.real(i, 1),
        xyz.real(i, 2));
    if (r.valid) offsets.push_back(std::hypot(r.dx, r.dy));
  }
  std::sort(offsets.begin(), offsets.end());
  check::is_true(offsets.size() > 1000, "enough reflections");
  check::close(offsets[offsets.size() / 2], 0.0, 1e-9, "the median must be exact");
  // The ninetieth percentile is 1.9e-9 px rather than zero: the iteration in
  // the forward map stops after three passes, and for reflections whose angle
  // is moderately sensitive that leaves a couple of nanopixels. A fourth pass
  // would remove it and buy nothing.
  check::close(offsets[offsets.size() * 9 / 10], 0.0, 1e-7,
               "and the ninetieth percentile to within the iteration limit");
  const double bad =
      static_cast<double>(std::count_if(offsets.begin(), offsets.end(),
                                        [](double d) { return d > 0.05; })) /
      static_cast<double>(offsets.size());
  // Measured at 4.3 per cent, and near-independent of how large the drift is,
  // which is what marks it as a property of those reflections rather than of
  // the amount of drift.
  check::is_true(bad < 0.08, "the ill-conditioned tail must stay small");
}

TEST(reflections_with_an_undetermined_rotation_angle_are_dropped) {
  // Waterman eqn (40) divides by the volume of the parallelepiped formed by
  // the rotation axis, the reciprocal lattice vector and the beam. Near the
  // rotation axis it goes to zero: the angle at which such a reflection
  // diffracts is arbitrarily sensitive to the model, and the Lorentz factor
  // has the same asymptote so its observed angular centroid is poor as well.
  // Fitting them puts noise into the rotation-angle residual that no model can
  // remove, which is why DIALS discards below 0.05 and why this now does.
  const Experiment truth = base_experiment();
  const Table t = observations_from(truth, 3.0);

  const auto run = [&](double cutoff) {
    ExperimentList list;
    list.experiments.push_back(truth);
    RefineOptions options;
    options.min_volume = cutoff;
    options.outlier_sigma = 0.0;
    options.macrocycles = 1;
    return refine(list, t, options);
  };

  const RefineResult all = run(0.0);
  const RefineResult cut = run(0.05);
  check::equal(static_cast<long long>(all.n_ill_conditioned), 0,
               "nothing dropped without a cutoff");
  check::is_true(cut.n_ill_conditioned > 0, "the cutoff drops something");
  check::is_true(cut.n_used < all.n_used, "and so fits fewer reflections");
  // The count is reported because the residual is averaged over what is left,
  // and a residual over a different set compares with nothing.
  check::equal(static_cast<long long>(all.n_used - cut.n_used),
               static_cast<long long>(cut.n_ill_conditioned),
               "the difference is exactly what was dropped");
}

TEST(the_cutoff_keeps_the_reflections_far_from_the_rotation_axis) {
  // It must remove a specific population, not simply the first few. Every
  // reflection kept has to have a volume above the cutoff and every one
  // dropped below it -- otherwise the filter is measuring something else.
  const Experiment truth = base_experiment();
  const Table t = observations_from(truth, 3.0);
  const Column &xyz = t.at("xyzobs.px.value");
  const Column &miller = t.at("miller_index");

  std::size_t below = 0, above = 0;
  for (std::size_t i = 0; i < t.nrows; ++i) {
    const PredictionState s = prediction_state(
        truth, 0, static_cast<int>(miller.integer(i, 0)),
        static_cast<int>(miller.integer(i, 1)),
        static_cast<int>(miller.integer(i, 2)), xyz.real(i, 2));
    if (!s.valid) continue;
    if (std::abs(s.volume) < 0.05) {
      ++below;
    } else {
      ++above;
    }
  }
  check::is_true(below > 0, "the synthetic data contains such reflections");
  check::is_true(above > below * 4, "and they are the minority");

  ExperimentList list;
  list.experiments.push_back(truth);
  RefineOptions options;
  options.min_volume = 0.05;
  options.outlier_sigma = 0.0;
  options.macrocycles = 1;
  const RefineResult result = refine(list, t, options);
  check::equal(static_cast<long long>(result.n_ill_conditioned),
               static_cast<long long>(below),
               "exactly the low-volume reflections are the ones dropped");
}

TEST(the_cell_and_the_detector_distance_are_degenerate_in_position) {
  // Scaling the cell and the detector distance together leaves the positions
  // on the detector almost unchanged. Measured here, against a truth that fits
  // exactly, both moved by two tenths of a per cent:
  //
  //     distance and cell together     0.062  0.072  0.318
  //     distance alone                 0.778  1.100  0.000
  //     cell alone                     0.836  1.152  0.318
  //
  // Twelve times smaller in position when they move together, which is the
  // degenerate direction. The rotation angle is NOT degenerate -- it sees the
  // cell and not the distance -- which is why a static refinement pins the
  // pair and a scan-varying one, with a hundred and sixty crystal parameters
  // free to absorb the angular residual by drifting the orientation, does not.
  //
  // On 1800 images of insulin at eighteen control points the distance drifted
  // 0.27 mm and the cell volume fell 0.59 per cent, for five thousandths of a
  // pixel. Holding the detector where the static pass put it kept the volume
  // within 0.04 per cent of what dials.refine reports.
  const Experiment truth = base_experiment();
  const Table t = observations_from(truth, 3.0);
  const Column &xyz = t.at("xyzobs.px.value");
  const Column &miller = t.at("miller_index");

  const auto rmsd = [&](double distance_scale, double cell_scale) {
    Experiment scaled = truth;
    const Vec3 normal =
        truth.detector[0].fast.cross(truth.detector[0].slow).normalized();
    const double distance = truth.detector[0].origin.dot(normal);
    scaled.detector.panels[0].origin =
        truth.detector[0].origin + normal * (distance * (distance_scale - 1.0));
    scaled.crystal->A = truth.crystal->A * (1.0 / cell_scale);

    double sx = 0.0, sz = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < t.nrows; ++i) {
      const Residual r = centroid_residual(
          scaled, 0, static_cast<int>(miller.integer(i, 0)),
          static_cast<int>(miller.integer(i, 1)),
          static_cast<int>(miller.integer(i, 2)), xyz.real(i, 0), xyz.real(i, 1),
          xyz.real(i, 2));
      if (!r.valid) continue;
      sx += r.dx * r.dx;
      sz += r.dz * r.dz;
      ++n;
    }
    check::is_true(n > 1000, "enough reflections");
    return std::make_pair(std::sqrt(sx / static_cast<double>(n)),
                          std::sqrt(sz / static_cast<double>(n)));
  };

  const auto exact = rmsd(1.0, 1.0);
  check::close(exact.first, 0.0, 1e-9, "the truth fits exactly");

  const auto together = rmsd(1.002, 1.002);
  const auto distance_only = rmsd(1.002, 1.0);

  check::is_true(together.first < 0.2 * distance_only.first,
                 "moving both together costs far less in position");
  // And the rotation angle does not join in: it responds to the cell and not
  // to the distance, which is what stops the pair drifting in a static fit.
  check::close(distance_only.second, 0.0, 1e-9,
               "the distance alone does not move the rotation angle");
  check::is_true(together.second > 0.1,
                 "but the cell does, so the angle breaks the degeneracy");
}
