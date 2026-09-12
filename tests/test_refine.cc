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
