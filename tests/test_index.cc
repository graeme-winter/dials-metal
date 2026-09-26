// Indexing, tested against ground truth.
//
// The closed loop is the point of having built prediction first: take a known
// crystal, predict where its reflections land, throw the indices away, and
// require the indexer to recover the cell from the centroids alone. Nothing
// about the answer is a judgement call -- the cell that went in is the cell
// that must come out.
//
// What the loop cannot test is anything that depends on the data being real:
// outliers, a second lattice, spots from ice. Those need real files, and the
// real-file result is recorded in README rather than asserted here, because
// the strong spot list is fifteen megabytes and does not belong in git.

#include <algorithm>
#include <cmath>
#include <map>
#include <tuple>
#include <complex>

#include "../src/expt.h"
#include "../src/fft.h"
#include "../src/index.h"
#include "../src/predict.h"
#include "../src/refl.h"
#include "check.h"

namespace mxi {

namespace {

Experiment synthetic(double cell = 78.0, double tilt = 0.7) {
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
  e.detector.panels.push_back(p);

  e.goniometer.axis = {1.0, 0.0, 0.0};
  e.scan.first_image = 1;
  e.scan.last_image = 900;
  e.scan.osc_start = 0.0;
  e.scan.osc_width = 0.1;

  // A general orientation, so nothing is accidentally axis-aligned.
  const Mat3 r = rotation({0.3, -0.5, 0.81}, tilt);
  e.crystal = Crystal::from_real_space(r * Vec3{cell, 0.0, 0.0},
                                       r * Vec3{0.0, cell, 0.0},
                                       r * Vec3{0.0, 0.0, cell});
  return e;
}

// Predict, then present the predictions as if they were found spots.
Table spots_from(const Experiment &e, double d_min) {
  PredictOptions options;
  options.d_min = d_min;
  const std::vector<Prediction> predictions = predict(e, options);

  Table t;
  t.nrows = predictions.size();
  Column &xyz = t.real_column("xyzobs.px.value", "vec3<double>", 3);
  Column &panel = t.int_column("panel", "std::size_t", 1);
  Column &id = t.int_column("id", "int", 1);
  for (std::size_t i = 0; i < predictions.size(); ++i) {
    xyz.reals[i * 3 + 0] = predictions[i].px_fast;
    xyz.reals[i * 3 + 1] = predictions[i].px_slow;
    xyz.reals[i * 3 + 2] = predictions[i].z;
    panel.ints[i] = static_cast<std::int64_t>(predictions[i].panel);
    id.ints[i] = 0;
  }
  return t;
}

// Two cells describe the same lattice if their metric tensors match after some
// unimodular change of basis. Comparing cell parameters directly would fail on
// a perfectly good answer in a different setting, so compare what cannot
// change: the volume, and the sorted edge lengths after reduction.
void check_same_cell(const UnitCell &got, const UnitCell &want,
                     double tolerance) {
  double a[3] = {got.a, got.b, got.c};
  double b[3] = {want.a, want.b, want.c};
  std::sort(a, a + 3);
  std::sort(b, b + 3);
  for (int i = 0; i < 3; ++i) {
    check::close(a[i], b[i], tolerance * b[i], "cell edge");
  }
  check::close(got.volume(), want.volume(), 3.0 * tolerance * want.volume(),
               "cell volume");
}

}  // namespace

TEST(index_recovers_a_cubic_cell_from_its_own_predictions) {
  Experiment truth = synthetic(78.0);
  Table spots = spots_from(truth, 3.0);
  check::is_true(spots.nrows > 2000, "the sweep should predict plenty");

  ExperimentList list;
  Experiment blank = truth;
  blank.crystal.reset();
  list.experiments.push_back(blank);

  IndexOptions options;
  const IndexResult result = index(list, spots, options);

  check::is_true(result.fraction_indexed() > 0.95,
                 "nearly every predicted spot should index");
  check_same_cell(result.crystal.cell(), truth.crystal->cell(), 0.005);
  // Predictions are exact, so the fit should be far better than the acceptance
  // tolerance. A large residual here means the basis is a near miss.
  check::is_true(result.rmsd_index < 0.02, "residual on exact data");
  check::is_true(list[0].crystal.has_value(), "the crystal is set on the list");
}

TEST(index_recovers_a_triclinic_cell) {
  // Cubic is the easy case: every candidate vector is the same length and
  // almost any three of them work. A cell with three different edges and no
  // right angles is what separates a working basis search from a lucky one.
  Experiment truth = synthetic();
  const Mat3 r = rotation({0.1, 0.9, -0.42}, 1.1);
  truth.crystal = Crystal::from_real_space(r * Vec3{41.0, 0.0, 0.0},
                                           r * Vec3{7.0, 52.0, 0.0},
                                           r * Vec3{-5.0, 9.0, 63.0});
  Table spots = spots_from(truth, 2.5);
  check::is_true(spots.nrows > 2000, "enough spots");

  ExperimentList list;
  Experiment blank = truth;
  blank.crystal.reset();
  list.experiments.push_back(blank);

  const IndexResult result = index(list, spots, {});
  check::is_true(result.fraction_indexed() > 0.9, "most spots should index");
  check_same_cell(result.crystal.cell(), truth.crystal->cell(), 0.01);
}

TEST(index_pools_reciprocal_space_across_sweeps) {
  // Two sweeps of one crystal at different goniometer settings, indexed
  // together. Each keeps its own goniometer, so the pooling into one crystal
  // frame is exact rather than an approximation -- which is the whole reason
  // joint indexing is the right first step for a multi-sweep experiment.
  Experiment first = synthetic(78.0);
  Experiment second = synthetic(78.0);
  second.goniometer.fixed = rotation({0.0, 1.0, 0.0}, Scan::radians(35.0));

  Table a = spots_from(first, 3.5);
  Table b = spots_from(second, 3.5);
  check::is_true(a.nrows > 500 && b.nrows > 500, "both sweeps give spots");

  Table both;
  both.nrows = a.nrows + b.nrows;
  Column &xyz = both.real_column("xyzobs.px.value", "vec3<double>", 3);
  Column &panel = both.int_column("panel", "std::size_t", 1);
  Column &id = both.int_column("id", "int", 1);
  for (std::size_t i = 0; i < a.nrows; ++i) {
    for (std::size_t k = 0; k < 3; ++k) {
      xyz.reals[i * 3 + k] = a.at("xyzobs.px.value").real(i, k);
    }
    panel.ints[i] = 0;
    id.ints[i] = 0;
  }
  for (std::size_t i = 0; i < b.nrows; ++i) {
    const std::size_t row = a.nrows + i;
    for (std::size_t k = 0; k < 3; ++k) {
      xyz.reals[row * 3 + k] = b.at("xyzobs.px.value").real(i, k);
    }
    panel.ints[row] = 0;
    id.ints[row] = 1;
  }

  ExperimentList list;
  Experiment x = first, y = second;
  x.crystal.reset();
  y.crystal.reset();
  list.experiments.push_back(x);
  list.experiments.push_back(y);

  const IndexResult result = index(list, both, {});
  check::is_true(result.fraction_indexed() > 0.9,
                 "both sweeps must index against one crystal");
  check_same_cell(result.crystal.cell(), first.crystal->cell(), 0.01);
  // Every experiment gets the crystal, which is what makes the basis common.
  check::is_true(list[0].crystal.has_value() && list[1].crystal.has_value(),
                 "the crystal is shared across the list");
}

TEST(index_writes_miller_indices_that_agree_with_the_truth) {
  // Up to a reindexing: the indexer has no way to prefer one basis of a cubic
  // lattice over another, so the recovered indices are the true ones
  // transformed by some unimodular matrix. What must hold is that the
  // transform is the SAME for every reflection, which is what makes the
  // solution a lattice rather than a coincidence.
  Experiment truth = synthetic(78.0);
  PredictOptions po;
  po.d_min = 3.5;
  const std::vector<Prediction> predictions = predict(truth, po);

  Table spots = spots_from(truth, 3.5);
  ExperimentList list;
  Experiment blank = truth;
  blank.crystal.reset();
  list.experiments.push_back(blank);
  const IndexResult result = index(list, spots, {});
  check::is_true(result.n_indexed > 0, "something indexed");

  // h_true = T h_found for a single integer T, recovered from three
  // independent reflections and then checked against all of them.
  const Column &miller = spots.at("miller_index");
  std::vector<std::size_t> rows;
  for (std::size_t i = 0; i < spots.nrows && rows.size() < 3; ++i) {
    const Vec3 h{static_cast<double>(miller.integer(i, 0)),
                 static_cast<double>(miller.integer(i, 1)),
                 static_cast<double>(miller.integer(i, 2))};
    if (h.norm_squared() == 0.0) continue;
    if (rows.size() == 2) {
      const Vec3 p{static_cast<double>(miller.integer(rows[0], 0)),
                   static_cast<double>(miller.integer(rows[0], 1)),
                   static_cast<double>(miller.integer(rows[0], 2))};
      const Vec3 q{static_cast<double>(miller.integer(rows[1], 0)),
                   static_cast<double>(miller.integer(rows[1], 1)),
                   static_cast<double>(miller.integer(rows[1], 2))};
      if (std::abs(Mat3::from_rows(p, q, h).determinant()) < 0.5) continue;
    }
    rows.push_back(i);
  }
  check::equal(static_cast<long long>(rows.size()), 3,
               "three independent indexed reflections");

  Mat3 found = Mat3::from_rows(
      {static_cast<double>(miller.integer(rows[0], 0)),
       static_cast<double>(miller.integer(rows[0], 1)),
       static_cast<double>(miller.integer(rows[0], 2))},
      {static_cast<double>(miller.integer(rows[1], 0)),
       static_cast<double>(miller.integer(rows[1], 1)),
       static_cast<double>(miller.integer(rows[1], 2))},
      {static_cast<double>(miller.integer(rows[2], 0)),
       static_cast<double>(miller.integer(rows[2], 1)),
       static_cast<double>(miller.integer(rows[2], 2))});
  Mat3 wanted = Mat3::from_rows(
      {static_cast<double>(predictions[rows[0]].h),
       static_cast<double>(predictions[rows[0]].k),
       static_cast<double>(predictions[rows[0]].l)},
      {static_cast<double>(predictions[rows[1]].h),
       static_cast<double>(predictions[rows[1]].k),
       static_cast<double>(predictions[rows[1]].l)},
      {static_cast<double>(predictions[rows[2]].h),
       static_cast<double>(predictions[rows[2]].k),
       static_cast<double>(predictions[rows[2]].l)});

  bool ok = false;
  const Mat3 transform = found.inverse(&ok) * wanted;
  check::is_true(ok, "the three reflections must be independent");
  // A change of basis between two descriptions of one lattice is unimodular.
  check::close(std::abs(transform.determinant()), 1.0, 1e-6,
               "the transform must be unimodular");

  std::size_t agree = 0, tested = 0;
  for (std::size_t i = 0; i < spots.nrows; ++i) {
    const Vec3 h{static_cast<double>(miller.integer(i, 0)),
                 static_cast<double>(miller.integer(i, 1)),
                 static_cast<double>(miller.integer(i, 2))};
    if (h.norm_squared() == 0.0) continue;
    ++tested;
    const Vec3 mapped = Mat3::from_rows(h, {0, 0, 0}, {0, 0, 0}) * Vec3{0, 0, 0};
    (void)mapped;
    const Vec3 t{h.dot(transform.column(0)), h.dot(transform.column(1)),
                 h.dot(transform.column(2))};
    if (std::abs(t.x - predictions[i].h) < 1e-6 &&
        std::abs(t.y - predictions[i].k) < 1e-6 &&
        std::abs(t.z - predictions[i].l) < 1e-6) {
      ++agree;
    }
  }
  check::is_true(tested > 500, "plenty to check");
  check::is_true(static_cast<double>(agree) / static_cast<double>(tested) > 0.99,
                 "one transform must explain nearly every index");
}

TEST(reduce_basis_shortens_without_changing_the_lattice) {
  const Mat3 rows = Mat3::from_rows({10.0, 0.0, 0.0}, {90.0, 10.0, 0.0},
                                    {30.0, 70.0, 10.0});
  const Mat3 reduced = reduce_basis(rows);
  // Volume is the invariant: reduction is a unimodular change of basis.
  check::close(std::abs(reduced.determinant()), std::abs(rows.determinant()),
               1e-9 * std::abs(rows.determinant()), "volume preserved");
  double before = 0.0, after = 0.0;
  for (std::size_t i = 0; i < 3; ++i) {
    before += rows.row(i).norm_squared();
    after += reduced.row(i).norm_squared();
  }
  check::is_true(after < before, "the reduced basis must be shorter");
  check::is_true(reduced.determinant() > 0.0, "and right-handed");
}

TEST(reduce_basis_leaves_an_already_reduced_cell_alone) {
  const Mat3 rows =
      Mat3::from_rows({10.0, 0.0, 0.0}, {0.0, 11.0, 0.0}, {0.0, 0.0, 12.0});
  const Mat3 reduced = reduce_basis(rows);
  double length[3];
  for (std::size_t i = 0; i < 3; ++i) length[i] = reduced.row(i).norm();
  std::sort(length, length + 3);
  check::close(length[0], 10.0, 1e-9, "a");
  check::close(length[1], 11.0, 1e-9, "b");
  check::close(length[2], 12.0, 1e-9, "c");
}

TEST(estimate_max_cell_recovers_a_known_spacing) {
  // A perfect reciprocal lattice of spacing 1/60, so the estimate should come
  // out near 60 -- above it, since the estimate carries a deliberate margin
  // and is a lower bound on the cell before that.
  std::vector<Vec3> points;
  const double spacing = 1.0 / 60.0;
  for (int h = -5; h <= 5; ++h) {
    for (int k = -5; k <= 5; ++k) {
      for (int l = -5; l <= 5; ++l) {
        points.push_back({h * spacing, k * spacing, l * spacing});
      }
    }
  }
  const double estimate = estimate_max_cell(points);
  check::is_true(estimate >= 60.0, "must not underestimate the cell");
  check::is_true(estimate < 120.0, "and must not be wild");
}

TEST(index_refuses_rather_than_inventing_a_lattice) {
  // Random points are not a lattice. An indexer that returns a cell for them
  // will return a cell for anything.
  Experiment e = synthetic();
  Table spots;
  spots.nrows = 400;
  Column &xyz = spots.real_column("xyzobs.px.value", "vec3<double>", 3);
  Column &panel = spots.int_column("panel", "std::size_t", 1);
  Column &id = spots.int_column("id", "int", 1);
  unsigned state = 12345;
  const auto next = [&state]() {
    state = state * 1664525u + 1013904223u;
    return static_cast<double>(state >> 8) / static_cast<double>(1u << 24);
  };
  for (std::size_t i = 0; i < spots.nrows; ++i) {
    xyz.reals[i * 3 + 0] = 100.0 + next() * 3900.0;
    xyz.reals[i * 3 + 1] = 100.0 + next() * 4100.0;
    xyz.reals[i * 3 + 2] = next() * 900.0;
    panel.ints[i] = 0;
    id.ints[i] = 0;
  }
  ExperimentList list;
  Experiment blank = e;
  blank.crystal.reset();
  list.experiments.push_back(blank);

  const IndexResult result = index(list, spots, {});
  check::is_true(result.fraction_indexed() < 0.5,
                 "random points must not index well");
}

TEST(macrocycles_keep_a_weak_population_from_dragging_the_model) {
  // The mechanism this exists for, reproduced. Take exact predictions from a
  // known crystal, then add a population of weak spots displaced by a fifth of
  // a pixel outwards -- small enough to index within tolerance, consistent
  // enough that no outlier rejection will find them, which is precisely the
  // situation on real data.
  //
  // Assigning once and fitting everything lets them pull the cell. Refining on
  // the strong half and re-assigning does not.
  Experiment truth = synthetic(78.0);
  Table good = spots_from(truth, 3.0);
  check::is_true(good.nrows > 2000, "enough good spots");

  const std::size_t n_bad = good.nrows / 4;
  Table all;
  all.nrows = good.nrows + n_bad;
  Column &xyz = all.real_column("xyzobs.px.value", "vec3<double>", 3);
  Column &panel = all.int_column("panel", "std::size_t", 1);
  Column &id = all.int_column("id", "int", 1);
  Column &signal = all.int_column("n_signal", "int", 1);
  const Column &from = good.at("xyzobs.px.value");

  const Vec3 centre{2074.0, 2181.0, 0.0};
  for (std::size_t i = 0; i < good.nrows; ++i) {
    for (std::size_t k = 0; k < 3; ++k) xyz.reals[i * 3 + k] = from.real(i, k);
    signal.ints[i] = 30;  // strong
  }
  for (std::size_t i = 0; i < n_bad; ++i) {
    const std::size_t row = good.nrows + i;
    const std::size_t src = (i * 3) % good.nrows;
    const double dx = from.real(src, 0) - centre.x;
    const double dy = from.real(src, 1) - centre.y;
    const double r = std::hypot(dx, dy);
    // Pushed outwards: a coherent bias, not noise.
    xyz.reals[row * 3 + 0] = from.real(src, 0) + 0.2 * dx / std::fmax(r, 1.0);
    xyz.reals[row * 3 + 1] = from.real(src, 1) + 0.2 * dy / std::fmax(r, 1.0);
    xyz.reals[row * 3 + 2] = from.real(src, 2);
    signal.ints[row] = 3;  // weak
  }
  (void)panel;
  (void)id;

  const auto recover = [&](int cycles, bool strong) {
    ExperimentList list;
    Experiment blank = truth;
    blank.crystal.reset();
    list.experiments.push_back(blank);
    Table copy = all;
    IndexOptions options;
    options.macrocycles = cycles;
    options.refine_on_strong = strong;
    const IndexResult r = index(list, copy, options);
    const UnitCell c = r.crystal.cell();
    const UnitCell t = truth.crystal->cell();
    return std::abs(c.volume() - t.volume()) / t.volume();
  };

  const double once = recover(0, false);
  const double cycled = recover(3, true);
  check::is_true(cycled < once,
                 "macrocycles on strong reflections must beat assigning once");
  check::is_true(cycled < 0.002, "and land within 0.2 per cent on volume");
}

TEST(the_transform_in_use_agrees_with_the_built_in) {
  // Says nothing unless a library is linked, and then says the thing that
  // matters: a library can be planned wrongly, normalised differently, or
  // linked against a build of another precision, and none of that shows in a
  // result that still looks like a lattice. Peaks in roughly the right places
  // are not evidence. Agreeing with a transform that is itself tested against
  // a direct summation is.
  for (std::size_t n : {std::size_t(8), std::size_t(16)}) {
    std::vector<std::complex<double>> a(n * n * n);
    std::uint64_t state = 987654321;
    const auto uniform = [&state]() {
      state = state * 6364136223846793005ULL + 1442695040888963407ULL;
      return static_cast<double>((state >> 11) & ((1ULL << 53) - 1)) /
                 static_cast<double>(1ULL << 53) -
             0.5;
    };
    for (auto &z : a) z = {uniform(), uniform()};
    std::vector<std::complex<double>> b = a;

    fft3d(a, n, +1);
    fft3d_builtin(b, n, +1);

    double worst = 0.0;
    double scale = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
      worst = std::fmax(worst, std::abs(a[i] - b[i]));
      scale = std::fmax(scale, std::abs(b[i]));
    }
    check::is_true(scale > 0.0, "the reference is not all zero");
    check::is_true(worst < 1e-9 * scale,
                   "the transform in use agrees with the built-in");
  }
}

TEST(max_cell_is_not_fooled_by_the_same_reflection_measured_again) {
  // A sweep of several full rotations measures every reflection once per turn,
  // and rotating s1 back by phi and by phi + 360 degrees gives the same vector
  // -- so each point has copies of itself, separated only by the error in its
  // centroid. Ungrouped, the nearest-neighbour spacing measures that error
  // instead of the lattice.
  //
  // The numbers below are the failure, not a hypothetical: a hundred-Angstrom
  // cubic lattice over ten turns with a 1e-4 scatter estimates at eleven
  // thousand Angstroms, which finds no candidate vectors at all.
  std::vector<Vec3> one;
  const double d = 1.0 / 67.0;
  for (int h = -8; h <= 8; ++h) {
    for (int k = -8; k <= 8; ++k) {
      for (int l = -8; l <= 8; ++l) {
        if (!h && !k && !l) continue;
        one.push_back({h * d, k * d, l * d});
      }
    }
  }
  const double truth = estimate_max_cell(one);
  check::close(truth, 100.5, 1.0, "one turn gives the right cell");

  // A deterministic scatter, so the test cannot pass or fail by luck.
  std::uint64_t state = 4321;
  const auto jitter = [&state](double scale) {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    const double u = static_cast<double>((state >> 11) & ((1ULL << 53) - 1)) /
                     static_cast<double>(1ULL << 53);
    return (u - 0.5) * 2.0 * scale;
  };

  for (double scale : {1e-4, 1e-3}) {
    std::vector<Vec3> many;
    std::vector<int> groups;
    for (int turn = 0; turn < 10; ++turn) {
      for (const Vec3 &p : one) {
        many.push_back({p.x + jitter(scale), p.y + jitter(scale),
                        p.z + jitter(scale)});
        groups.push_back(turn);
      }
    }
    const double ungrouped = estimate_max_cell(many);
    const double grouped = estimate_max_cell(many, groups);
    // Ungrouped is wrong here only when the copies land outside the 1e-4 floor
    // the estimator already has. With this scatter some fall inside it and are
    // skipped, which is why the floor looked like it worked for so long: it
    // rescues some datasets and not others, and nothing says which.
    check::is_true(ungrouped >= grouped,
                   "grouping never makes the estimate larger");
    check::is_true(grouped < 1.5 * truth,
                   "grouped stays close to the one-turn answer");
  }
}

TEST(a_sweep_shorter_than_a_block_is_all_one_group) {
  // The grouping must not split a sweep smaller than one block: doing so would
  // compare each point with a fraction of the lattice and estimate a cell that
  // is too large, which is the bug it was written to fix.
  Experiment e;
  e.beam.direction = {0.0, 0.0, 1.0};
  e.beam.wavelength = 1.0;
  e.goniometer.axis = {1.0, 0.0, 0.0};
  e.scan.first_image = 1;
  e.scan.last_image = 20;
  e.scan.osc_start = 0.0;
  e.scan.osc_width = 1.0;  // 20 degrees in total, inside one 30 degree block
  ExperimentList list;
  list.experiments.push_back(e);

  Table t;
  t.nrows = 5;
  Column &obs = t.real_column("xyzobs.px.value", "vec3<double>", 3);
  for (std::size_t i = 0; i < t.nrows; ++i) {
    obs.reals[i * 3 + 2] = static_cast<double>(i) * 4.0;  // 0 to 16 degrees
  }
  const std::vector<int> groups = observation_groups(list, t);
  for (std::size_t i = 1; i < groups.size(); ++i) {
    check::equal(static_cast<long long>(groups[i]),
                 static_cast<long long>(groups[0]),
                 "a sweep inside one block is one group");
  }

  // And a full rotation is twelve blocks of thirty degrees, so the two ends of
  // it are not compared with each other. That is the point: a reflection's two
  // passages through the Ewald sphere are far apart in phi and map to the same
  // place in the crystal frame.
  e.scan.last_image = 360;
  list.experiments[0] = e;
  for (std::size_t i = 0; i < t.nrows; ++i) {
    obs.reals[i * 3 + 2] = static_cast<double>(i) * 70.0;  // 0 to 280 degrees
  }
  const std::vector<int> spread = observation_groups(list, t);
  check::is_true(spread.front() != spread.back(),
                 "the two ends of a rotation are different blocks");
  check::equal(static_cast<long long>(spread[0]),
               static_cast<long long>(observation_groups(list, t, 360.0)[0]),
               "and a 360 degree block puts them together again");
}

TEST(a_reflection_is_seen_twice_in_one_rotation_and_the_blocks_know_it) {
  // The failure that grouping by turn did not fix. Over a full rotation every
  // reciprocal lattice point crosses the Ewald sphere twice, entering and
  // leaving, at two values of phi far apart -- and both observations rotate
  // back to the same place in the crystal frame. One turn therefore holds two
  // copies of everything, and grouping by turn keeps them together.
  //
  // Ten turns with two passages each, on a 90 Angstrom lattice: 6517 Angstroms
  // ungrouped, 2244 grouped by turn, 141 grouped by thirty degrees against a
  // true 135. Those are the numbers this reproduces.
  std::vector<Vec3> one;
  const double d = 1.0 / 90.0;
  for (int h = -10; h <= 10; ++h) {
    for (int k = -10; k <= 10; ++k) {
      for (int l = -10; l <= 10; ++l) {
        if (!h && !k && !l) continue;
        one.push_back({h * d, k * d, l * d});
      }
    }
  }
  const double truth = estimate_max_cell(one);

  std::uint64_t state = 99;
  const auto jitter = [&state](double scale) {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    const double u = static_cast<double>((state >> 11) & ((1ULL << 53) - 1)) /
                     static_cast<double>(1ULL << 53);
    return (u - 0.5) * 2.0 * scale;
  };

  std::vector<Vec3> many;
  std::vector<int> by_turn;
  std::vector<int> by_block;
  for (int turn = 0; turn < 10; ++turn) {
    // The two passages of one reflection, a hundred and fifty degrees apart.
    for (double phi : {20.0, 170.0}) {
      for (const Vec3 &p : one) {
        many.push_back({p.x + jitter(5e-4), p.y + jitter(5e-4),
                        p.z + jitter(5e-4)});
        const double angle = turn * 360.0 + phi;
        by_turn.push_back(static_cast<int>(angle / 360.0));
        by_block.push_back(static_cast<int>(angle / 30.0));
      }
    }
  }

  const double by_turn_estimate = estimate_max_cell(many, by_turn);
  const double by_block_estimate = estimate_max_cell(many, by_block);
  check::is_true(by_turn_estimate > 5.0 * truth,
                 "grouping by turn is not enough, which is the whole point");
  check::is_true(by_block_estimate < 1.2 * truth,
                 "grouping by thirty degrees recovers the cell");
}

TEST(a_sweep_of_several_turns_predicts_each_reflection_once_per_turn) {
  // A reciprocal lattice point that crosses the Ewald sphere at phi crosses it
  // again at phi + 2 pi. A sweep of ten full rotations therefore records every
  // reflection ten times, which is the entire reason for collecting one.
  //
  // Predicting it once left a tenth of the reflections: on a real ten rotation
  // sweep, 749786 predictions where DIALS made 7516507. The ratio was 10.02,
  // which is the number of turns and is what gave it away.
  const auto sweep = [](int turns) {
    Experiment e;
    e.beam.direction = {0.0, 0.0, 1.0};
    e.beam.wavelength = 0.9537;
    Panel p;
    p.fast = {1.0, 0.0, 0.0};
    p.slow = {0.0, -1.0, 0.0};
    p.pixel_size[0] = p.pixel_size[1] = 0.075;
    p.image_size[0] = 2068;
    p.image_size[1] = 2162;
    p.origin = {-77.5, 81.1, -168.5};
    e.detector.panels.push_back(p);
    e.goniometer.axis = {1.0, 0.0, 0.0};
    e.scan.first_image = 1;
    e.scan.last_image = 3600 * turns;
    e.scan.osc_start = 0.0;
    e.scan.osc_width = 0.1;
    const Mat3 r = rotation({0.3, -0.5, 0.81}, 0.7);
    e.crystal = Crystal::from_real_space(r * Vec3{78, 0, 0}, r * Vec3{0, 78, 0},
                                         r * Vec3{0, 0, 78});
    return e;
  };

  PredictOptions options;
  options.d_min = 3.0;
  const std::size_t one = predict(sweep(1), options).size();
  check::is_true(one > 1000, "a single turn predicts something");

  for (int turns : {2, 5, 10}) {
    const Experiment e = sweep(turns);
    const std::vector<Prediction> v = predict(e, options);
    check::equal(static_cast<long long>(v.size()),
                 static_cast<long long>(one) * turns,
                 "exactly one prediction per reflection per turn");
    // And every one of them is on the scan. A full turn or more used to skip
    // the wrap entirely, on the grounds that everything is inside it --
    // everything is, but only after wrapping, and without it z_from_phi put
    // 378278 of 749786 predictions off the end of the scan.
    const double images = static_cast<double>(e.scan.num_images());
    std::size_t off = 0;
    for (const Prediction &p : v) {
      if (p.z < -1.0 || p.z > images + 1.0) ++off;
    }
    check::equal(static_cast<long long>(off), 0, "and none is off the scan");
  }
}

TEST(a_scan_varying_crystal_also_predicts_once_per_turn) {
  // The test above uses a scan-static crystal and passed while this was
  // broken. converge_root answers with a crossing from ewald_intersections,
  // which is in a principal 2 pi interval whatever it was seeded with, so the
  // turn was lost on the way out -- and build then wrapped phi into the first
  // turn anyway, because wrap_from returns a value in [lo, lo + 2 pi).
  //
  // On ten rotations of real data every shoebox landed in the first 3627
  // frames of 36000 and the same pixels were integrated ten times over. What
  // caught it was a diagnostic, not a test: "frames read (N wanted by a
  // shoebox, each read 7.55 times)".
  //
  // So this asserts what the other test asserts, on a crystal that varies --
  // which is every real one.
  const auto sweep = [](int turns) {
    Experiment e;
    e.beam.direction = {0.0, 0.0, 1.0};
    e.beam.wavelength = 0.9537;
    Panel p;
    p.fast = {1.0, 0.0, 0.0};
    p.slow = {0.0, -1.0, 0.0};
    p.pixel_size[0] = p.pixel_size[1] = 0.075;
    p.image_size[0] = 2068;
    p.image_size[1] = 2162;
    p.origin = {-77.5, 81.1, -168.5};
    e.detector.panels.push_back(p);
    e.goniometer.axis = {1.0, 0.0, 0.0};
    e.scan.first_image = 1;
    e.scan.last_image = 3600 * turns;
    e.scan.osc_start = 0.0;
    e.scan.osc_width = 0.1;
    const Mat3 r = rotation({0.3, -0.5, 0.81}, 0.7);
    e.crystal = Crystal::from_real_space(r * Vec3{78, 0, 0}, r * Vec3{0, 78, 0},
                                         r * Vec3{0, 0, 78});
    // A scan-varying crystal: the same A at every point, so the geometry is
    // identical to the static case and only the code path differs. Any
    // difference in the answer is then the code path and nothing else.
    e.crystal->A_points.assign(static_cast<std::size_t>(3600 * turns) + 1,
                               e.crystal->A);
    check::is_true(e.crystal->scan_varying(), "the crystal really does vary");
    return e;
  };

  PredictOptions options;
  options.d_min = 3.0;
  const std::size_t one = predict(sweep(1), options).size();
  check::is_true(one > 1000, "a single turn predicts something");

  for (int turns : {2, 5, 10}) {
    const Experiment e = sweep(turns);
    const std::vector<Prediction> v = predict(e, options);
    check::equal(static_cast<long long>(v.size()),
                 static_cast<long long>(one) * turns,
                 "one prediction per reflection per turn, scan-varying");

    // And they are spread over the whole scan rather than piled into its
    // first turn, which is the failure this is here for. Every turn should
    // hold about the same number.
    std::vector<std::size_t> per_turn(static_cast<std::size_t>(turns), 0);
    for (const Prediction &p : v) {
      const int turn = static_cast<int>(p.z / 3600.0);
      if (turn >= 0 && turn < turns) ++per_turn[static_cast<std::size_t>(turn)];
    }
    for (std::size_t t = 0; t < per_turn.size(); ++t) {
      check::is_true(per_turn[t] > one / 2,
                     "every turn has its own share of the predictions");
    }
  }
}

TEST(a_crystal_that_moves_does_not_collapse_its_turns_together) {
  // The test above uses the same A at every scan point and passed while this
  // was broken, because with an A that does not change the geometry is the
  // scan-static one and only the code path differs. Every real crystal moves.
  //
  // converge_root compared an angle in turn k against two roots in turn zero,
  // both about 2 pi k away, so which was nearer depended on which side of zero
  // each fell -- and it could take the OTHER crossing. Once a pass had
  // assigned a principal angle, the next evaluated the setting matrix at turn
  // zero's frame. On a four turn sweep that put 218640 of 583040 predictions
  // within five frames of another observation of the same reflection with the
  // same flag: the right count, in the wrong places.
  //
  // Found by mxeq explain on a real ten rotation comparison, which said that
  // 85.7 per cent of our unpartnered observations had a same-index,
  // same-flag partner within the radius -- the signature of duplicates.
  const int turns = 4;
  Experiment e;
  e.beam.direction = {0.0, 0.0, 1.0};
  e.beam.wavelength = 0.9537;
  Panel p;
  p.fast = {1.0, 0.0, 0.0};
  p.slow = {0.0, -1.0, 0.0};
  p.pixel_size[0] = p.pixel_size[1] = 0.075;
  p.image_size[0] = 2068;
  p.image_size[1] = 2162;
  p.origin = {-77.5, 81.1, -168.5};
  e.detector.panels.push_back(p);
  e.goniometer.axis = {1.0, 0.0, 0.0};
  e.scan.first_image = 1;
  e.scan.last_image = 3600 * turns;
  e.scan.osc_start = 0.0;
  e.scan.osc_width = 0.1;
  const Mat3 r = rotation({0.3, -0.5, 0.81}, 0.7);
  e.crystal = Crystal::from_real_space(r * Vec3{78, 0, 0}, r * Vec3{0, 78, 0},
                                       r * Vec3{0, 0, 78});
  // A crystal that genuinely moves through the scan.
  const std::size_t points = static_cast<std::size_t>(3600 * turns) + 1;
  e.crystal->A_points.resize(points);
  for (std::size_t k = 0; k < points; ++k) {
    const Mat3 wobble =
        rotation({0.0, 0.0, 1.0}, 1e-4 * std::sin(static_cast<double>(k) * 1e-3));
    e.crystal->A_points[k] = wobble * e.crystal->A;
  }

  PredictOptions options;
  options.d_min = 3.0;
  const std::vector<Prediction> v = predict(e, options);

  // No two observations of one reflection with the same flag within five
  // frames: those are a whole turn apart, or they are the same observation
  // predicted twice.
  std::map<std::tuple<int, int, int, bool>, std::vector<double>> groups;
  for (const Prediction &q : v) groups[{q.h, q.k, q.l, q.entering}].push_back(q.z);
  std::size_t collapsed = 0;
  for (auto &entry : groups) {
    std::vector<double> &z = entry.second;
    std::sort(z.begin(), z.end());
    for (std::size_t i = 1; i < z.size(); ++i) {
      if (z[i] - z[i - 1] < 5.0) ++collapsed;
    }
  }
  check::equal(static_cast<long long>(collapsed), 0,
               "no turn collapses onto another");
}

TEST(every_prediction_is_on_the_sphere_under_the_crystal_at_its_own_position) {
  // The invariant a scan-varying prediction has to satisfy, whatever else it
  // gets right: rotate A(z) h to the predicted phi, where z is the predicted
  // image, and it lands on the Ewald sphere. A prediction made with the crystal
  // from somewhere else in the scan does not.
  //
  // ewald_intersections returned a leaving root of -259.9 degrees for a
  // reflection at 100.1 on a scan from 0 to 180, and converge_root looked the
  // crystal up at z_from_phi of that: image -2600, clamped to image zero. So
  // it used the crystal from the start of the scan for a reflection in the
  // middle. On a real 1800 image sweep, leaving reflections were exact to
  // image 900 and up to 0.96 images late after it, and had been since the
  // first scan-varying version. Every earlier test passed, because none used
  // a crystal that moved enough for its start to differ from its middle.
  Experiment e;
  e.beam.direction = {0.0, 0.0, 1.0};
  e.beam.wavelength = 0.9537;
  Panel p;
  p.fast = {1.0, 0.0, 0.0};
  p.slow = {0.0, -1.0, 0.0};
  p.pixel_size[0] = p.pixel_size[1] = 0.075;
  p.image_size[0] = 2068;
  p.image_size[1] = 2162;
  p.origin = {-77.5, 81.1, -168.5};
  e.detector.panels.push_back(p);
  e.goniometer.axis = {1.0, 0.0, 0.0};
  e.scan.first_image = 1;
  e.scan.last_image = 1800;
  e.scan.osc_start = 0.0;
  e.scan.osc_width = 0.1;
  const Mat3 r = rotation({0.3, -0.5, 0.81}, 0.7);
  e.crystal = Crystal::from_real_space(r * Vec3{78, 0, 0}, r * Vec3{0, 78, 0},
                                       r * Vec3{0, 0, 78});
  // A crystal that turns steadily through the scan, half a degree end to end,
  // so that its model at image zero is not its model at image 900.
  const std::size_t points = 1801;
  e.crystal->A_points.resize(points);
  for (std::size_t k = 0; k < points; ++k) {
    const double turn = 0.5 * (M_PI / 180.0) * static_cast<double>(k) / 1800.0;
    e.crystal->A_points[k] = rotation({0.0, 1.0, 0.0}, turn) * e.crystal->A;
  }

  PredictOptions options;
  options.d_min = 2.5;
  const std::vector<Prediction> v = predict(e, options);
  check::is_true(v.size() > 1000, "something was predicted");

  const Vec3 s0 = e.beam.s0();
  const double radius = s0.norm();
  std::size_t off = 0;
  double worst = 0.0;
  for (const Prediction &q : v) {
    const Vec3 hkl{static_cast<double>(q.h), static_cast<double>(q.k),
                   static_cast<double>(q.l)};
    const Vec3 s1 = s0 + e.goniometer.rotation_at(q.phi) * (e.setting_at(q.z) * hkl);
    const double miss = std::fabs(s1.norm() - radius) / radius;
    worst = std::max(worst, miss);
    if (miss > 1e-5) ++off;
  }
  check::equal(static_cast<long long>(off), 0,
               "every prediction is on the sphere under its own crystal");
}

}  // namespace mxi
