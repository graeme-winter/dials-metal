// The goniometer decomposition, on the only data that can test it.
//
// Every setting angle in the insulin sweep is zero, so its fixed and setting
// rotations both collapse to the identity and any arrangement of them passes.
// That was recorded as an open item; this closes it.
//
// l-cysteine, four sweeps from one crystal on a fixed-chi goniometer. What
// makes each variant of getting it wrong detectable, and by which sweep, is
// the point of the tests at the bottom of this file. They are not decoration:
// a test that passes under the wrong composition as readily as the right one
// tests nothing, and the only way to know which this is, is to try the wrong
// ones.

#include <algorithm>
#include <cmath>
#include <vector>

#include "../src/geometry.hh"
#include "check.hh"
#include "real_cysteine.hh"
#include "real_threeaxis.hh"

namespace mxi {

namespace {

Vec3 v(const std::array<double, 3> &a) { return {a[0], a[1], a[2]}; }
Vec3 v(const double (&a)[3]) { return {a[0], a[1], a[2]}; }

// Both datasets are emitted by the same generator, so one template builds an
// Experiment from either. The axis count is whatever the file had -- fixing it
// at two would have silently truncated the three-axis goniometer.
template <typename Sweep>
Experiment build(const Sweep &s, const double (&a)[3], const double (&b)[3],
                 const double (&c)[3]) {
  Experiment e;
  e.beam.direction = v(s.beam_direction);
  e.beam.wavelength = s.wavelength;

  Panel p;
  p.fast = v(s.fast);
  p.slow = v(s.slow);
  p.origin = v(s.origin);
  p.pixel_size[0] = s.pixel_size[0];
  p.pixel_size[1] = s.pixel_size[1];
  p.image_size[0] = s.image_size[0];
  p.image_size[1] = s.image_size[1];
  p.parallax = true;
  p.mu = s.mu;
  p.thickness = s.thickness;
  e.detector.panels.push_back(p);

  std::vector<Vec3> axes;
  for (const std::array<double, 3> &axis : s.axes)
    axes.push_back(v(axis));
  e.goniometer = Goniometer::from_axes(axes, s.angles,
                                       static_cast<std::size_t>(s.scan_axis));

  e.scan.first_image = s.image_range[0];
  e.scan.last_image = s.image_range[1];
  e.scan.osc_start = s.osc_start;
  e.scan.osc_width = s.osc_width;
  e.crystal = Crystal::from_real_space(v(a), v(b), v(c));
  return e;
}

Experiment build(const cysteine::Sweep &s) {
  return build(s, cysteine::kRealSpaceA, cysteine::kRealSpaceB,
               cysteine::kRealSpaceC);
}
Experiment build(const threeaxis::Sweep &s) {
  return build(s, threeaxis::kRealSpaceA, threeaxis::kRealSpaceB,
               threeaxis::kRealSpaceC);
}

double difference(const Mat3 &a, const Mat3 &b) {
  double worst = 0.0;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      worst = std::fmax(worst, std::abs(a(i, j) - b(i, j)));
    }
  }
  return worst;
}

// Median disagreement with DIALS' own rlp column for one sweep, optionally
// with the goniometer replaced by a deliberately wrong one.
template <typename Sweep>
double disagreement(const Sweep &s, const Goniometer *replacement) {
  Experiment e = build(s);
  if (replacement)
    e.goniometer = *replacement;
  std::vector<double> residual;
  for (const auto &r : s.rows) {
    const Vec3 got =
        reciprocal_lattice_point(e, 0, r.px_fast, r.px_slow, r.px_z);
    residual.push_back((got - v(r.rlp)).norm());
  }
  std::sort(residual.begin(), residual.end());
  return residual[residual.size() / 2];
}

} // namespace

TEST(cysteine_has_the_four_sweeps_that_make_it_useful) {
  const std::vector<cysteine::Sweep> &all = cysteine::sweeps();
  check::equal(static_cast<long long>(all.size()), 4, "four sweeps");
  // Two sweeps with a non-zero setting angle below the scan axis, which is
  // what exercises the fixed rotation at all.
  int with_setting = 0, scan_axis_zero = 0;
  for (const cysteine::Sweep &s : all) {
    if (s.scan_axis == 1 && s.angles[0] != 0.0)
      ++with_setting;
    if (s.scan_axis == 0)
      ++scan_axis_zero;
  }
  check::equal(with_setting, 2, "sweeps with a non-identity fixed rotation");
  // And one that scans a different axis, so nothing may assume the scan axis
  // is the last one.
  check::equal(scan_axis_zero, 1, "sweeps scanning the first axis");
}

TEST(cysteine_all_sweeps_map_onto_the_dials_reciprocal_lattice) {
  for (const cysteine::Sweep &s : cysteine::sweeps()) {
    // Measured at 1.8e-5 to 3.8e-5, against a reciprocal cell edge of 0.082.
    check::is_true(disagreement(s, nullptr) < 1e-4,
                   "rlp must agree with DIALS on every sweep");
  }
}

TEST(cysteine_residual_against_the_lattice_is_the_models_not_ours) {
  // The residual against A h here is a hundredfold larger than for insulin:
  // 3e-3 to 9e-3 on a reciprocal cell edge of 0.082.
  //
  // That is NOT an error in this code, and it is not "what chemical data looks
  // like" either. It is the cost of a constraint that is provably false for
  // this experiment: one UB matrix and perfect goniometry shared across four
  // sweeps. The goniometer does not return to precisely the same place, so no
  // single matrix can fit all four. Breaking that constraint -- splitting the
  // sweeps and refining them separately about a common bulk matrix -- is the
  // first thing refinement does, and it improves these residuals greatly.
  //
  // Asserted so the distinction between a model constraint and a mapping error
  // is recorded rather than remembered.
  const Experiment e = build(cysteine::sweeps()[0]);
  std::vector<double> from_lattice, from_dials;
  for (const cysteine::Row &r : cysteine::sweeps()[0].rows) {
    const Vec3 ideal =
        e.crystal->A * Vec3{static_cast<double>(r.h), static_cast<double>(r.k),
                            static_cast<double>(r.l)};
    from_lattice.push_back((v(r.rlp) - ideal).norm());
    const Vec3 got =
        reciprocal_lattice_point(e, 0, r.px_fast, r.px_slow, r.px_z);
    from_dials.push_back((got - v(r.rlp)).norm());
  }
  std::sort(from_lattice.begin(), from_lattice.end());
  std::sort(from_dials.begin(), from_dials.end());
  const double model = from_lattice[from_lattice.size() / 2];
  const double ours = from_dials[from_dials.size() / 2];
  check::is_true(model > 1e-3, "the model residual is large on this data");
  check::is_true(ours * 20.0 < model,
                 "our disagreement with DIALS must be far below the model's");
}

// --------------------------------------------------------------------------
// The wrong compositions, and which sweep catches each
// --------------------------------------------------------------------------

TEST(cysteine_swapping_fixed_and_setting_is_caught) {
  // Breaks exactly the two sweeps with a non-zero phi setting, and leaves the
  // other two untouched -- which is why the insulin data could never have
  // caught this.
  int broken = 0;
  for (const cysteine::Sweep &s : cysteine::sweeps()) {
    Goniometer g = build(s).goniometer;
    std::swap(g.fixed, g.setting);
    if (disagreement(s, &g) > 0.01)
      ++broken;
  }
  check::equal(broken, 2, "two sweeps must detect the swap");
}

TEST(cysteine_dropping_the_fixed_rotation_is_caught) {
  int broken = 0;
  for (const cysteine::Sweep &s : cysteine::sweeps()) {
    Goniometer g = build(s).goniometer;
    g.fixed = Mat3::identity();
    if (disagreement(s, &g) > 0.01)
      ++broken;
  }
  check::equal(broken, 2, "two sweeps must detect the missing fixed rotation");
}

TEST(cysteine_folding_in_the_scan_axis_angle_is_caught) {
  // The rule that the scan axis's own entry in `angles` must be ignored,
  // because it is the scan start and the scan already supplies it. All three
  // omega sweeps carry -145 there, and their scans begin at -145.
  //
  // This is the variant that sweep 0 catches and the other wrong compositions
  // do not: its fixed rotation is otherwise the identity, so it is blind to
  // the swap and to the dropped rotation, and sensitive only to this.
  int broken = 0;
  bool sweep_zero_caught = false;
  const std::vector<cysteine::Sweep> &all = cysteine::sweeps();
  for (std::size_t i = 0; i < all.size(); ++i) {
    const cysteine::Sweep &s = all[i];
    Goniometer g = build(s).goniometer;
    g.fixed = g.fixed * rotation(v(s.axes[s.scan_axis]),
                                 Scan::radians(s.angles[s.scan_axis]));
    if (disagreement(s, &g) > 0.01) {
      ++broken;
      if (i == 0)
        sweep_zero_caught = true;
    }
  }
  check::equal(broken, 3, "the three omega sweeps must detect this");
  check::is_true(sweep_zero_caught, "sweep 0 is the one that pins this rule");
}

TEST(cysteine_scan_start_equals_the_scan_axis_setting_angle) {
  // Direct evidence for the rule above, independent of any mapping: the angle
  // recorded for the scan axis IS the scan's start angle. Composing it into
  // the fixed rotation would apply the same rotation twice.
  for (const cysteine::Sweep &s : cysteine::sweeps()) {
    check::close(s.angles[s.scan_axis], s.osc_start, 1e-9,
                 "scan-axis angle is the scan start");
  }
}

TEST(cysteine_parallax_uses_each_sweeps_own_sensor) {
  // A different detector from the insulin data: 320 micron silicon, mu 1.415,
  // 172 micron pixels, and one of the two panels at a 2theta offset. If the
  // sensor parameters were hard-coded anywhere this would move.
  for (const cysteine::Sweep &s : cysteine::sweeps()) {
    check::close(s.thickness, 0.32, 1e-9, "320 micron sensor");
    check::is_true(s.mu > 1.0 && s.mu < 2.0, "mu for this detector");
  }
  // The two detectors are genuinely different, so the sweeps are not all
  // secretly running through the same geometry.
  const Vec3 a = v(cysteine::sweeps()[0].origin);
  const Vec3 b = v(cysteine::sweeps()[3].origin);
  check::is_true((a - b).norm() > 10.0, "the two detectors differ");
}

// --------------------------------------------------------------------------
// Composition order, which no data to hand can settle
// --------------------------------------------------------------------------
//
// The axes run from the sample outwards to the laboratory, so an axis further
// out applies later and multiplies on the left. The l-cysteine goniometer has
// two axes, so at most one ever lies below the scan axis and the order is
// unobservable: every test above passes under either convention. These tests
// pin the implementation against the physical arrangement instead, which is
// the best that can be done until a three-circle instrument turns up.

namespace {

// Three mutually non-commuting axes, so that a wrong order is a different
// matrix rather than the same one.
const Vec3 kInner{1.0, 0.0, 0.0};
const Vec3 kMiddle{0.0, 1.0, 0.0};
const Vec3 kOuter = Vec3{0.3, 0.4, 0.866}.normalized();

} // namespace

TEST(axes_further_from_the_sample_apply_later) {
  // Sample on kInner, which sits on kMiddle, which sits on the scanned axis.
  const Goniometer g =
      Goniometer::from_axes({kInner, kMiddle, kOuter}, {30.0, 50.0, 0.0}, 2);
  const Mat3 inner = rotation(kInner, Scan::radians(30.0));
  const Mat3 middle = rotation(kMiddle, Scan::radians(50.0));

  check::close(difference(g.fixed, middle * inner), 0.0, 1e-12,
               "fixed must be outer-times-inner");
  // And the inside-out order must be a genuinely different matrix, or this
  // test would pass under either convention and prove nothing.
  check::is_true(difference(inner * middle, middle * inner) > 0.1,
                 "the two orders must actually differ");
}

TEST(the_inner_axis_is_itself_carried_by_the_outer_one) {
  // The physical statement behind the order: the direction of the axis nearest
  // the sample is not fixed in the laboratory -- the axis outside it moves it.
  // So the fixed rotation applied to the inner axis direction must equal the
  // outer rotation applied to it, the inner axis being invariant under its own
  // rotation.
  const Goniometer g =
      Goniometer::from_axes({kInner, kMiddle, kOuter}, {30.0, 50.0, 0.0}, 2);
  const Vec3 carried = g.fixed * kInner;
  const Vec3 expected = rotation(kMiddle, Scan::radians(50.0)) * kInner;
  check::close((carried - expected).norm(), 0.0, 1e-12,
               "the inner axis is carried by the outer one");
}

TEST(setting_rotation_composes_in_the_same_direction) {
  // Scanning the innermost axis, with two axes above it.
  const Goniometer g =
      Goniometer::from_axes({kInner, kMiddle, kOuter}, {0.0, 50.0, 20.0}, 0);
  const Mat3 middle = rotation(kMiddle, Scan::radians(50.0));
  const Mat3 outer = rotation(kOuter, Scan::radians(20.0));
  check::close(difference(g.setting, outer * middle), 0.0, 1e-12,
               "setting must be outer-times-middle");
  check::close(rotation_angle(g.fixed), 0.0, 1e-12,
               "nothing lies below the innermost axis");
}

// --------------------------------------------------------------------------
// Three axes -- and why that is still not enough
// --------------------------------------------------------------------------
//
// Insulin again, four 360 degree omega sweeps with chi stepped 0, 10, 20, 30.
// Three axes: phi on the sample, chi carrying it, omega carrying that and
// scanned. Two axes therefore lie below the scan axis, which is the
// arrangement the composition order was said to need.
//
// It is still not enough, and the reason is worth stating plainly: phi is zero
// in all four sweeps, so its rotation is the identity and the two orders
// produce the same matrix to the last digit. Checked directly, not assumed --
// see the test at the end.
//
// What is actually required is two axes below the scan axis at SIMULTANEOUSLY
// non-zero angles. "Three circles" is not the condition; a non-zero inner
// angle is.

TEST(threeaxis_has_three_axes_two_of_them_below_the_scan_axis) {
  const std::vector<threeaxis::Sweep> &all = threeaxis::sweeps();
  check::equal(static_cast<long long>(all.size()), 4, "four sweeps");
  for (const threeaxis::Sweep &s : all) {
    check::equal(static_cast<long long>(s.axes.size()), 3, "three axes");
    check::equal(s.scan_axis, 2, "the outermost axis is scanned");
  }
}

TEST(threeaxis_all_sweeps_map_onto_the_dials_reciprocal_lattice) {
  for (const threeaxis::Sweep &s : threeaxis::sweeps()) {
    // Measured at 1.1e-5 to 2.8e-5 on a reciprocal cell edge of 0.0148.
    check::is_true(disagreement(s, nullptr) < 1e-4,
                   "rlp must agree with DIALS on every sweep");
  }
}

TEST(threeaxis_middle_axis_setting_varies_and_reaches_the_fixed_rotation) {
  // chi is stepped 0, 10, 20, 30 across the sweeps on one crystal, so the
  // fixed rotation is different for each and all four still map correctly.
  // That is stronger than l-cysteine, which had two distinct settings.
  std::vector<double> angle;
  for (const threeaxis::Sweep &s : threeaxis::sweeps()) {
    angle.push_back(Scan::degrees(rotation_angle(build(s).goniometer.fixed)));
  }
  check::equal(static_cast<long long>(angle.size()), 4, "four settings");
  check::close(angle[0], 0.0, 1e-9, "chi = 0");
  check::close(angle[1], 10.0, 1e-9, "chi = 10");
  check::close(angle[2], 20.0, 1e-9, "chi = 20");
  check::close(angle[3], 30.0, 1e-9, "chi = 30");
}

TEST(threeaxis_dropping_the_fixed_rotation_is_caught_by_three_sweeps) {
  int broken = 0;
  for (const threeaxis::Sweep &s : threeaxis::sweeps()) {
    Goniometer g = build(s).goniometer;
    g.fixed = Mat3::identity();
    if (disagreement(s, &g) > 0.001)
      ++broken;
  }
  // All but the chi = 0 sweep, whose fixed rotation is the identity anyway.
  check::equal(broken, 3, "three sweeps must detect the missing rotation");
}

TEST(threeaxis_cannot_distinguish_the_composition_order) {
  // The honest result. Two axes below the scan axis is necessary but not
  // sufficient: with phi at zero the orders coincide exactly, so this dataset
  // says nothing about which is right, and neither does any other to hand.
  //
  // Asserted rather than noted, so that if a dataset with a non-zero phi is
  // ever substituted here this test fails and says so.
  for (const threeaxis::Sweep &s : threeaxis::sweeps()) {
    const Goniometer outward = build(s).goniometer;

    Mat3 inside_out = Mat3::identity();
    for (std::size_t i = 0; i < static_cast<std::size_t>(s.scan_axis); ++i) {
      inside_out =
          inside_out * rotation(v(s.axes[i]), Scan::radians(s.angles[i]));
    }
    check::close(difference(outward.fixed, inside_out), 0.0, 1e-15,
                 "the two orders coincide here, so this data cannot decide");
  }
  // The inner angle is the reason.
  for (const threeaxis::Sweep &s : threeaxis::sweeps()) {
    check::close(s.angles[0], 0.0, 1e-12, "phi is zero throughout");
  }
}

} // namespace mxi
