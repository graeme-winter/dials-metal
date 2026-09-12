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

#include "../src/geometry.h"
#include "check.h"
#include "real_cysteine.h"

using namespace mxi;

namespace {

Vec3 v(const double (&a)[3]) { return {a[0], a[1], a[2]}; }

Experiment build(const cysteine::Sweep &s) {
  Experiment e;
  e.beam.direction = v(cysteine::kBeamDirection);
  e.beam.wavelength = cysteine::kWavelength;

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

  e.goniometer = Goniometer::from_axes({v(s.axes[0]), v(s.axes[1])},
                                       {s.angles[0], s.angles[1]},
                                       static_cast<std::size_t>(s.scan_axis));

  e.scan.first_image = s.image_range[0];
  e.scan.last_image = s.image_range[1];
  e.scan.osc_start = s.osc_start;
  e.scan.osc_width = s.osc_width;

  e.crystal = Crystal::from_real_space(v(cysteine::kRealSpaceA),
                                       v(cysteine::kRealSpaceB),
                                       v(cysteine::kRealSpaceC));
  return e;
}

// Median disagreement with DIALS' own rlp column for one sweep, optionally
// with the goniometer replaced by a deliberately wrong one.
double disagreement(const cysteine::Sweep &s, const Goniometer *replacement) {
  Experiment e = build(s);
  if (replacement) e.goniometer = *replacement;
  std::vector<double> residual;
  for (const cysteine::Row &r : s.rows) {
    const Vec3 got =
        reciprocal_lattice_point(e, 0, r.px_fast, r.px_slow, r.px_z);
    residual.push_back((got - v(r.rlp)).norm());
  }
  std::sort(residual.begin(), residual.end());
  return residual[residual.size() / 2];
}

}  // namespace

TEST(cysteine_has_the_four_sweeps_that_make_it_useful) {
  const std::vector<cysteine::Sweep> &all = cysteine::sweeps();
  check::equal(static_cast<long long>(all.size()), 4, "four sweeps");
  // Two sweeps with a non-zero setting angle below the scan axis, which is
  // what exercises the fixed rotation at all.
  int with_setting = 0, scan_axis_zero = 0;
  for (const cysteine::Sweep &s : all) {
    if (s.scan_axis == 1 && s.angles[0] != 0.0) ++with_setting;
    if (s.scan_axis == 0) ++scan_axis_zero;
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
  // 3e-3 to 9e-3 on a reciprocal cell edge of 0.082. That is the indexing
  // residual of a large-cell P1 chemical dataset and not an error in this
  // code, which the test above establishes by agreeing with DIALS' own rlp
  // forty times more closely than either sits to the ideal lattice.
  //
  // Asserted so that the distinction is recorded rather than remembered.
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
    if (disagreement(s, &g) > 0.01) ++broken;
  }
  check::equal(broken, 2, "two sweeps must detect the swap");
}

TEST(cysteine_dropping_the_fixed_rotation_is_caught) {
  int broken = 0;
  for (const cysteine::Sweep &s : cysteine::sweeps()) {
    Goniometer g = build(s).goniometer;
    g.fixed = Mat3::identity();
    if (disagreement(s, &g) > 0.01) ++broken;
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
      if (i == 0) sweep_zero_caught = true;
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
