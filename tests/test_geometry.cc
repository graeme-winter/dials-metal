// Geometry and prediction, tested against closed forms and against each other.
//
// The strongest test here is the round trip: predict a reflection from a known
// crystal, take the predicted detector position and scan angle, map it back
// into reciprocal space through the same geometry, and require that it lands
// on A h to within a tiny tolerance. That exercises the whole forward and
// reverse chain -- panel intersection, Ewald solution, goniometer rotation,
// scan angle -- and any sign error in one of them that is not exactly
// cancelled by a matching error in its inverse shows up immediately.
//
// What it cannot catch is a sign error that IS exactly cancelled, which is
// what a wrong convention shared between the forward and reverse maps looks
// like. Only a real .expt can catch that. See docs/conventions.md.

#include <array>
#include <cmath>

#include "../src/geometry.h"
#include "../src/linalg.h"
#include "../src/predict.h"
#include "check.h"

using namespace mxi;

namespace {

constexpr double kPi = 3.14159265358979323846;

// Insulin-like: cubic, a = 78 A, on a 4148 x 4362 detector at 200 mm.
Experiment insulin(double cell = 78.0) {
  Experiment e;
  e.beam.direction = {0.0, 0.0, 1.0};
  e.beam.wavelength = 0.9537;

  Panel p;
  p.name = "panel0";
  p.fast = {1.0, 0.0, 0.0};
  p.slow = {0.0, -1.0, 0.0};
  p.pixel_size[0] = p.pixel_size[1] = 0.075;
  p.image_size[0] = 4148;
  p.image_size[1] = 4362;
  // Beam centre near the middle of the panel.
  p.origin = {-0.075 * 2074.0, 0.075 * 2181.0, 200.0};
  e.detector.panels.push_back(p);

  e.goniometer.axis = {1.0, 0.0, 0.0};
  e.scan.first_image = 1;
  e.scan.last_image = 600;
  e.scan.osc_start = 0.0;
  e.scan.osc_width = 0.1;

  e.crystal = Crystal::from_real_space({cell, 0.0, 0.0}, {0.0, cell, 0.0},
                                       {0.0, 0.0, cell});
  return e;
}

// A crystal in a general orientation, so that nothing is accidentally
// diagonal. A cubic cell aligned with the axes will pass an inverse-transpose
// bug; this will not.
Crystal tilted_crystal(double cell = 78.0) {
  const Mat3 r = rotation({0.3, -0.5, 0.81}, 0.7);
  return Crystal::from_real_space(r * Vec3{cell, 0.0, 0.0},
                                  r * Vec3{0.0, cell, 0.0},
                                  r * Vec3{0.0, 0.0, cell});
}

}  // namespace

// --------------------------------------------------------------------------
// linear algebra
// --------------------------------------------------------------------------

TEST(mat3_inverse_round_trips) {
  const Mat3 m{1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 0.0};
  const Mat3 p = m * m.inverse();
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      check::close(p(i, j), i == j ? 1.0 : 0.0, 1e-12, "m * m^-1");
    }
  }
}

TEST(mat3_inverse_flags_a_singular_matrix) {
  const Mat3 singular{1, 2, 3, 2, 4, 6, 7, 8, 9};
  bool ok = true;
  singular.inverse(&ok);
  check::is_true(!ok, "singular matrix must be reported, not inverted");
}

TEST(rotation_is_orthogonal_and_has_the_right_angle) {
  const Vec3 axis{0.2, 0.9, -0.3};
  for (double angle : {0.0, 0.1, 1.0, 2.5, kPi - 1e-6}) {
    const Mat3 r = rotation(axis, angle);
    check::close(r.determinant(), 1.0, 1e-12, "det of a rotation");
    check::close(rotation_angle(r), angle, 1e-9, "recovered rotation angle");
    // Its own transpose must be its inverse.
    const Mat3 p = r * r.transpose();
    check::close(p(0, 0) + p(1, 1) + p(2, 2), 3.0, 1e-12, "R R^T = I");
  }
}

TEST(rotation_leaves_its_axis_alone) {
  const Vec3 axis = Vec3{1.0, 2.0, -3.0}.normalized();
  const Vec3 turned = rotation(axis, 1.234) * axis;
  check::close((turned - axis).norm(), 0.0, 1e-12, "axis is fixed");
}

TEST(rotation_is_right_handed) {
  // A quarter turn about z must take x to y, not to -y. This is the test that
  // fails if the sign convention in Rodrigues' formula is flipped, and every
  // downstream angle would then be negated.
  const Vec3 turned = rotation({0.0, 0.0, 1.0}, kPi / 2.0) * Vec3{1.0, 0.0, 0.0};
  check::close(turned.x, 0.0, 1e-12, "x component");
  check::close(turned.y, 1.0, 1e-12, "y component");
}

TEST(spd_solver_recovers_a_known_solution) {
  // A = L L^T for a known L, so the system is positive definite by
  // construction and the answer is exact rather than fitted.
  double a[9] = {4, 2, 1, 2, 5, 3, 1, 3, 6};
  const double x[3] = {1.0, -2.0, 3.0};
  double b[3] = {a[0] * x[0] + a[1] * x[1] + a[2] * x[2],
                 a[3] * x[0] + a[4] * x[1] + a[5] * x[2],
                 a[6] * x[0] + a[7] * x[1] + a[8] * x[2]};
  check::is_true(solve_spd(a, b, 3), "factorisation should succeed");
  for (int i = 0; i < 3; ++i) check::close(b[i], x[i], 1e-12, "solution");
}

TEST(spd_solver_refuses_an_indefinite_matrix) {
  double a[4] = {1, 2, 2, 1};  // eigenvalues 3 and -1
  double b[2] = {1, 1};
  check::is_true(!solve_spd(a, b, 2), "indefinite matrix must be refused");
}

// --------------------------------------------------------------------------
// crystal
// --------------------------------------------------------------------------

TEST(setting_matrix_is_the_inverse_not_the_inverse_transpose) {
  // A triclinic cell, where the two differ. On a cubic cell both are diagonal
  // and the wrong one passes, which is why this cell is deliberately nasty.
  const Crystal c = Crystal::from_real_space({10.0, 0.0, 0.0}, {3.0, 20.0, 0.0},
                                             {1.0, 2.0, 30.0});
  const Mat3 real = Mat3::from_rows({10.0, 0.0, 0.0}, {3.0, 20.0, 0.0},
                                    {1.0, 2.0, 30.0});
  const Mat3 p = real * c.A;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      check::close(p(i, j), i == j ? 1.0 : 0.0, 1e-12, "real . A = I");
    }
  }
}

TEST(cell_parameters_come_back_from_the_setting_matrix) {
  const Crystal c = Crystal::from_real_space({10.0, 0.0, 0.0}, {3.0, 20.0, 0.0},
                                             {1.0, 2.0, 30.0});
  const UnitCell u = c.cell();
  check::close(u.a, 10.0, 1e-10, "a");
  check::close(u.b, std::sqrt(9.0 + 400.0), 1e-10, "b");
  check::close(u.c, std::sqrt(1.0 + 4.0 + 900.0), 1e-10, "c");
  // Volume from the parameters must match the determinant of the real matrix.
  check::close(u.volume(), std::abs(c.A.inverse().determinant()), 1e-6,
               "volume");
}

TEST(cell_survives_a_general_rotation) {
  const UnitCell u = tilted_crystal().cell();
  check::close(u.a, 78.0, 1e-10, "a");
  check::close(u.alpha, 90.0, 1e-10, "alpha");
  check::close(u.gamma, 90.0, 1e-10, "gamma");
}

TEST(d_spacing_of_a_known_reflection) {
  const Crystal c = Crystal::from_real_space({78.0, 0.0, 0.0}, {0.0, 78.0, 0.0},
                                             {0.0, 0.0, 78.0});
  check::close(c.d_spacing(1, 0, 0), 78.0, 1e-10, "d(100)");
  check::close(c.d_spacing(2, 0, 0), 39.0, 1e-10, "d(200)");
  check::close(c.d_spacing(1, 1, 0), 78.0 / std::sqrt(2.0), 1e-10, "d(110)");
}

// --------------------------------------------------------------------------
// detector
// --------------------------------------------------------------------------

TEST(panel_lab_coordinate_and_intersection_are_inverses) {
  const Experiment e = insulin();
  const Panel &p = e.detector[0];
  for (double x : {0.0, 1.5, 2074.0, 4147.0}) {
    for (double y : {0.0, 900.25, 2181.0, 4361.0}) {
      const Vec3 lab = p.lab_coord(x, y);
      const auto back = p.intersect(lab.normalized() / e.beam.wavelength);
      check::is_true(back.has_value(), "ray from a panel point must hit it");
      check::close(back->first, x, 1e-9, "fast");
      check::close(back->second, y, 1e-9, "slow");
    }
  }
}

TEST(panel_accepts_a_ray_on_its_exact_corner) {
  // Solved in millimetres and divided back into pixels, the exact corner comes
  // out at -4e-13 rather than zero. Rejecting it loses reflections predicted on
  // an edge, and on a tiled detector loses rays that strike a seam.
  const Experiment e = insulin();
  const Panel &p = e.detector[0];
  const auto corner = p.intersect(p.lab_coord(0.0, 0.0).normalized());
  check::is_true(corner.has_value(), "the exact corner must be on the panel");
  check::close(corner->first, 0.0, 1e-9, "fast at the corner");
  check::close(corner->second, 0.0, 1e-9, "slow at the corner");
}

TEST(panel_rejects_rays_that_miss_or_go_backwards) {
  const Experiment e = insulin();
  const Panel &p = e.detector[0];
  // Straight back towards the source.
  check::is_true(!p.intersect(Vec3{0.0, 0.0, -1.0}).has_value(),
                 "backward ray must not hit");
  // Just past the far edge, which must still be rejected: the tolerance is a
  // nanopixel, not a licence to widen the panel.
  check::is_true(!p.intersect(p.lab_coord(4148.5, 100.0).normalized()).has_value(),
                 "a ray past the last pixel must not hit");
  // A ray that would hit the plane far outside the panel bounds.
  const Vec3 far = p.lab_coord(-500.0, -500.0);
  check::is_true(!p.intersect(far.normalized()).has_value(),
                 "ray outside the panel must not hit");
  // Parallel to the panel.
  check::is_true(!p.intersect(p.fast).has_value(), "parallel ray must not hit");
}

// --------------------------------------------------------------------------
// scan
// --------------------------------------------------------------------------

TEST(scan_angle_and_z_are_inverses) {
  Scan s;
  s.first_image = 1;
  s.last_image = 600;
  s.osc_start = 12.5;
  s.osc_width = 0.1;
  for (double z : {0.0, 0.5, 1.0, 300.0, 599.99}) {
    check::close(s.z_from_phi(s.phi_from_z(z)), z, 1e-9, "z round trip");
  }
  check::close(Scan::degrees(s.phi_from_z(0.0)), 12.5, 1e-12, "start angle");
  check::close(Scan::degrees(s.phi_from_z(600.0)), 72.5, 1e-12, "end angle");
}

TEST(scan_z_offset_is_visible_not_hidden) {
  // The known DIALS wrinkle: z is anchored to image number one, so a scan
  // starting at image 0 is out by one image unless z_offset says otherwise.
  // The point of this test is that the choice is a parameter and not a
  // silently baked-in constant.
  Scan s;
  s.first_image = 0;
  s.osc_start = 0.0;
  s.osc_width = 0.1;
  check::close(Scan::degrees(s.phi_from_z(10.0)), 1.0, 1e-12, "DIALS default");
  s.z_offset = -1.0;
  check::close(Scan::degrees(s.phi_from_z(10.0)), 1.1, 1e-12, "offset applied");
}

// --------------------------------------------------------------------------
// prediction
// --------------------------------------------------------------------------

TEST(a_point_beyond_the_ewald_diameter_never_diffracts) {
  const Experiment e = insulin();
  const double limit = 2.0 / e.beam.wavelength;
  const Intersections none = ewald_intersections(e, Vec3{limit * 1.1, 0.0, 0.0});
  check::is_true(!none.any, "beyond 2/lambda there is no solution");
}

TEST(ewald_solutions_actually_lie_on_the_sphere) {
  // The direct check of the mathematics: at each returned angle, the rotated
  // point plus s0 must have length exactly 1/lambda.
  Experiment e = insulin();
  e.crystal = tilted_crystal();
  const double radius = e.beam.s0().norm();
  int tested = 0;
  for (int h = -6; h <= 6; ++h) {
    for (int k = -6; k <= 6; ++k) {
      for (int l = -6; l <= 6; ++l) {
        if (!h && !k && !l) continue;
        const Vec3 r0 = e.crystal->A * Vec3{static_cast<double>(h),
                                            static_cast<double>(k),
                                            static_cast<double>(l)};
        const Intersections x = ewald_intersections(e, r0);
        if (!x.any) continue;
        for (int i = 0; i < 2; ++i) {
          const Vec3 s1 = e.beam.s0() + e.goniometer.rotation_at(x.phi[i]) * r0;
          check::close(s1.norm(), radius, 1e-12 * radius, "on the Ewald sphere");
          ++tested;
        }
      }
    }
  }
  check::is_true(tested > 100, "the test should have found solutions to check");
}

TEST(the_two_solutions_are_one_entering_and_one_exiting) {
  Experiment e = insulin();
  e.crystal = tilted_crystal();
  int checked = 0;
  for (int h = -5; h <= 5; ++h) {
    for (int k = -5; k <= 5; ++k) {
      for (int l = -5; l <= 5; ++l) {
        if (!h && !k && !l) continue;
        const Vec3 r0 = e.crystal->A * Vec3{static_cast<double>(h),
                                            static_cast<double>(k),
                                            static_cast<double>(l)};
        const Intersections x = ewald_intersections(e, r0);
        if (!x.any) continue;
        // A grazing intersection has both roots coincident and is allowed to
        // give the same flag twice; anything else must give one of each.
        if (std::abs(x.phi[0] - x.phi[1]) < 1e-6) continue;
        check::is_true(x.entering[0] != x.entering[1],
                       "one root entering, one exiting");
        ++checked;
      }
    }
  }
  check::is_true(checked > 50, "should have checked a decent number");
}

TEST(prediction_round_trips_through_reciprocal_space) {
  // The central test. Predict, then invert the whole chain and require the
  // recovered reciprocal lattice point to be A h again.
  Experiment e = insulin();
  e.crystal = tilted_crystal();
  PredictOptions options;
  options.d_min = 3.0;
  const std::vector<Prediction> predictions = predict(e, options);
  check::is_true(predictions.size() > 500,
                 "a 60 degree sweep of insulin to 3 A should give plenty");

  double worst = 0.0;
  for (const Prediction &p : predictions) {
    const Vec3 expected = e.crystal->A * Vec3{static_cast<double>(p.h),
                                              static_cast<double>(p.k),
                                              static_cast<double>(p.l)};
    const Vec3 recovered =
        reciprocal_lattice_point(e, p.panel, p.px_fast, p.px_slow, p.z);
    worst = std::fmax(worst, (recovered - expected).norm());
  }
  // In reciprocal Angstrom; the reciprocal cell edge is 1/78, so this is one
  // part in 10^8 of a lattice spacing.
  check::close(worst, 0.0, 1e-10, "worst round trip error");
}

TEST(predictions_land_inside_the_scan_and_the_detector) {
  Experiment e = insulin();
  e.crystal = tilted_crystal();
  PredictOptions options;
  options.d_min = 2.5;
  for (const Prediction &p : predict(e, options)) {
    check::is_true(p.z >= -1e-9 && p.z <= 600.0 + 1e-9, "z inside the scan");
    check::is_true(p.px_fast >= 0.0 && p.px_fast < 4148.0, "x on the panel");
    check::is_true(p.px_slow >= 0.0 && p.px_slow < 4362.0, "y on the panel");
  }
}

TEST(predictions_respect_the_resolution_limit) {
  Experiment e = insulin();
  e.crystal = tilted_crystal();
  PredictOptions options;
  options.d_min = 4.0;
  for (const Prediction &p : predict(e, options)) {
    const double d = e.crystal->d_spacing(p.h, p.k, p.l);
    check::is_true(d >= 4.0 - 1e-9, "no reflection beyond d_min");
  }
}

TEST(a_finer_resolution_limit_predicts_strictly_more) {
  Experiment e = insulin();
  e.crystal = tilted_crystal();
  PredictOptions coarse, fine;
  coarse.d_min = 4.0;
  fine.d_min = 3.0;
  check::is_true(predict(e, fine).size() > predict(e, coarse).size(),
                 "higher resolution must predict more reflections");
}

TEST(index_bounds_cover_the_requested_resolution) {
  const Crystal c = tilted_crystal();
  const std::array<int, 3> bounds = index_bounds(c, 3.0, 500);
  // 78 / 3 = 26, plus the rounding margin.
  for (int b : bounds) check::is_true(b >= 26 && b <= 30, "bound near 78/d_min");
  // Nothing inside the resolution shell may fall outside the box.
  for (int h = -40; h <= 40; ++h) {
    for (int k = -40; k <= 40; ++k) {
      for (int l = -40; l <= 40; ++l) {
        if (!h && !k && !l) continue;
        if (c.d_spacing(h, k, l) < 3.0) continue;
        check::is_true(std::abs(h) <= bounds[0] && std::abs(k) <= bounds[1] &&
                           std::abs(l) <= bounds[2],
                       "a reflection inside d_min fell outside the index box");
      }
    }
  }
}

TEST(predict_indices_agrees_with_the_full_enumeration) {
  Experiment e = insulin();
  e.crystal = tilted_crystal();
  PredictOptions options;
  options.d_min = 5.0;
  const std::vector<Prediction> all = predict(e, options);

  std::vector<std::array<int, 3>> asked;
  for (const Prediction &p : all) asked.push_back({p.h, p.k, p.l});
  const std::vector<Prediction> some = predict_indices(e, asked, options);
  check::is_true(some.size() >= all.size(),
                 "asking for the same indices must find at least as many");
}

TEST(a_crystal_free_experiment_predicts_nothing) {
  Experiment e = insulin();
  e.crystal.reset();
  check::equal(static_cast<long long>(predict(e).size()), 0,
               "no crystal, no predictions");
}
