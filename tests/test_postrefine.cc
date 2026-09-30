#include <cmath>
#include <limits>

#include "../src/postrefine.hh"
#include "../src/refl.hh"
#include "check.hh"

namespace mxi {

namespace {

//: Five integrated reflections, one for each reason post-refinement keeps a row
//: or leaves it out.
Table five_integrated() {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  Table t;
  t.nrows = 5;
  Column &miller = t.int_column("miller_index", "cctbx::miller::index<>", 3);
  Column &flags = t.int_column("flags", "std::size_t", 1);
  Column &res = t.real_column("xyzres.px.value", "vec3<double>", 3);
  Column &obs = t.real_column("xyzobs.px.value", "vec3<double>", 3);
  const std::int64_t hkl[5][3] = {
      {1, 2, 3}, {2, 0, 0}, {0, 0, 0}, {3, 1, 1}, {1, 1, 4}};
  const std::int64_t f[5] = {
      flag::kIntegratedSum,                        // kept
      flag::kIntegratedSum,                        // no centre of mass
      flag::kIntegratedSum,                        // not indexed
      flag::kForegroundIncludesBadPixels,          // not summed: crosses a gap
      flag::kIntegratedSum | flag::kIntegratedPrf, // kept
  };
  const double r[5] = {0.1, nan, 0.2, 0.3, -0.1};
  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t k = 0; k < 3; ++k) {
      miller.ints[i * 3 + k] = hkl[i][k];
      res.reals[i * 3 + k] = r[i];
      obs.reals[i * 3 + k] = static_cast<double>(10 * i + k);
    }
    flags.ints[i] = f[i];
  }
  return t;
}

} // namespace

TEST(post_refinement_keeps_the_summed_the_centred_and_the_indexed) {
  const std::vector<std::size_t> rows =
      rows_for_postrefinement(five_integrated());
  check::equal(static_cast<long long>(rows.size()), 2, "two of the five");
  check::equal(static_cast<long long>(rows[0]), 0, "the plain one");
  check::equal(static_cast<long long>(rows[1]), 4, "and the one also fitted");
}

TEST(selected_rows_carry_every_component_in_order) {
  const Table t = five_integrated();
  const Table s = select_rows(t, {4, 0});
  check::equal(static_cast<long long>(s.nrows), 2, "two rows");
  check::close(s.at("xyzobs.px.value").real(0, 0), 40.0, 0.0, "row 4 first");
  check::close(s.at("xyzobs.px.value").real(0, 2), 42.0, 0.0,
               "all three components");
  check::close(s.at("xyzobs.px.value").real(1, 1), 1.0, 0.0, "then row 0");
  check::equal(static_cast<long long>(s.at("miller_index").integer(0, 2)), 4,
               "integer columns too");
  bool refused = false;
  try {
    select_rows(t, {5});
  } catch (const std::exception &) {
    refused = true;
  }
  check::is_true(refused, "a row past the end is an error, not garbage");
}

TEST(post_refinement_points_follow_the_rotation) {
  Scan scan;
  scan.first_image = 1;
  scan.osc_width = 0.1;
  scan.last_image = 300; // 30 degrees: the five the effect was measured with
  check::equal(static_cast<long long>(postrefinement_points(scan)), 5,
               "30 degrees");
  // One per 10 degrees since 30 September 2026, as mxi_refine --scan-varying:
  // it was one per 36 and two more, 7 and 12 below.
  scan.last_image = 1800; // 180
  check::equal(static_cast<long long>(postrefinement_points(scan)), 18,
               "180 degrees");
  scan.last_image = 3600; // 360
  check::equal(static_cast<long long>(postrefinement_points(scan)), 36,
               "360 degrees");
}

TEST(a_scan_varying_crystal_has_a_control_point_per_ten_degrees) {
  // The default for mxi_refine --scan-varying and --postrefine alike: one per
  // 10 degrees, as the reference profiles' scan blocks, and at least five --
  // by the angle, not the number of images.
  const auto points = [](int images, double width) {
    Scan s;
    s.first_image = 1;
    s.last_image = images;
    s.osc_width = width;
    return static_cast<long long>(scan_varying_points(s));
  };
  check::equal(points(3600, 0.1), 36, "a full turn");
  check::equal(points(1800, 0.1), 18, "half a turn");
  check::equal(points(300, 0.1), 5, "30 degrees, the least");
  check::equal(points(7200, 0.05), 36, "a full turn in finer images");
}

} // namespace mxi
