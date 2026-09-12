// mxi_residuals: decompose centroid residuals radially and tangentially.
//
//   mxi_residuals refined.expt refined.refl
//
// Why radial and tangential rather than x and y: a detector distance error, a
// wrong pixel size and a mis-modelled absorption depth are all radial, while a
// panel rotation about the beam is tangential. In x and y each of those is a
// pattern that depends on where you look; in radius and tangent they separate.
//
// The residual this was written to chase is a monotonic radial gradient of a
// few tenths of a pixel, present in DIALS' own output as well as here, and not
// explained. It is not the detector distance -- refinement is free to move that
// and does -- and it is not the absorption depth being the unconditional
// rather than the conditional mean, which was tested and makes no difference.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include <cstdlib>
#include <string>

#include "expt.h"
#include "linalg.h"
#include "refl.h"

using namespace mxi;

namespace {

struct Point {
  double radius = 0.0, radial = 0.0, tangential = 0.0, dx = 0.0, dy = 0.0;
  double x = 0.0, y = 0.0;  // relative to the beam centre
};

double median(std::vector<double> v) {
  if (v.empty()) return 0.0;
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  return v[v.size() / 2];
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::printf("usage: mxi_residuals EXPT REFL [--bins N] [--clip PX]\n");
    return 2;
  }
  int bins = 8;
  double clip = 1.5;
  for (int i = 3; i + 1 < argc; ++i) {
    if (std::string(argv[i]) == "--bins") bins = std::atoi(argv[i + 1]);
    if (std::string(argv[i]) == "--clip") clip = std::atof(argv[i + 1]);
  }

  try {
    const ExperimentList experiments = read_experiments(argv[1]);
    const Table t = read_reflections(argv[2]);
    if (!t.has("xyzcal.px")) {
      std::fprintf(stderr,
                   "mxi_residuals: no xyzcal.px; run mxi_refine first\n");
      return 1;
    }
    const Column &obs = t.at("xyzobs.px.value");
    const Column &cal = t.at("xyzcal.px");
    const Column &miller = t.at("miller_index");
    const bool has_id = t.has("id");

    std::vector<Point> points;
    for (std::size_t i = 0; i < t.nrows; ++i) {
      if (!miller.integer(i, 0) && !miller.integer(i, 1) && !miller.integer(i, 2)) {
        continue;
      }
      if (cal.real(i, 0) == 0.0 && cal.real(i, 1) == 0.0) continue;
      const std::size_t id =
          has_id ? static_cast<std::size_t>(std::max<std::int64_t>(0, t.at("id").integer(i)))
                 : 0;
      if (id >= experiments.size()) continue;
      const Experiment &e = experiments[id];
      const Panel &p = e.detector[0];

      // Beam centre: where the incident beam meets the panel.
      const auto centre = p.intersect(e.beam.s0());
      if (!centre) continue;

      const double dx = obs.real(i, 0) - cal.real(i, 0);
      const double dy = obs.real(i, 1) - cal.real(i, 1);
      if (std::abs(dx) > clip || std::abs(dy) > clip) continue;

      const double rx = cal.real(i, 0) - centre->first;
      const double ry = cal.real(i, 1) - centre->second;
      const double r = std::hypot(rx, ry);
      if (r < 1.0) continue;
      Point q;
      q.radius = r;
      q.radial = (dx * rx + dy * ry) / r;
      q.tangential = (-dx * ry + dy * rx) / r;
      q.dx = dx;
      q.dy = dy;
      q.x = rx;
      q.y = ry;
      points.push_back(q);
    }
    if (points.size() < 50) {
      std::fprintf(stderr, "mxi_residuals: only %zu usable reflections\n",
                   points.size());
      return 1;
    }

    std::sort(points.begin(), points.end(),
              [](const Point &a, const Point &b) { return a.radius < b.radius; });
    std::printf("%zu reflections, clipped at %.2f px\n\n", points.size(), clip);
    std::printf("  radius (px)        n     radial   tangential        dx        dy\n");

    const std::size_t per = points.size() / static_cast<std::size_t>(bins);
    for (int b = 0; b < bins; ++b) {
      const std::size_t from = static_cast<std::size_t>(b) * per;
      const std::size_t to =
          (b + 1 == bins) ? points.size() : from + per;
      std::vector<double> radial, tangential, dx, dy;
      for (std::size_t i = from; i < to; ++i) {
        radial.push_back(points[i].radial);
        tangential.push_back(points[i].tangential);
        dx.push_back(points[i].dx);
        dy.push_back(points[i].dy);
      }
      std::printf("  %5.0f - %5.0f  %7zu  %+9.4f    %+9.4f %+9.4f %+9.4f\n",
                  points[from].radius, points[to - 1].radius, to - from,
                  median(radial), median(tangential), median(dx), median(dy));
    }

    // The affine part of the residual field, which is everything the detector
    // model can express: two scales, a rotation about the beam, a shear, and a
    // translation. Fitting and removing it answers the only question that
    // matters about a residual pattern -- whether refinement could have
    // removed it and did not, or whether the model has no parameter for it.
    {
      double n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
      double tx = 0, txx = 0, txy = 0, ty = 0, tyx = 0, tyy = 0;
      for (const Point &q : points) {
        n += 1;
        sx += q.x;
        sy += q.y;
        sxx += q.x * q.x;
        sxy += q.x * q.y;
        syy += q.y * q.y;
        tx += q.dx;
        txx += q.dx * q.x;
        txy += q.dx * q.y;
        ty += q.dy;
        tyx += q.dy * q.x;
        tyy += q.dy * q.y;
      }
      double a[9] = {n, sx, sy, sx, sxx, sxy, sy, sxy, syy};
      double bx[3] = {tx, txx, txy};
      double by[3] = {ty, tyx, tyy};
      double a2[9];
      for (int i = 0; i < 9; ++i) a2[i] = a[i];
      if (solve_spd(a, bx, 3) && solve_spd(a2, by, 3)) {
        const double scale_fast = bx[1], scale_slow = by[2];
        const double turn = 0.5 * (bx[2] - by[1]);
        const double shear = 0.5 * (bx[2] + by[1]);
        std::printf("\n  affine part, all of which the detector model can express:\n");
        std::printf("    scale fast   %+8.4f %%\n", 100 * scale_fast);
        std::printf("    scale slow   %+8.4f %%\n", 100 * scale_slow);
        std::printf("    rotation     %+8.4f mrad about the beam\n", 1000 * turn);
        std::printf("    shear        %+8.4f mrad\n", 1000 * shear);
        std::printf("\n  radial residual AFTER removing the affine part:\n");
        for (int b = 0; b < bins; ++b) {
          const std::size_t from = static_cast<std::size_t>(b) * per;
          const std::size_t to = (b + 1 == bins) ? points.size() : from + per;
          std::vector<double> left;
          for (std::size_t i = from; i < to; ++i) {
            const Point &q = points[i];
            const double ex = q.dx - (bx[0] + bx[1] * q.x + bx[2] * q.y);
            const double ey = q.dy - (by[0] + by[1] * q.x + by[2] * q.y);
            left.push_back((ex * q.x + ey * q.y) / q.radius);
          }
          std::printf("    %5.0f - %5.0f  %+9.4f\n", points[from].radius,
                      points[to - 1].radius, median(left));
        }
      }
    }

    // A straight line through the radial residual. Its slope is what a
    // detector distance error would look like, so a non-zero slope after
    // refinement means the gradient is not a distance error -- refinement was
    // free to remove it and did not.
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    const double n = static_cast<double>(points.size());
    for (const Point &q : points) {
      sx += q.radius;
      sy += q.radial;
      sxx += q.radius * q.radius;
      sxy += q.radius * q.radial;
    }
    const double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    std::printf("\n  radial slope %+.4e px per px over %.0f px = %+.3f px\n",
                slope, points.back().radius - points.front().radius,
                slope * (points.back().radius - points.front().radius));
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_residuals: %s\n", e.what());
    return 1;
  }
}
