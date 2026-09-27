#include "postrefine.hh"

#include <algorithm>
#include <cmath>

namespace mxi {

std::vector<std::size_t> rows_for_postrefinement(const Table &integrated) {
  std::vector<std::size_t> rows;
  if (!integrated.has("flags") || !integrated.has("xyzres.px.value") ||
      !integrated.has("miller_index"))
    return rows;
  const Column &flags = integrated.at("flags");
  const Column &res = integrated.at("xyzres.px.value");
  const Column &miller = integrated.at("miller_index");
  for (std::size_t i = 0; i < integrated.nrows; ++i) {
    if ((flags.ints[i] & flag::kIntegratedSum) == 0)
      continue;
    if (!std::isfinite(res.real(i, 0)) || !std::isfinite(res.real(i, 1)) ||
        !std::isfinite(res.real(i, 2)))
      continue;
    if (miller.integer(i, 0) == 0 && miller.integer(i, 1) == 0 &&
        miller.integer(i, 2) == 0)
      continue;
    rows.push_back(i);
  }
  return rows;
}

std::size_t postrefinement_points(const Scan &scan) {
  const double degrees =
      std::abs(scan.osc_width) * static_cast<double>(scan.num_images());
  const std::size_t by_rotation =
      static_cast<std::size_t>(std::ceil(degrees / 36.0)) + 2;
  return std::max<std::size_t>(5, by_rotation);
}

PostrefineResult postrefine(ExperimentList &experiments,
                            const Table &integrated, std::size_t points) {
  PostrefineResult out;
  out.candidates = integrated.nrows;
  const std::vector<std::size_t> rows = rows_for_postrefinement(integrated);
  out.selected = rows.size();
  if (rows.empty() || experiments.size() == 0)
    return out;
  if (points == 0)
    points = postrefinement_points(experiments[0].scan);
  // mxi_refine's defaults, with analytic derivatives, as the effect was
  // measured: the crystal and the detector statically, then the crystal
  // scan-varying with the detector held where the static pass left it.
  RefineOptions options;
  options.analytic = true;
  out.refinement = refine_in_two_passes(
      experiments, select_rows(integrated, rows), options, points);
  return out;
}

} // namespace mxi
