#include "refine.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "predict.h"

namespace mxi {

Residual centroid_residual(const Experiment &e, std::size_t panel, int h, int k,
                           int l, double px_fast, double px_slow, double z) {
  Residual out;
  if (!e.crystal) return out;
  if (panel >= e.detector.size()) return out;

  // The setting matrix at the observation's own scan position. DIALS predicts
  // scan-varying reflections iteratively, solving for the frame and then
  // re-evaluating A there; using the observed z instead is an approximation
  // whose error is dA/dz times the residual in z, which is a third of an image
  // out of hundreds. It is second order and it is not free -- it means this
  // residual is not exactly DIALS' residual for a scan-varying model.
  const Vec3 r0 = e.setting_at(z) * Vec3{static_cast<double>(h),
                                         static_cast<double>(k),
                                         static_cast<double>(l)};
  const Intersections cross = ewald_intersections(e, r0);
  if (!cross.any) return out;

  // Which of the two roots. Taking the nearer to the observed angle is what
  // makes this robust to a model that is still far out: the alternative, using
  // the `entering` flag, assumes the flag was computed under a model close
  // enough to be trusted, which at the start of refinement it is not.
  const double phi_obs = e.scan.phi_from_z(z);
  const double two_pi = 2.0 * 3.14159265358979323846;
  double best = 0.0;
  double best_gap = 1e30;
  for (int i = 0; i < 2; ++i) {
    double gap = std::fmod(cross.phi[i] - phi_obs, two_pi);
    if (gap > 0.5 * two_pi) gap -= two_pi;
    if (gap < -0.5 * two_pi) gap += two_pi;
    if (std::abs(gap) < best_gap) {
      best_gap = std::abs(gap);
      best = phi_obs + gap;
    }
  }

  const Vec3 s1 = e.beam.s0() + e.goniometer.rotation_at(best) * r0;
  const auto hit = e.detector[panel].intersect(s1);
  if (!hit) return out;

  out.valid = true;
  out.dx = px_fast - hit->first;
  out.dy = px_slow - hit->second;
  out.dz = z - e.scan.z_from_phi(best);
  return out;
}

namespace {

// One reflection, reduced to what the target function needs. Named TargetRow
// rather than Observation because geometry.h already has an Observation and
// the two would be ambiguous at namespace scope.
struct TargetRow {
  std::size_t experiment = 0;
  std::size_t panel = 0;
  int h = 0, k = 0, l = 0;
  double px_fast = 0.0, px_slow = 0.0, z = 0.0;
  double weight[3] = {1.0, 1.0, 1.0};
  bool active = true;
};

std::vector<TargetRow> gather(const ExperimentList &experiments,
                                const Table &reflections,
                                const RefineOptions &options) {
  std::vector<TargetRow> out;
  if (!reflections.has("miller_index")) return out;
  const Column &miller = reflections.at("miller_index");
  const Column &xyz = reflections.at("xyzobs.px.value");
  const bool has_id = reflections.has("id");
  const bool has_panel = reflections.has("panel");
  const bool has_var = reflections.has("xyzobs.px.variance");

  for (std::size_t i = 0; i < reflections.nrows; ++i) {
    TargetRow o;
    o.h = static_cast<int>(miller.integer(i, 0));
    o.k = static_cast<int>(miller.integer(i, 1));
    o.l = static_cast<int>(miller.integer(i, 2));
    if (!o.h && !o.k && !o.l) continue;
    o.experiment = has_id ? static_cast<std::size_t>(
                                std::max<std::int64_t>(0, reflections.at("id").integer(i)))
                          : 0;
    if (o.experiment >= experiments.size()) continue;
    o.panel = has_panel ? static_cast<std::size_t>(reflections.at("panel").integer(i)) : 0;
    o.px_fast = xyz.real(i, 0);
    o.px_slow = xyz.real(i, 1);
    o.z = xyz.real(i, 2);
    if (has_var && !options.unit_weights) {
      const Column &v = reflections.at("xyzobs.px.variance");
      for (std::size_t k = 0; k < 3; ++k) {
        const double variance = v.real(i, k);
        // A zero variance would carry infinite weight and silently decide the
        // refinement on its own, so it falls back to unit weight instead.
        o.weight[k] = variance > 0.0 ? 1.0 / variance : 1.0;
      }
    }
    o.weight[2] *= options.z_weight;
    out.push_back(o);
  }
  return out;
}

// The parameter vector is laid out as: the crystal blocks first (one, or one
// per experiment), then six detector parameters per experiment, then two beam
// parameters per experiment.
struct Layout {
  bool crystal = false, detector = false, beam = false;
  bool shared_crystal = true;
  std::size_t n_experiments = 0;
  std::size_t crystal_blocks = 0;
  std::size_t points = 1;  // control points per crystal

  std::size_t block() const { return 9 * points; }
  std::size_t size() const {
    return crystal_blocks * block() + (detector ? 6 * n_experiments : 0) +
           (beam ? 2 * n_experiments : 0);
  }
  std::size_t crystal_at(std::size_t experiment) const {
    return (shared_crystal ? 0 : experiment) * block();
  }
  std::size_t detector_at(std::size_t experiment) const {
    return crystal_blocks * block() + 6 * experiment;
  }
  std::size_t beam_at(std::size_t experiment) const {
    return crystal_blocks * block() + (detector ? 6 * n_experiments : 0) +
           2 * experiment;
  }

  // Which experiment a parameter belongs to, or -1 for one shared by all.
  // With a scan-varying model there can be two hundred parameters, and a
  // numerical Jacobian that re-evaluated every reflection for each of them
  // would be quadratic in the number of sweeps for no reason: a crystal
  // parameter of sweep two cannot move a reflection of sweep three.
  long owner(std::size_t p) const {
    if (p < crystal_blocks * block()) {
      return shared_crystal ? -1 : static_cast<long>(p / block());
    }
    std::size_t q = p - crystal_blocks * block();
    if (detector && q < 6 * n_experiments) return static_cast<long>(q / 6);
    if (detector) q -= 6 * n_experiments;
    return static_cast<long>(q / 2);
  }
};

// Apply a parameter shift to a copy of the models. Shifts are always applied
// to the state at the start of the step, never accumulated, so a rejected
// Levenberg step leaves nothing behind.
void apply(const ExperimentList &base, const Layout &layout,
           const std::vector<double> &shift, ExperimentList *out) {
  *out = base;
  for (std::size_t i = 0; i < out->size(); ++i) {
    Experiment &e = (*out)[i];
    if (layout.crystal && e.crystal) {
      const std::size_t at = layout.crystal_at(i);
      if (layout.points < 2) {
        for (std::size_t k = 0; k < 9; ++k) e.crystal->A.m[k] += shift[at + k];
      } else {
        for (std::size_t c = 0; c < layout.points; ++c) {
          for (std::size_t k = 0; k < 9; ++k) {
            e.crystal->A_points[c].m[k] += shift[at + c * 9 + k];
          }
        }
      }
    }
    if (layout.detector) {
      const std::size_t at = layout.detector_at(i);
      const Vec3 translation{shift[at], shift[at + 1], shift[at + 2]};
      const Vec3 turn{shift[at + 3], shift[at + 4], shift[at + 5]};
      const double angle = turn.norm();
      for (Panel &p : e.detector.panels) {
        if (angle > 0.0) {
          // Rotate about the panel centre, not the laboratory origin: a
          // rotation about a point two hundred millimetres away is mostly a
          // translation, and the two parameters would then be so correlated
          // that the normal matrix is near singular.
          const Vec3 centre = p.lab_coord_mm(
              0.5 * static_cast<double>(p.image_size[0]) * p.pixel_size[0],
              0.5 * static_cast<double>(p.image_size[1]) * p.pixel_size[1]);
          const Mat3 r = rotation(turn / angle, angle);
          p.fast = r * p.fast;
          p.slow = r * p.slow;
          p.origin = centre + r * (p.origin - centre);
        }
        p.origin += translation;
      }
    }
    if (layout.beam) {
      const std::size_t at = layout.beam_at(i);
      // Two tilts, about axes perpendicular to the beam, which is all the
      // freedom a direction has. A third would be a rotation about the beam
      // itself and would do nothing at all.
      const Vec3 d = e.beam.direction.normalized();
      Vec3 u = Vec3{0.0, 0.0, 1.0}.cross(d);
      if (u.norm() < 1e-6) u = Vec3{1.0, 0.0, 0.0}.cross(d);
      u = u.normalized();
      const Vec3 v = d.cross(u);
      e.beam.direction = (d + u * shift[at] + v * shift[at + 1]).normalized();
    }
  }
}

// `only` restricts evaluation to one experiment, leaving the rest of `out`
// alone. Used for the numerical Jacobian, where a parameter of one sweep
// cannot move another sweep's reflections.
double residuals_of(const ExperimentList &experiments,
                    const std::vector<TargetRow> &observations,
                    std::vector<double> *out, long only = -1) {
  if (only < 0) out->assign(observations.size() * 3, 0.0);
  double total = 0.0;
  for (std::size_t i = 0; i < observations.size(); ++i) {
    const TargetRow &o = observations[i];
    if (!o.active) continue;
    if (only >= 0 && o.experiment != static_cast<std::size_t>(only)) continue;
    const Residual r = centroid_residual(experiments[o.experiment], o.panel, o.h,
                                         o.k, o.l, o.px_fast, o.px_slow, o.z);
    if (!r.valid) continue;
    (*out)[i * 3 + 0] = r.dx;
    (*out)[i * 3 + 1] = r.dy;
    (*out)[i * 3 + 2] = r.dz;
    total += o.weight[0] * r.dx * r.dx + o.weight[1] * r.dy * r.dy +
             o.weight[2] * r.dz * r.dz;
  }
  return total;
}

std::vector<double> step_sizes(const ExperimentList &experiments,
                               const Layout &layout) {
  std::vector<double> step(layout.size(), 0.0);
  for (std::size_t i = 0; i < experiments.size(); ++i) {
    if (layout.crystal && experiments[i].crystal) {
      const std::size_t at = layout.crystal_at(i);
      const std::size_t n = layout.block();
      // Scaled to the matrix itself: the elements of A are around 1/60 for a
      // protein and 1/5 for a small molecule, and one absolute step cannot
      // suit both.
      double scale = 0.0;
      for (double v : experiments[i].crystal->A.m) scale = std::fmax(scale, std::abs(v));
      for (std::size_t k = 0; k < n; ++k) step[at + k] = 1e-6 * std::fmax(scale, 1e-6);
    }
    if (layout.detector) {
      const std::size_t at = layout.detector_at(i);
      for (std::size_t k = 0; k < 3; ++k) step[at + k] = 1e-4;      // mm
      for (std::size_t k = 3; k < 6; ++k) step[at + k] = 1e-6;      // radians
    }
    if (layout.beam) {
      const std::size_t at = layout.beam_at(i);
      step[at] = step[at + 1] = 1e-7;
    }
  }
  return step;
}

double robust_spread(std::vector<double> values) {
  if (values.empty()) return 0.0;
  std::nth_element(values.begin(), values.begin() + values.size() / 2, values.end());
  const double median = values[values.size() / 2];
  std::vector<double> deviation;
  deviation.reserve(values.size());
  for (double v : values) deviation.push_back(std::abs(v - median));
  std::nth_element(deviation.begin(), deviation.begin() + deviation.size() / 2,
                   deviation.end());
  // Scaled so it is comparable with a standard deviation for normal data.
  return 1.4826 * deviation[deviation.size() / 2];
}

}  // namespace

RefineResult refine(ExperimentList &experiments, const Table &reflections,
                    const RefineOptions &options) {
  RefineResult result;
  std::vector<TargetRow> observations = gather(experiments, reflections, options);
  if (observations.size() < 20) return result;

  Layout layout;
  layout.crystal = options.crystal;
  layout.detector = options.detector;
  layout.beam = options.beam;
  layout.shared_crystal = options.shared_crystal;
  layout.n_experiments = experiments.size();
  layout.crystal_blocks =
      options.crystal ? (options.shared_crystal ? 1 : experiments.size()) : 0;
  layout.points = std::max<std::size_t>(1, options.scan_points);

  // Seed the control points from the static matrix, so a scan-varying run
  // starts exactly where a static one would and can only improve on it.
  if (layout.points > 1) {
    for (Experiment &e : experiments) {
      if (!e.crystal) continue;
      if (e.crystal->A_points.size() != layout.points) {
        e.crystal->A_points.assign(layout.points, e.crystal->A);
      }
    }
  }

  const std::size_t n = layout.size();
  if (n == 0) return result;

  // A shared crystal must start from one crystal, or the first shift would be
  // applied to several different starting matrices and mean different things
  // for each.
  if (options.crystal && options.shared_crystal && experiments[0].crystal) {
    for (Experiment &e : experiments) e.crystal = experiments[0].crystal;
  }

  std::vector<double> residual;
  double lambda = 1e-3;

  for (int macro = 0; macro < std::max(1, options.macrocycles); ++macro) {
    double previous = residuals_of(experiments, observations, &residual);

    for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
      const std::vector<double> step = step_sizes(experiments, layout);

      // Numerical Jacobian: one forward difference per parameter. Central
      // differences would cost twice as much for an accuracy the Gauss-Newton
      // step does not need, since the step is recomputed every iteration.
      std::vector<std::vector<double>> jacobian(n);
      std::vector<double> shift(n, 0.0);
      ExperimentList trial;
      std::vector<double> moved(residual.size());
      for (std::size_t p = 0; p < n; ++p) {
        std::fill(shift.begin(), shift.end(), 0.0);
        shift[p] = step[p];
        apply(experiments, layout, shift, &trial);
        const long only = layout.owner(p);
        moved = residual;
        residuals_of(trial, observations, &moved, only);
        jacobian[p].resize(moved.size());
        for (std::size_t i = 0; i < moved.size(); ++i) {
          // d(residual)/d(parameter); residual is observed minus calculated,
          // so this is the negative of the derivative of the prediction.
          jacobian[p][i] = (moved[i] - residual[i]) / step[p];
        }
      }

      std::vector<double> normal(n * n, 0.0);
      std::vector<double> rhs(n, 0.0);
      for (std::size_t i = 0; i < observations.size(); ++i) {
        if (!observations[i].active) continue;
        for (std::size_t k = 0; k < 3; ++k) {
          const std::size_t row = i * 3 + k;
          const double w = observations[i].weight[k];
          for (std::size_t a = 0; a < n; ++a) {
            const double ja = jacobian[a][row];
            if (ja == 0.0) continue;
            rhs[a] -= w * ja * residual[row];
            for (std::size_t b = 0; b <= a; ++b) {
              normal[a * n + b] += w * ja * jacobian[b][row];
            }
          }
        }
      }
      for (std::size_t a = 0; a < n; ++a) {
        for (std::size_t b = a + 1; b < n; ++b) normal[a * n + b] = normal[b * n + a];
      }

      bool stepped = false;
      for (int attempt = 0; attempt < 8; ++attempt) {
        std::vector<double> damped = normal;
        for (std::size_t a = 0; a < n; ++a) damped[a * n + a] *= (1.0 + lambda);
        std::vector<double> solution = rhs;
        if (!solve_spd(damped.data(), solution.data(), n)) {
          lambda *= 10.0;
          continue;
        }
        // The sign: residual is observed minus calculated, and the Jacobian is
        // of the residual, so the normal equations already solve for the shift
        // that reduces it. Negating here was the first thing that made this
        // diverge smoothly.
        apply(experiments, layout, solution, &trial);
        std::vector<double> trial_residual;
        const double value = residuals_of(trial, observations, &trial_residual);
        if (value < previous) {
          experiments = trial;
          residual = trial_residual;
          lambda = std::fmax(lambda * 0.3, 1e-9);
          const double improvement = (previous - value) / std::fmax(previous, 1e-30);
          previous = value;
          stepped = true;
          ++result.iterations;
          if (improvement < options.convergence) {
            result.converged = true;
          }
          break;
        }
        lambda *= 10.0;
      }
      if (!stepped || result.converged) break;
    }
    result.converged = false;

    // Outlier rejection between macrocycles, never inside one: rejecting while
    // the model is still moving throws away reflections for being far from a
    // prediction that was wrong.
    if (options.outlier_sigma > 0.0 && macro + 1 < options.macrocycles) {
      std::vector<double> dx, dy, dz;
      for (std::size_t i = 0; i < observations.size(); ++i) {
        if (!observations[i].active) continue;
        dx.push_back(residual[i * 3 + 0]);
        dy.push_back(residual[i * 3 + 1]);
        dz.push_back(residual[i * 3 + 2]);
      }
      const double sx = robust_spread(dx), sy = robust_spread(dy), sz = robust_spread(dz);
      std::size_t rejected = 0;
      for (std::size_t i = 0; i < observations.size(); ++i) {
        if (!observations[i].active) continue;
        const bool bad =
            (sx > 0 && std::abs(residual[i * 3 + 0]) > options.outlier_sigma * sx) ||
            (sy > 0 && std::abs(residual[i * 3 + 1]) > options.outlier_sigma * sy) ||
            (sz > 0 && std::abs(residual[i * 3 + 2]) > options.outlier_sigma * sz);
        if (bad) {
          observations[i].active = false;
          ++rejected;
        }
      }
      result.n_rejected += rejected;
      if (options.verbose) {
        std::printf("  macrocycle %d: rejected %zu outliers\n", macro + 1, rejected);
      }
    }

    if (options.verbose) {
      double sx = 0, sy = 0, sz = 0;
      std::size_t count = 0;
      for (std::size_t i = 0; i < observations.size(); ++i) {
        if (!observations[i].active) continue;
        sx += residual[i * 3 + 0] * residual[i * 3 + 0];
        sy += residual[i * 3 + 1] * residual[i * 3 + 1];
        sz += residual[i * 3 + 2] * residual[i * 3 + 2];
        ++count;
      }
      if (count) {
        std::printf("  macrocycle %d: %zu refl, rmsd %.4f %.4f %.4f px,px,images\n",
                    macro + 1, count, std::sqrt(sx / static_cast<double>(count)),
                    std::sqrt(sy / static_cast<double>(count)),
                    std::sqrt(sz / static_cast<double>(count)));
      }
    }
  }

  // Keep the static matrix consistent with the scan-varying one, taking the
  // middle of the scan. Otherwise `A` is whatever it was before the control
  // points moved, and every cell reported from it is stale -- which is how a
  // refinement comes to print a cell it does not believe.
  for (Experiment &e : experiments) {
    if (e.crystal && e.crystal->scan_varying()) e.crystal->A = e.crystal->A_at(0.5);
  }

  double sx = 0, sy = 0, sz = 0;
  std::size_t count = 0;
  for (std::size_t i = 0; i < observations.size(); ++i) {
    if (!observations[i].active) continue;
    sx += residual[i * 3 + 0] * residual[i * 3 + 0];
    sy += residual[i * 3 + 1] * residual[i * 3 + 1];
    sz += residual[i * 3 + 2] * residual[i * 3 + 2];
    ++count;
  }
  result.n_used = count;
  if (count) {
    result.rmsd_x = std::sqrt(sx / static_cast<double>(count));
    result.rmsd_y = std::sqrt(sy / static_cast<double>(count));
    result.rmsd_z = std::sqrt(sz / static_cast<double>(count));
  }
  return result;
}

void update_predictions(const ExperimentList &experiments, Table &reflections) {
  if (!reflections.has("miller_index")) return;
  const Column &miller = reflections.at("miller_index");
  const Column &xyz = reflections.at("xyzobs.px.value");
  const bool has_id = reflections.has("id");
  const bool has_panel = reflections.has("panel");

  Column &cal = reflections.real_column("xyzcal.px", "vec3<double>", 3);
  for (std::size_t i = 0; i < reflections.nrows; ++i) {
    const int h = static_cast<int>(miller.integer(i, 0));
    const int k = static_cast<int>(miller.integer(i, 1));
    const int l = static_cast<int>(miller.integer(i, 2));
    if (!h && !k && !l) continue;
    const std::size_t id =
        has_id ? static_cast<std::size_t>(
                     std::max<std::int64_t>(0, reflections.at("id").integer(i)))
               : 0;
    if (id >= experiments.size()) continue;
    const std::size_t panel =
        has_panel ? static_cast<std::size_t>(reflections.at("panel").integer(i)) : 0;
    const Residual r =
        centroid_residual(experiments[id], panel, h, k, l, xyz.real(i, 0),
                          xyz.real(i, 1), xyz.real(i, 2));
    if (!r.valid) continue;
    cal.reals[i * 3 + 0] = xyz.real(i, 0) - r.dx;
    cal.reals[i * 3 + 1] = xyz.real(i, 1) - r.dy;
    cal.reals[i * 3 + 2] = xyz.real(i, 2) - r.dz;
  }
}

}  // namespace mxi
