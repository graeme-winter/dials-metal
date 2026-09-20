#include "refine.h"

#include <algorithm>
#include <chrono>
#include <thread>
#include <cmath>
#include <cstdio>

#include "derivatives.h"
#include "predict.h"

namespace mxi {

// Where refinement's time goes. File scope rather than a return value because
// refine() is on a public header and its signature is not worth changing for
// this. Placed by position, not by pattern: the last timer added by pattern
// landed inside a lambda that runs four hundred and fifty million times.
//: Threads for the Jacobian. Zero means hardware_concurrency, one means none.
//: A knob because a machine running several of these at once wants fewer, and
//: because a threading change that cannot be turned off cannot be measured
//: against the version without it.
std::size_t g_jacobian_threads = 0;
double g_jacobian_seconds = 0.0;
double g_normal_seconds = 0.0;

namespace {
double now_seconds() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}
}  // namespace


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
  //: Which row of the reflection table this came from, so that what the
  //: refinement used can be reported back.
  std::size_t row = 0;
  std::size_t experiment = 0;
  std::size_t panel = 0;
  int h = 0, k = 0, l = 0;
  double px_fast = 0.0, px_slow = 0.0, z = 0.0;
  double weight[3] = {1.0, 1.0, 1.0};
  bool active = true;
};

std::vector<TargetRow> gather(const ExperimentList &experiments,
                                const Table &reflections,
                                const RefineOptions &options,
                                std::size_t *ill_conditioned = nullptr) {
  std::vector<TargetRow> out;
  if (!reflections.has("miller_index")) return out;
  const Column &miller = reflections.at("miller_index");
  const Column &xyz = reflections.at("xyzobs.px.value");
  const bool has_id = reflections.has("id");
  const bool has_panel = reflections.has("panel");
  const bool has_var = reflections.has("xyzobs.px.variance");

  // Strength threshold, measured against this dataset's own median so it
  // travels between detectors and spot finders.
  double threshold = -1.0;
  const char *strength = reflections.has("n_signal") ? "n_signal"
                         : reflections.has("intensity.sum.value")
                             ? "intensity.sum.value"
                             : nullptr;
  if (options.strong_only && strength) {
    const Column &v = reflections.at(strength);
    std::vector<double> values;
    for (std::size_t i = 0; i < reflections.nrows; ++i) {
      values.push_back(v.integral ? static_cast<double>(v.integer(i)) : v.real(i));
    }
    std::nth_element(values.begin(), values.begin() + values.size() / 2,
                     values.end());
    threshold = values[values.size() / 2];
  }

  for (std::size_t i = 0; i < reflections.nrows; ++i) {
    if (threshold >= 0.0) {
      const Column &v = reflections.at(strength);
      const double value =
          v.integral ? static_cast<double>(v.integer(i)) : v.real(i);
      if (value < threshold) continue;
    }
    TargetRow o;
    o.row = i;
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

    // The rotation angle has to be determined by the data before it is worth
    // fitting. Measured on the geometry as it stands, not on the refined one:
    // the volume moves very little under refinement, and recomputing it each
    // macrocycle would let the set of reflections churn.
    if (options.min_volume > 0.0 && o.experiment < experiments.size()) {
      const PredictionState state =
          prediction_state(experiments[o.experiment], o.panel, o.h, o.k, o.l, o.z);
      if (!state.valid || std::abs(state.volume) < options.min_volume) {
        if (ill_conditioned != nullptr) ++*ill_conditioned;
        continue;
      }
    }
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
      const double panel_shift[6] = {shift[at],     shift[at + 1],
                                     shift[at + 2], shift[at + 3],
                                     shift[at + 4], shift[at + 5]};
      // perturb_panel lives beside its own derivative, so the two cannot
      // drift apart. A derivative that does not match how the model actually
      // moves is a silent failure: the refinement takes confident steps in a
      // direction that means nothing.
      for (Panel &p : e.detector.panels) p = perturb_panel(p, panel_shift);
    }
    if (layout.beam) {
      const std::size_t at = layout.beam_at(i);
      const double beam_shift[2] = {shift[at], shift[at + 1]};
      e.beam = perturb_beam(e.beam, beam_shift);
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

// The Jacobian of the residual, analytically. The residual is observed minus
// calculated, so every entry is the NEGATIVE of the derivative of the
// prediction -- which matches what the finite-difference branch produces,
// since it differences the residual and not the prediction.
//: How many threads to build the Jacobian with. One below a threshold, where
//: the work is smaller than the cost of starting threads.
std::size_t jacobian_thread_count(std::size_t observations) {
  if (g_jacobian_threads == 1) return 1;
  std::size_t wanted = g_jacobian_threads;
  if (wanted == 0) {
    wanted = std::thread::hardware_concurrency();
    if (wanted == 0) wanted = 1;
  }
  // Below a few thousand reflections the threads cost more than they save.
  const std::size_t by_work = observations / 2000;
  return std::max<std::size_t>(1, std::min(wanted, by_work));
}

void build_analytic_jacobian(const ExperimentList &experiments,
                             const std::vector<TargetRow> &observations,
                             const Layout &layout,
                             std::vector<std::vector<double>> *jacobian) {
  const std::size_t n = layout.size();
  jacobian->assign(n, std::vector<double>(observations.size() * 3, 0.0));

  // One reflection per unit of work. Every thread writes only the three
  // entries belonging to its own reflection, in every parameter's row, so no
  // two threads touch the same double and no locking is needed. The rows are
  // sized above, before any thread starts, because resizing a vector another
  // thread is reading from is not something a lock would fix.
  //
  // The reduction into the normal equations is NOT threaded: it is five per
  // cent of refinement, it sums into one small matrix, and a parallel sum in
  // a different order would change the last bits of the answer for nothing.
  const std::size_t threads = jacobian_thread_count(observations.size());
  const auto chunk = [&](std::size_t from, std::size_t to) {
  for (std::size_t i = from; i < to; ++i) {
    const TargetRow &o = observations[i];
    if (!o.active) continue;
    const Experiment &e = experiments[o.experiment];
    const PredictionState s =
        prediction_state(e, o.panel, o.h, o.k, o.l, o.z);
    if (!s.valid || s.volume == 0.0) continue;

    // The residual is in pixels and images; the derivatives are in millimetres
    // and radians. The parallax correction sits between the two, so the
    // conversion is a 2x2 Jacobian rather than a division by the pixel size.
    const Panel &p = e.detector[o.panel];
    double J[4];
    p.mm_to_px_jacobian(s.v.x / s.v.z, s.v.y / s.v.z, J);
    const double per_image =
        e.scan.osc_width != 0.0 ? 1.0 / Scan::radians(e.scan.osc_width) : 0.0;

    const auto place = [&](std::size_t parameter, const CentroidDerivative &d,
                           double weight) {
      const double dpx_fast = weight * (J[0] * d.dX + J[1] * d.dY);
      const double dpx_slow = weight * (J[2] * d.dX + J[3] * d.dY);
      const double dz = weight * d.dphi * per_image;
      (*jacobian)[parameter][i * 3 + 0] = -dpx_fast;
      (*jacobian)[parameter][i * 3 + 1] = -dpx_slow;
      (*jacobian)[parameter][i * 3 + 2] = -dz;
    };

    if (layout.crystal && e.crystal) {
      const auto d = crystal_derivatives(s, o.h, o.k, o.l);
      const std::size_t at = layout.crystal_at(o.experiment);
      if (layout.points < 2) {
        for (std::size_t k = 0; k < 9; ++k) place(at + k, d[k], 1.0);
      } else {
        // Only four control points are touched. This is the banding, and it
        // is the reason the B-spline was chosen over an interpolating spline.
        const SplineWeights w = spline_weights(e, o.z);
        for (std::size_t c = 0; c < w.count; ++c) {
          for (std::size_t k = 0; k < 9; ++k) {
            const std::size_t parameter = at + w.index[c] * 9 + k;
            const double dpx_fast = w.weight[c] * (J[0] * d[k].dX + J[1] * d[k].dY);
            const double dpx_slow = w.weight[c] * (J[2] * d[k].dX + J[3] * d[k].dY);
            const double dz = w.weight[c] * d[k].dphi * per_image;
            // Accumulated, not assigned: the padding repeats a control point
            // at the ends of the scan and both contributions are real.
            (*jacobian)[parameter][i * 3 + 0] -= dpx_fast;
            (*jacobian)[parameter][i * 3 + 1] -= dpx_slow;
            (*jacobian)[parameter][i * 3 + 2] -= dz;
          }
        }
      }
    }
    if (layout.detector) {
      const auto d = detector_derivatives(s, p);
      const std::size_t at = layout.detector_at(o.experiment);
      for (std::size_t k = 0; k < 6; ++k) place(at + k, d[k], 1.0);
    }
    if (layout.beam) {
      const auto d = beam_derivatives(s, e.beam);
      const std::size_t at = layout.beam_at(o.experiment);
      for (std::size_t k = 0; k < 2; ++k) place(at + k, d[k], 1.0);
    }
  }
  };

  if (threads <= 1) {
    chunk(0, observations.size());
    return;
  }
  std::vector<std::thread> pool;
  pool.reserve(threads - 1);
  const std::size_t each = (observations.size() + threads - 1) / threads;
  for (std::size_t t = 1; t < threads; ++t) {
    const std::size_t from = std::min(t * each, observations.size());
    const std::size_t to = std::min(from + each, observations.size());
    if (from < to) pool.emplace_back(chunk, from, to);
  }
  // This thread takes the first chunk rather than waiting for the others.
  chunk(0, std::min(each, observations.size()));
  for (std::thread &t : pool) t.join();
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
  std::size_t ill_conditioned = 0;
  std::vector<TargetRow> observations =
      gather(experiments, reflections, options, &ill_conditioned);
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
      const double t_jacobian = now_seconds();
      std::vector<std::vector<double>> jacobian(n);
      ExperimentList trial;
      if (options.analytic) {
        build_analytic_jacobian(experiments, observations, layout, &jacobian);
      } else {
        std::vector<double> shift(n, 0.0);
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
      }

      g_jacobian_seconds += now_seconds() - t_jacobian;
      const double t_normal = now_seconds();
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
      g_normal_seconds += now_seconds() - t_normal;
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

  result.n_ill_conditioned = ill_conditioned;
  result.rows_used.clear();
  result.rows_rejected.clear();
  for (const TargetRow &o : observations) {
    if (o.active) {
      result.rows_used.push_back(o.row);
    } else {
      result.rows_rejected.push_back(o.row);
    }
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

JacobianComparison compare_jacobians(const ExperimentList &experiments,
                                     const Table &reflections,
                                     const RefineOptions &options) {
  JacobianComparison out;
  std::size_t ill_conditioned = 0;
  std::vector<TargetRow> observations =
      gather(experiments, reflections, options, &ill_conditioned);
  if (observations.empty()) return out;

  Layout layout;
  layout.crystal = options.crystal;
  layout.detector = options.detector;
  layout.beam = options.beam;
  layout.shared_crystal = options.shared_crystal;
  layout.n_experiments = experiments.size();
  layout.crystal_blocks =
      options.crystal ? (options.shared_crystal ? 1 : experiments.size()) : 0;
  layout.points = std::max<std::size_t>(1, options.scan_points);
  const std::size_t n = layout.size();
  if (n == 0) return out;

  ExperimentList base = experiments;
  if (layout.points > 1) {
    for (Experiment &e : base) {
      if (e.crystal && e.crystal->A_points.size() != layout.points) {
        e.crystal->A_points.assign(layout.points, e.crystal->A);
      }
    }
  }

  std::vector<std::vector<double>> analytic;
  build_analytic_jacobian(base, observations, layout, &analytic);

  std::vector<double> residual;
  residuals_of(base, observations, &residual);
  const std::vector<double> step = step_sizes(base, layout);

  std::vector<double> relative;
  std::vector<double> shift(n, 0.0);
  ExperimentList trial;
  for (std::size_t p = 0; p < n; ++p) {
    std::fill(shift.begin(), shift.end(), 0.0);
    shift[p] = step[p];
    apply(base, layout, shift, &trial);
    std::vector<double> plus = residual;
    residuals_of(trial, observations, &plus, layout.owner(p));
    std::fill(shift.begin(), shift.end(), 0.0);
    shift[p] = -step[p];
    apply(base, layout, shift, &trial);
    std::vector<double> minus = residual;
    residuals_of(trial, observations, &minus, layout.owner(p));

    for (std::size_t i = 0; i < plus.size(); ++i) {
      const double numeric = (plus[i] - minus[i]) / (2.0 * step[p]);
      const double exact = analytic[p][i];
      const double size = std::fmax(std::abs(numeric), std::abs(exact));
      // Entries where both are essentially zero say nothing about agreement.
      if (size < 1e-3) continue;
      relative.push_back(std::abs(numeric - exact) / size);
    }
  }
  if (relative.empty()) return out;
  std::sort(relative.begin(), relative.end());
  out.compared = relative.size();
  out.median_relative = relative[relative.size() / 2];
  out.percentile_99 = relative[relative.size() * 99 / 100];
  out.percentile_999 = relative[relative.size() * 999 / 1000];
  out.worst_relative = relative.back();
  for (double v : relative) {
    if (v > 0.5) ++out.grossly_different;
  }
  return out;
}

void set_indexed_flags(Table &reflections) {
  if (!reflections.has("miller_index")) return;
  const Column &miller = reflections.at("miller_index");
  // Modified, not replaced: the strong bit dials.find_spots set has to survive.
  Column &flags = reflections.modify_int_column("flags", "std::size_t", 1);
  for (std::size_t i = 0; i < reflections.nrows; ++i) {
    const bool indexed = miller.integer(i, 0) != 0 || miller.integer(i, 1) != 0 ||
                         miller.integer(i, 2) != 0;
    // Cleared as well as set: a reflection that was indexed on an earlier pass
    // and is not any more must stop claiming to be.
    if (indexed) {
      flags.ints[i] |= flag::kIndexed;
    } else {
      flags.ints[i] &= ~flag::kIndexed;
    }
  }
}

void set_refinement_flags(const RefineResult &result, Table &reflections) {
  Column &flags = reflections.modify_int_column("flags", "std::size_t", 1);
  for (std::size_t i = 0; i < reflections.nrows; ++i) {
    flags.ints[i] &= ~flag::kUsedInRefinement;
  }
  for (std::size_t i = 0; i < reflections.nrows; ++i) {
    flags.ints[i] &= ~flag::kCentroidOutlier;
  }
  for (std::size_t row : result.rows_used) {
    if (row < reflections.nrows) flags.ints[row] |= flag::kUsedInRefinement;
  }
  for (std::size_t row : result.rows_rejected) {
    if (row < reflections.nrows) flags.ints[row] |= flag::kCentroidOutlier;
  }
}

void add_observed_columns(const ExperimentList &experiments, Table &reflections) {
  if (!reflections.has("xyzobs.px.value")) return;
  const Column &xyz = reflections.at("xyzobs.px.value");
  const bool has_variance = reflections.has("xyzobs.px.variance");
  const bool has_id = reflections.has("id");
  const bool has_panel = reflections.has("panel");

  Column &mm = reflections.real_column("xyzobs.mm.value", "vec3<double>", 3);
  Column &mm_variance =
      reflections.real_column("xyzobs.mm.variance", "vec3<double>", 3);
  for (std::size_t i = 0; i < reflections.nrows; ++i) {
    const std::size_t id =
        has_id ? static_cast<std::size_t>(
                     std::max<std::int64_t>(0, reflections.at("id").integer(i)))
               : 0;
    if (id >= experiments.size()) continue;
    const Experiment &e = experiments[id];
    const std::size_t panel =
        has_panel ? static_cast<std::size_t>(reflections.at("panel").integer(i)) : 0;
    if (panel >= e.detector.size()) continue;
    const Panel &p = e.detector[panel];

    const double px_fast = xyz.real(i, 0);
    const double px_slow = xyz.real(i, 1);
    const double z = xyz.real(i, 2);

    // Millimetres through the parallax-corrected mapping, not a division by
    // the pixel size: on a real Eiger the two differ in the fourth digit, and
    // DIALS uses the corrected one.
    const auto position = p.px_to_mm(px_fast, px_slow);
    const double phi = e.scan.phi_from_z(z);
    mm.reals[i * 3 + 0] = position.first;
    mm.reals[i * 3 + 1] = position.second;
    mm.reals[i * 3 + 2] = phi;

    if (has_variance) {
      // The same factors, squared. A variance in pixels becomes one in
      // millimetres and a variance in images becomes one in radians.
      const Column &v = reflections.at("xyzobs.px.variance");
      const double width = Scan::radians(e.scan.osc_width);
      mm_variance.reals[i * 3 + 0] = v.real(i, 0) * p.pixel_size[0] * p.pixel_size[0];
      mm_variance.reals[i * 3 + 1] = v.real(i, 1) * p.pixel_size[1] * p.pixel_size[1];
      mm_variance.reals[i * 3 + 2] = v.real(i, 2) * width * width;
    }

  }
}

void add_reciprocal_columns(const ExperimentList &experiments,
                            Table &reflections) {
  if (!reflections.has("xyzobs.px.value")) return;
  const Column &xyz = reflections.at("xyzobs.px.value");
  const bool has_id = reflections.has("id");
  const bool has_panel = reflections.has("panel");

  Column &s1_column = reflections.real_column("s1", "vec3<double>", 3);
  Column &rlp_column = reflections.real_column("rlp", "vec3<double>", 3);
  Column &entering = reflections.int_column("entering", "bool", 1);
  Column &imageset = reflections.int_column("imageset_id", "int", 1);

  for (std::size_t i = 0; i < reflections.nrows; ++i) {
    const std::size_t id =
        has_id ? static_cast<std::size_t>(
                     std::max<std::int64_t>(0, reflections.at("id").integer(i)))
               : 0;
    if (id >= experiments.size()) continue;
    const Experiment &e = experiments[id];
    const std::size_t panel =
        has_panel ? static_cast<std::size_t>(reflections.at("panel").integer(i)) : 0;
    if (panel >= e.detector.size()) continue;
    const Panel &p = e.detector[panel];

    const double px_fast = xyz.real(i, 0);
    const double px_slow = xyz.real(i, 1);
    const double z = xyz.real(i, 2);
    const auto position = p.px_to_mm(px_fast, px_slow);

    // The scattering vector the observation implies, of the same length as the
    // incident beam because the scattering is elastic.
    const Vec3 lab = p.lab_coord_mm(position.first, position.second);
    const double length = lab.norm();
    if (!(length > 0.0)) continue;
    const Vec3 s0 = e.beam.s0();
    const Vec3 s1 = lab * (s0.norm() / length);
    for (std::size_t k = 0; k < 3; ++k) s1_column.reals[i * 3 + k] = s1[k];

    // And the reciprocal lattice point it came from, rotated back into the
    // crystal's frame at the start of the scan.
    const Vec3 rlp = reciprocal_lattice_point(e, panel, px_fast, px_slow, z);
    for (std::size_t k = 0; k < 3; ++k) rlp_column.reals[i * 3 + k] = rlp[k];

    // Entering or exiting the Ewald sphere. Part of a reflection's identity,
    // not a detail: the same Miller index can be recorded both ways in one
    // scan, and a join on the index alone would merge them.
    //
    // The sign is GREATER than zero, the opposite of the test as usually
    // written, because s0 here points from the source towards the sample and
    // dxtbx's beam direction points the other way -- the same convention
    // difference that once cost 2.1 inverse Angstroms of residual. Determined
    // against DIALS' own flags rather than reasoned about: 13760 of 13766
    // agreed, the rest having a triple product within rounding of zero.
    entering.ints[i] = s1.dot(e.goniometer.lab_axis().cross(s0)) > 0.0 ? 1 : 0;
    imageset.ints[i] = static_cast<std::int64_t>(id);
  }
}

void update_predictions(const ExperimentList &experiments, Table &reflections) {
  if (!reflections.has("miller_index")) return;
  const Column &miller = reflections.at("miller_index");
  const Column &xyz = reflections.at("xyzobs.px.value");
  const bool has_id = reflections.has("id");
  const bool has_panel = reflections.has("panel");

  Column &cal = reflections.real_column("xyzcal.px", "vec3<double>", 3);
  Column &cal_mm = reflections.real_column("xyzcal.mm", "vec3<double>", 3);
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
    // The same prediction in the units DIALS refines in.
    const Panel &p = experiments[id].detector[panel];
    const auto mm = p.px_to_mm(cal.reals[i * 3 + 0], cal.reals[i * 3 + 1]);
    cal_mm.reals[i * 3 + 0] = mm.first;
    cal_mm.reals[i * 3 + 1] = mm.second;
    cal_mm.reals[i * 3 + 2] =
        experiments[id].scan.phi_from_z(cal.reals[i * 3 + 2]);
  }
}

}  // namespace mxi
