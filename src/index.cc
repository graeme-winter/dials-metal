#include "index.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <numeric>

#include "fft.hh"
#include "refine.hh"

namespace mxi {

namespace {

// Wall clock in seconds. Wall rather than CPU: the question is how long
// someone waits, and a thread count that changes the answer is part of it.
//: The last transform and peak search, in seconds. A file-scope pair rather
//: than a return value because find_candidate_vectors is on a public header
//: and its signature is not worth changing to answer one question about where
//: four seconds go.
std::size_t g_triples_scored = 0;
std::size_t g_triples_skipped = 0;
double g_last_fft_seconds = 0.0;
double g_last_peak_seconds = 0.0;

double now_seconds() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

} // namespace

std::vector<Vec3> reciprocal_lattice_points(const ExperimentList &experiments,
                                            const Table &reflections) {
  const Column &xyz = reflections.at("xyzobs.px.value");
  const bool has_id = reflections.has("id");
  const bool has_panel = reflections.has("panel");

  std::vector<Vec3> out;
  out.reserve(reflections.nrows);
  for (std::size_t i = 0; i < reflections.nrows; ++i) {
    const std::size_t id =
        has_id ? static_cast<std::size_t>(
                     std::max<std::int64_t>(0, reflections.at("id").integer(i)))
               : 0;
    if (id >= experiments.size()) {
      throw ReflError("reflection " + std::to_string(i) +
                      " has experiment id " + std::to_string(id) +
                      " but there are only " +
                      std::to_string(experiments.size()) + " experiments");
    }
    const std::size_t panel =
        has_panel ? static_cast<std::size_t>(reflections.at("panel").integer(i))
                  : 0;
    out.push_back(reciprocal_lattice_point(experiments[id], panel,
                                           xyz.real(i, 0), xyz.real(i, 1),
                                           xyz.real(i, 2)));
  }
  return out;
}

std::vector<int> observation_groups(const ExperimentList &experiments,
                                    const Table &reflections,
                                    double block_degrees) {
  std::vector<int> groups(reflections.nrows, 0);
  if (!reflections.has("xyzobs.px.value"))
    return groups;
  const Column &obs = reflections.at("xyzobs.px.value");
  const bool has_id = reflections.has("id");
  for (std::size_t i = 0; i < reflections.nrows; ++i) {
    std::size_t which = 0;
    if (has_id) {
      const std::int64_t id = reflections.at("id").integer(i);
      if (id >= 0 && static_cast<std::size_t>(id) < experiments.size()) {
        which = static_cast<std::size_t>(id);
      }
    }
    const Scan &scan = experiments[which].scan;
    const double phi = Scan::degrees(scan.phi_from_z(obs.real(i, 2)));
    // Which angular block, counted from zero. A sweep shorter than one block
    // puts everything in one group and nothing changes.
    const double width = block_degrees > 0.0 ? block_degrees : 30.0;
    const int block = static_cast<int>(std::floor(phi / width));
    groups[i] = static_cast<int>(which) * 100000 + block;
  }
  return groups;
}

double estimate_max_cell(const std::vector<Vec3> &points,
                         const std::vector<int> &groups) {
  if (points.size() < 4)
    return 0.0;
  if (groups.size() != points.size())
    return estimate_max_cell(points);
  const std::size_t sample = std::min<std::size_t>(500, points.size());
  const std::size_t step = std::max<std::size_t>(1, points.size() / sample);
  constexpr double kCoincident = 1e-4;

  std::vector<double> nearest;
  for (std::size_t i = 0; i < points.size(); i += step) {
    double best = 1e30;
    for (std::size_t j = 0; j < points.size(); ++j) {
      if (i == j || groups[i] != groups[j])
        continue;
      const double d = (points[i] - points[j]).norm_squared();
      if (d > kCoincident * kCoincident && d < best)
        best = d;
    }
    if (best < 1e29)
      nearest.push_back(std::sqrt(best));
  }
  if (nearest.empty())
    return 0.0;
  std::nth_element(nearest.begin(), nearest.begin() + nearest.size() / 2,
                   nearest.end());
  const double median = nearest[nearest.size() / 2];
  if (!(median > 0.0))
    return 0.0;
  return 1.5 / median;
}

double estimate_max_cell(const std::vector<Vec3> &points) {
  if (points.size() < 4)
    return 0.0;
  // The nearest-neighbour distance between reciprocal lattice points is, for a
  // well-sampled lattice, the shortest reciprocal cell edge -- so its
  // reciprocal is the longest real-space edge. The median is used rather than
  // the minimum because a handful of near-coincident spots (a split peak, two
  // lattices) would otherwise set the cell for the whole dataset.
  //
  // Sampled rather than exhaustive: this is O(sample * n), and the median of a
  // few hundred nearest-neighbour distances is as good as the median of tens
  // of thousands for a number that only has to be roughly right.
  const std::size_t sample = std::min<std::size_t>(500, points.size());
  const std::size_t step = std::max<std::size_t>(1, points.size() / sample);

  // Coincident points are skipped, and this is not a guard against bad data --
  // it is a property of the input. Pooling several sweeps of one crystal puts
  // the same reflection, measured at different goniometer settings, at the
  // same place in the crystal frame. That is the whole point of measuring more
  // than one sweep. Counting those pairs makes the nearest-neighbour spacing
  // vanish and the estimated cell diverge: on two sweeps of insulin it came
  // out at 9e15 Angstrom and no candidate vectors were found at all.
  //
  // The floor corresponds to a ten-thousand Angstrom cell, which is an order
  // of magnitude beyond anything this will ever see, so nothing real is lost.
  constexpr double kCoincident = 1e-4;

  std::vector<double> nearest;
  for (std::size_t i = 0; i < points.size(); i += step) {
    double best = 1e30;
    for (std::size_t j = 0; j < points.size(); ++j) {
      if (i == j)
        continue;
      const double d = (points[i] - points[j]).norm_squared();
      if (d > kCoincident * kCoincident && d < best)
        best = d;
    }
    if (best < 1e29)
      nearest.push_back(std::sqrt(best));
  }
  if (nearest.empty())
    return 0.0;
  std::nth_element(nearest.begin(), nearest.begin() + nearest.size() / 2,
                   nearest.end());
  const double median = nearest[nearest.size() / 2];
  if (!(median > 0.0))
    return 0.0;
  // A margin, because the estimate is a lower bound on the cell: a lattice
  // seen at less than full sampling gives neighbours further apart than the
  // true cell spacing, never closer.
  return 1.5 / median;
}

namespace {

struct Peak {
  Vec3 position;
  double height = 0.0;
};

} // namespace

std::vector<Vec3> find_candidate_vectors(const std::vector<Vec3> &points,
                                         const IndexOptions &options,
                                         double d_min, double max_cell,
                                         std::size_t grid) {
  const std::size_t n = grid;
  const double spacing = 1.0 / (2.0 * max_cell); // reciprocal grid spacing
  const double q_max = 1.0 / d_min;

  std::vector<std::complex<double>> f(n * n * n,
                                      std::complex<double>(0.0, 0.0));
  std::size_t used = 0;
  for (const Vec3 &r : points) {
    if (r.norm() > q_max)
      continue;
    // Nearest grid point, wrapped. Wrapping is not a nicety: the transform is
    // periodic, so a point placed at index -3 belongs at n-3 and anywhere else
    // is a different reflection.
    long index[3];
    bool ok = true;
    for (int k = 0; k < 3; ++k) {
      const long g = std::lround(r[static_cast<std::size_t>(k)] / spacing);
      if (std::abs(g) > static_cast<long>(n) / 2) {
        ok = false;
        break;
      }
      index[k] = (g + static_cast<long>(n)) % static_cast<long>(n);
    }
    if (!ok)
      continue;
    f[(static_cast<std::size_t>(index[0]) * n +
       static_cast<std::size_t>(index[1])) *
          n +
      static_cast<std::size_t>(index[2])] += 1.0;
    ++used;
  }
  if (used < 10)
    return {};

  const double t_fft = now_seconds();
  fft3d(f, n, +1);
  g_last_fft_seconds = now_seconds() - t_fft;

  // Peak search on the modulus. A grid point is a peak if it is at least as
  // large as all twenty-six of its neighbours; ties are broken by taking the
  // first, which only matters for exactly flat regions that are not peaks.
  const double real_spacing = 1.0 / (static_cast<double>(n) * spacing);
  const auto position = [&](long i, long j, long k) {
    const auto wrap = [&](long v) {
      return v > static_cast<long>(n) / 2 ? v - static_cast<long>(n) : v;
    };
    return Vec3{static_cast<double>(wrap(i)) * real_spacing,
                static_cast<double>(wrap(j)) * real_spacing,
                static_cast<double>(wrap(k)) * real_spacing};
  };

  std::vector<double> modulus(f.size());
  const double t_peaks = now_seconds();
  for (std::size_t i = 0; i < f.size(); ++i)
    modulus[i] = std::abs(f[i]);

  const long ln = static_cast<long>(n);
  const auto at = [&](long i, long j, long k) -> double {
    const auto w = [&](long v) { return ((v % ln) + ln) % ln; };
    return modulus[(static_cast<std::size_t>(w(i)) * n +
                    static_cast<std::size_t>(w(j))) *
                       n +
                   static_cast<std::size_t>(w(k))];
  };

  std::vector<Peak> peaks;
  for (long i = 0; i < ln; ++i) {
    for (long j = 0; j < ln; ++j) {
      for (long k = 0; k < ln; ++k) {
        const double v = at(i, j, k);
        if (v <= 0.0)
          continue;
        const Vec3 x = position(i, j, k);
        const double length = x.norm();
        if (length < options.min_cell || length > max_cell)
          continue;
        bool best = true;
        for (long di = -1; di <= 1 && best; ++di) {
          for (long dj = -1; dj <= 1 && best; ++dj) {
            for (long dk = -1; dk <= 1; ++dk) {
              if (!di && !dj && !dk)
                continue;
              if (at(i + di, j + dj, k + dk) > v) {
                best = false;
                break;
              }
            }
          }
        }
        if (!best)
          continue;
        // Centroid over the immediate neighbourhood. The grid step is around
        // 0.8 Angstrom, which is an eighth of the smallest cell edge worth
        // finding, so an uncentroided peak is not a usable basis vector.
        Vec3 centre{0.0, 0.0, 0.0};
        double weight = 0.0;
        for (long di = -1; di <= 1; ++di) {
          for (long dj = -1; dj <= 1; ++dj) {
            for (long dk = -1; dk <= 1; ++dk) {
              const double w = at(i + di, j + dj, k + dk);
              centre += position(i + di, j + dj, k + dk) * w;
              weight += w;
            }
          }
        }
        peaks.push_back({weight > 0.0 ? centre / weight : x, v});
      }
    }
  }

  std::sort(peaks.begin(), peaks.end(),
            [](const Peak &a, const Peak &b) { return a.height > b.height; });

  // A peak and its negation describe the same lattice vector, and so do two
  // peaks a fraction of a grid step apart. Keep the stronger of each.
  std::vector<Vec3> out;
  for (const Peak &p : peaks) {
    bool duplicate = false;
    for (const Vec3 &kept : out) {
      const double tolerance = 0.5 * real_spacing + 0.02 * kept.norm();
      if ((p.position - kept).norm() < tolerance ||
          (p.position + kept).norm() < tolerance) {
        duplicate = true;
        break;
      }
    }
    if (duplicate)
      continue;
    out.push_back(p.position);
    if (out.size() >= options.n_candidates)
      break;
  }
  g_last_peak_seconds = now_seconds() - t_peaks;
  return out;
}

// Count reflections whose fractional indices are all within `tolerance` of an
// integer, and accumulate the squared miss for the RMS.
static std::size_t score_basis(const Mat3 &real_rows,
                               const std::vector<Vec3> &points,
                               double tolerance, double *sum_squared) {
  std::size_t count = 0;
  double total = 0.0;
  for (const Vec3 &r : points) {
    // h = M r, M having the real-space basis vectors as its rows, because
    // h_i is the dot product of the i-th real vector with the scattering
    // vector. No matrix inversion is needed anywhere in the scoring loop.
    const Vec3 h = real_rows * r;
    double worst = 0.0;
    double miss = 0.0;
    for (std::size_t k = 0; k < 3; ++k) {
      // rint, not round. They differ only on an exact tie at .5, which is
      // outside every tolerance used here -- a tolerance of 0.5 or more would
      // accept every reflection -- so the result cannot change. The cost does:
      // round has round-half-away-from-zero semantics that no instruction
      // implements, so it compiles to a sequence, while rint is the hardware
      // rounding instruction where there is one. Measured at 3.51 against 1.22
      // nanoseconds a value, and there are three of them per reflection per
      // triple scored.
      const double d = h[k] - std::rint(h[k]);
      worst = std::fmax(worst, std::abs(d));
      miss += d * d;
    }
    if (worst < tolerance) {
      ++count;
      total += miss;
    }
  }
  if (sum_squared)
    *sum_squared = total;
  return count;
}

Mat3 reduce_basis(const Mat3 &real_space_rows) {
  // Repeatedly shorten one basis vector by a multiple of another. This is
  // Buerger reduction by the obvious greedy route rather than Niggli's
  // algorithm: it reaches a reduced cell for everything encountered here and
  // is twenty lines rather than two hundred. It does not guarantee the Niggli
  // setting, so it is a reduction and should not be described as more.
  Vec3 v[3] = {real_space_rows.row(0), real_space_rows.row(1),
               real_space_rows.row(2)};
  for (int pass = 0; pass < 100; ++pass) {
    bool changed = false;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        if (i == j)
          continue;
        const double denominator = v[j].norm_squared();
        if (denominator <= 0.0)
          continue;
        const double multiple = std::round(v[i].dot(v[j]) / denominator);
        if (multiple == 0.0)
          continue;
        const Vec3 shortened = v[i] - v[j] * multiple;
        if (shortened.norm_squared() < v[i].norm_squared() - 1e-9) {
          v[i] = shortened;
          changed = true;
        }
      }
    }
    // Keep them in increasing length, which makes the output stable between
    // runs and puts a and c where a reader expects them.
    std::sort(v, v + 3, [](const Vec3 &a, const Vec3 &b) {
      return a.norm_squared() < b.norm_squared();
    });
    if (!changed)
      break;
  }
  // Sign convention. A reduced cell should have all three angles acute or all
  // three obtuse -- Niggli's type I and type II -- and which of the two is
  // reachable is a property of the lattice, not a choice. Negating one vector
  // flips the sign of two of the three dot products, so the right combination
  // of signs always reaches it.
  //
  // Without this the cell comes out mixed, as 67.45 67.45 67.08 with angles
  // 109.5 70.7 70.9, which describes the same lattice as all-obtuse 109.45 but
  // compares equal to nothing and looks wrong to every reader.
  Mat3 best = Mat3::from_rows(v[0], v[1], v[2]);
  int best_score = -1;
  for (int mask = 0; mask < 8; ++mask) {
    Vec3 w[3] = {v[0], v[1], v[2]};
    for (int k = 0; k < 3; ++k) {
      if (mask & (1 << k))
        w[static_cast<std::size_t>(k)] = -w[static_cast<std::size_t>(k)];
    }
    const Mat3 candidate = Mat3::from_rows(w[0], w[1], w[2]);
    // A left-handed basis indexes as well as a right-handed one and is still
    // wrong, so it is not a candidate at all.
    if (candidate.determinant() <= 0.0)
      continue;
    const double bc = w[1].dot(w[2]);
    const double ac = w[0].dot(w[2]);
    const double ab = w[0].dot(w[1]);
    const int positive = (bc > 0) + (ac > 0) + (ab > 0);
    // Three of one kind is a proper reduced cell; two is worse than three and
    // better than a mixture with no majority.
    const int score = positive == 3 || positive == 0 ? 2 : 0;
    if (score > best_score) {
      best_score = score;
      best = candidate;
    }
  }
  return best;
}

bool choose_basis(const std::vector<Vec3> &candidates,
                  const std::vector<Vec3> &points, double tolerance,
                  Crystal *crystal, std::size_t *n_indexed) {
  const std::size_t n = candidates.size();
  if (n < 3)
    return false;

  std::size_t best_count = 0;
  double best_volume = 0.0;
  Mat3 best_rows = Mat3::identity();
  bool found = false;

  for (std::size_t a = 0; a < n; ++a) {
    for (std::size_t b = a + 1; b < n; ++b) {
      for (std::size_t c = b + 1; c < n; ++c) {
        const Mat3 rows =
            Mat3::from_rows(candidates[a], candidates[b], candidates[c]);
        const double volume = std::abs(rows.determinant());
        // Near-coplanar triples have a tiny volume and index everything
        // badly while scoring plausibly, so they are rejected on geometry
        // before they are scored.
        const double scale =
            candidates[a].norm() * candidates[b].norm() * candidates[c].norm();
        if (volume < 0.05 * scale) {
          ++g_triples_skipped;
          continue;
        }
        ++g_triples_scored;

        const std::size_t count = score_basis(rows, points, tolerance, nullptr);
        // More reflections indexed wins. On a tie the smaller cell wins,
        // because a supercell indexes everything its sublattice does and
        // would otherwise be chosen half the time by whichever came first.
        if (count > best_count ||
            (count == best_count && found && volume < best_volume - 1e-6)) {
          best_count = count;
          best_volume = volume;
          best_rows = rows;
          found = true;
        }
      }
    }
  }
  if (!found || best_count == 0)
    return false;

  const Mat3 reduced = reduce_basis(best_rows);
  // Reduction must not lose reflections. If it does, the reduction found a
  // different lattice and the unreduced basis is the safer answer.
  const std::size_t after = score_basis(reduced, points, tolerance, nullptr);
  const Mat3 chosen = after >= best_count ? reduced : best_rows;

  crystal->A = chosen.inverse();
  if (n_indexed)
    *n_indexed = std::max(best_count, after);
  return true;
}

namespace {

// Least squares fit of A to the reflections that currently index, solving
// r = A h for A. The normal equations are 3x3: A = (sum r h^T)(sum h h^T)^-1,
// so this costs one pass over the data and one small inverse.
//
// It matters more than it looks. The basis vectors out of the transform are
// peak centroids on a grid whose step is around 0.8 Angstrom, good to a few
// parts in a thousand, and a few parts in a thousand is a tenth of an index
// out at the edge of the detector. Fitting to the data is what turns a
// lattice that is recognisably right into one that is numerically right.
bool refit(const std::vector<Vec3> &points, double tolerance, Mat3 *A,
           std::size_t *n_used) {
  const Mat3 rows = A->inverse();
  double rh[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  double hh[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::size_t used = 0;

  for (const Vec3 &r : points) {
    const Vec3 f = rows * r;
    Vec3 h;
    double worst = 0.0;
    for (std::size_t k = 0; k < 3; ++k) {
      h[k] = std::round(f[k]);
      worst = std::fmax(worst, std::abs(f[k] - h[k]));
    }
    if (worst >= tolerance)
      continue;
    if (h.norm_squared() == 0.0)
      continue;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        rh[i * 3 + j] += r[i] * h[j];
        hh[i * 3 + j] += h[i] * h[j];
      }
    }
    ++used;
  }
  if (used < 10)
    return false;

  Mat3 sum_rh, sum_hh;
  for (std::size_t i = 0; i < 9; ++i) {
    sum_rh.m[i] = rh[i];
    sum_hh.m[i] = hh[i];
  }
  bool ok = false;
  const Mat3 inverse = sum_hh.inverse(&ok);
  // A rank-deficient h h^T means the indexed reflections lie in a plane, so
  // one direction of the cell is unconstrained. Refusing is the only honest
  // response; fitting anyway would produce a confident cell with an invented
  // third axis.
  if (!ok)
    return false;
  *A = sum_rh * inverse;
  if (n_used)
    *n_used = used;
  return true;
}

} // namespace

IndexResult index(ExperimentList &experiments, Table &reflections,
                  const IndexOptions &options) {
  IndexResult result;
  const double t_total = now_seconds();
  const double t_points = now_seconds();
  const std::vector<Vec3> points =
      reciprocal_lattice_points(experiments, reflections);
  result.timing.reciprocal_points = now_seconds() - t_points;
  result.n_total = points.size();
  if (points.size() < 10)
    return result;

  double d_min = options.d_min;
  if (d_min <= 0.0) {
    // The resolution of the data itself, taken at the 95th percentile of
    // |r| rather than the maximum: a few outliers at the corner of the
    // detector should not set the grid size for everything.
    std::vector<double> lengths;
    lengths.reserve(points.size());
    for (const Vec3 &r : points)
      lengths.push_back(r.norm());
    std::sort(lengths.begin(), lengths.end());
    const double q =
        lengths[static_cast<std::size_t>(0.95 * (lengths.size() - 1))];
    d_min = q > 0.0 ? 1.0 / q : 0.0;
  }
  if (!(d_min > 0.0))
    return result;

  const double t_max_cell = now_seconds();
  double max_cell = options.max_cell;
  if (max_cell <= 0.0) {
    // Grouped by sweep and turn, so that a reflection measured again on the
    // next rotation is not mistaken for a neighbour of itself.
    max_cell =
        estimate_max_cell(points, observation_groups(experiments, reflections));
  }
  if (!(max_cell > 0.0))
    return result;
  result.timing.max_cell = now_seconds() - t_max_cell;

  std::size_t grid = options.grid;
  if (grid == 0) {
    grid = next_power_of_two(
        static_cast<std::size_t>(std::ceil(4.0 * max_cell / d_min)));
    grid = std::min<std::size_t>(std::max<std::size_t>(grid, 64), 256);
  }

  result.d_min = d_min;
  result.max_cell = max_cell;
  result.grid = grid;

  const double t_candidates = now_seconds();
  result.candidates =
      find_candidate_vectors(points, options, d_min, max_cell, grid);
  result.timing.candidate_vectors = now_seconds() - t_candidates;
  result.timing.fft = g_last_fft_seconds;
  result.timing.peak_search = g_last_peak_seconds;
  if (options.verbose) {
    std::printf("  %zu candidate basis vectors\n", result.candidates.size());
    for (std::size_t i = 0;
         i < std::min<std::size_t>(6, result.candidates.size()); ++i) {
      std::printf("    |v| = %8.3f\n", result.candidates[i].norm());
    }
  }

  const double t_choose = now_seconds();
  g_triples_scored = 0;
  g_triples_skipped = 0;
  std::size_t n_indexed = 0;
  if (!choose_basis(result.candidates, points, options.tolerance,
                    &result.crystal, &n_indexed)) {
    return result;
  }
  result.timing.choose_basis = now_seconds() - t_choose;
  result.timing.triples_scored = g_triples_scored;
  result.timing.triples_skipped = g_triples_skipped;
  const double t_fit = now_seconds();

  // Fit, reduce, fit again. The first fit pulls the transform's peak centroids
  // onto the data; reduction then changes basis, which is exact and needs no
  // refitting, but the reassignment it implies can bring in reflections that
  // were outside tolerance before, so a second round is worth its cost.
  // Tighten as it converges. Fitting at the acceptance tolerance throughout
  // lets reflections that are a quarter of an index out -- misindexed, or
  // from a second lattice -- pull the basis with the same weight as ones that
  // are spot on. Each round rejects more, and because the basis improves each
  // time, the reflections that survive are the ones that belong.
  const double schedule[] = {options.tolerance, options.tolerance * 0.5,
                             options.tolerance * 0.25,
                             options.tolerance * 0.15};
  for (double tolerance : schedule) {
    std::size_t used = 0;
    Mat3 trial = result.crystal.A;
    if (!refit(points, tolerance, &trial, &used))
      break;
    result.crystal.A = trial;
    if (options.verbose) {
      std::printf("  fit at tol %.3f: %zu reflections\n", tolerance, used);
    }
  }
  if (options.verbose) {
    const UnitCell u = result.crystal.cell();
    std::printf("  after fit:   %.3f %.3f %.3f  %.2f %.2f %.2f  V=%.0f\n", u.a,
                u.b, u.c, u.alpha, u.beta, u.gamma, u.volume());
  }
  const Mat3 reduced = reduce_basis(result.crystal.A.inverse());
  if (options.verbose) {
    Crystal t;
    t.A = reduced.inverse();
    const UnitCell u = t.cell();
    std::printf("  after reduce:%.3f %.3f %.3f  %.2f %.2f %.2f  V=%.0f\n", u.a,
                u.b, u.c, u.alpha, u.beta, u.gamma, u.volume());
  }
  // Accept the reduction on the invariant that actually holds -- the cell
  // volume -- and not on the number of reflections indexed.
  //
  // The count is basis-dependent and therefore useless for this comparison.
  // Acceptance asks whether every fractional index is within tolerance of an
  // integer, and a unimodular change of basis mixes the three components: a
  // reflection sitting at (0.2, 0.2, 0.2) in one basis can be at (0.4, 0, 0.2)
  // in another that describes exactly the same lattice. An earlier version
  // guarded the reduction with this count and duly rejected a cell of
  // 67.11 67.47 67.47, 109.56 109.23 109.29 -- which is the right answer -- in
  // favour of the unreduced 67.47 77.92 77.87 it started from.
  const double volume_before =
      std::abs(result.crystal.A.inverse().determinant());
  const double volume_after = std::abs(reduced.determinant());
  if (volume_before > 0.0 &&
      std::abs(volume_after - volume_before) < 1e-6 * volume_before) {
    result.crystal.A = reduced.inverse();
  }
  for (double tolerance : schedule) {
    Mat3 trial = result.crystal.A;
    if (!refit(points, tolerance, &trial, nullptr))
      break;
    result.crystal.A = trial;
  }

  for (Experiment &e : experiments)
    e.crystal = result.crystal;

  // Assign indices under the current model. Returns how many took.
  const auto assign = [&](const std::vector<Vec3> &rlp) {
    Column &miller =
        reflections.int_column("miller_index", "cctbx::miller::index<>", 3);
    const Mat3 rows = result.crystal.A.inverse();
    double sum_squared = 0.0;
    std::size_t count = 0;
    for (std::size_t i = 0; i < rlp.size(); ++i) {
      const Vec3 h = rows * rlp[i];
      double worst = 0.0;
      double miss = 0.0;
      long hkl[3];
      for (std::size_t k = 0; k < 3; ++k) {
        const double rounded = std::round(h[k]);
        hkl[k] = static_cast<long>(rounded);
        const double d = h[k] - rounded;
        worst = std::fmax(worst, std::abs(d));
        miss += d * d;
      }
      if (worst >= options.tolerance || (!hkl[0] && !hkl[1] && !hkl[2]))
        continue;
      for (std::size_t k = 0; k < 3; ++k)
        miller.ints[i * 3 + k] = hkl[k];
      sum_squared += miss;
      ++count;
    }
    result.n_indexed = count;
    result.rmsd_index =
        count ? std::sqrt(sum_squared / static_cast<double>(count)) : 0.0;
  };

  assign(points);

  // Macrocycles: refine on the strong reflections, then assign to all again
  // under the improved model. The reciprocal lattice points move when the
  // detector does, so they are recomputed every cycle rather than reused.
  result.timing.fit_and_reduce = now_seconds() - t_fit;
  const double t_cycles = now_seconds();
  for (int cycle = 0; cycle < options.macrocycles; ++cycle) {
    const double t_subset = now_seconds();
    Table subset = reflections;
    std::size_t n_strong = reflections.nrows;
    if (options.refine_on_strong) {
      // Strength measured against this dataset's own median rather than an
      // absolute count, so it travels between detectors and spot finders.
      const char *column = reflections.has("n_signal") ? "n_signal"
                           : reflections.has("intensity.sum.value")
                               ? "intensity.sum.value"
                               : nullptr;
      if (column) {
        std::vector<double> values;
        for (std::size_t i = 0; i < reflections.nrows; ++i) {
          const Column &v = reflections.at(column);
          values.push_back(v.integral ? static_cast<double>(v.integer(i))
                                      : v.real(i));
        }
        std::vector<double> sorted = values;
        std::nth_element(sorted.begin(), sorted.begin() + sorted.size() / 2,
                         sorted.end());
        const double threshold = sorted[sorted.size() / 2];
        Column &m =
            subset.int_column("miller_index", "cctbx::miller::index<>", 3);
        const Column &original = reflections.at("miller_index");
        n_strong = 0;
        for (std::size_t i = 0; i < reflections.nrows; ++i) {
          if (values[i] < threshold)
            continue;
          for (std::size_t k = 0; k < 3; ++k) {
            m.ints[i * 3 + k] = original.integer(i, k);
          }
          if (original.integer(i, 0) || original.integer(i, 1) ||
              original.integer(i, 2)) {
            ++n_strong;
          }
        }
      }
    }

    result.timing.subset_copy += now_seconds() - t_subset;

    const double t_refine = now_seconds();
    RefineOptions refinement;
    refinement.outlier_sigma = 3.0;
    refinement.macrocycles = 2;
    refinement.verbose = false;
    // Analytical derivatives. They are validated against the finite
    // differences elsewhere and are the whole reason that code exists;
    // indexing had been using finite differences all along, which costs one
    // full residual evaluation per parameter per iteration instead of one
    // pass, and which no amount of threading the analytical path could help
    // because the analytical path was never called.
    refinement.analytic = true;
    const RefineResult r = refine(experiments, subset, refinement);
    result.timing.refine += now_seconds() - t_refine;
    if (r.n_used == 0)
      break;
    result.n_refined_on = n_strong;
    result.cycles_run = cycle + 1;
    if (experiments[0].crystal)
      result.crystal = *experiments[0].crystal;

    const double t_reassign = now_seconds();
    assign(reciprocal_lattice_points(experiments, reflections));
    result.timing.reassign += now_seconds() - t_reassign;
    IndexCycle record;
    record.refined_on = n_strong;
    record.rejected = r.n_rejected;
    record.indexed = result.n_indexed;
    record.rmsd_x = r.rmsd_x;
    record.rmsd_y = r.rmsd_y;
    record.rmsd_z = r.rmsd_z;
    record.rmsd_index = result.rmsd_index;
    result.cycles.push_back(record);
  }

  for (Experiment &e : experiments)
    e.crystal = result.crystal;
  result.timing.macrocycles = now_seconds() - t_cycles;
  result.timing.total = now_seconds() - t_total;
  return result;
}

} // namespace mxi
