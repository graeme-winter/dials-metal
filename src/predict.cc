#include "predict.hh"

#include <algorithm>
#include <atomic>
#include <thread>
#include <cmath>

namespace mxi {

namespace {

constexpr double kPi = 3.14159265358979323846;

// Put an angle into the half-open interval starting at `from`, so that a
// predicted phi can be compared against the scan range without either end
// needing a special case.
double wrap_from(double phi, double from) {
  const double two_pi = 2.0 * kPi;
  double d = std::fmod(phi - from, two_pi);
  if (d < 0.0) d += two_pi;
  return from + d;
}

}  // namespace

Intersections ewald_intersections(const Experiment &e, const Vec3 &r0) {
  Intersections out;

  const Vec3 m2 = e.goniometer.axis.normalized();
  // Work in the frame the rotation is about: move the incident beam through
  // the inverse setting rotation rather than rotating the axis, so the
  // decomposition below is about a fixed vector.
  const Vec3 s0p = e.goniometer.setting.transpose() * e.beam.s0();
  const Vec3 u = e.goniometer.fixed * r0;

  const double r_squared = u.norm_squared();
  if (r_squared <= 0.0) return out;
  // Beyond twice the Ewald radius nothing can reach the sphere at any angle.
  if (r_squared > 4.0 * e.beam.s0().norm_squared()) return out;

  const double u_parallel = u.dot(m2);
  const double s_parallel = s0p.dot(m2);
  const double a0 = u_parallel * s_parallel;
  const double b = s0p.dot(u) - a0;
  const double c = s0p.dot(m2.cross(u));
  const double d = -0.5 * r_squared - a0;

  const double amplitude = std::sqrt(b * b + c * c);
  if (amplitude <= 0.0) return out;
  const double cosine = d / amplitude;
  // A point whose component along the axis never brings it to the sphere.
  // The equality case is a tangential grazing, which is a measure-zero event
  // that acos handles without special-casing.
  if (cosine < -1.0 || cosine > 1.0) return out;

  const double centre = std::atan2(c, b);
  const double offset = std::acos(cosine);

  out.any = true;
  out.phi[0] = centre - offset;
  out.phi[1] = centre + offset;

  // Entering or exiting the sphere, from the sign of d/dphi of the distance to
  // the sphere centre. The derivative of s0 . R r at the crossing is
  // s0' . (m2 x R r), and the point is entering while that is negative.
  //
  // CONVENTION: this is a belief about DIALS' definition and is unvalidated;
  // if it is backwards, every `entering` flag is inverted and every keyed join
  // against DIALS output pairs the wrong two observations. See
  // docs/conventions.md.
  for (int i = 0; i < 2; ++i) {
    const Vec3 rotated = rotation(m2, out.phi[i]) * u;
    out.entering[i] = s0p.dot(m2.cross(rotated)) < 0.0;
  }
  return out;
}

std::array<int, 3> index_bounds(const Crystal &crystal, double d_min,
                                int max_index) {
  const Mat3 real = crystal.A.inverse();
  std::array<int, 3> bounds{};
  for (int i = 0; i < 3; ++i) {
    const double length = real.row(static_cast<std::size_t>(i)).norm();
    const double limit = d_min > 0.0 ? length / d_min : 0.0;
    bounds[static_cast<std::size_t>(i)] =
        std::min(max_index, static_cast<int>(std::ceil(limit)) + 1);
  }
  return bounds;
}

namespace {

// For a scan-varying crystal the setting matrix depends on where the
// reflection lands, and where it lands depends on the setting matrix. Solve it
// by iteration.
//
// The iteration must be seeded PER ROOT, with that root's own angle. A
// reflection crosses the Ewald sphere twice, at angles that can be most of the
// scan apart, and the crystal is a different shape at each. Seeding both from
// a single guess -- the start of the scan, say -- makes the iteration converge
// on whichever root is nearer that guess, so the second root is evaluated with
// the first one's setting matrix.
//
// That is not a subtle error. It made the forward and reverse maps disagree by
// 0.62 px on a half-degree drift, and because refinement's target uses the
// reverse map while the data came from the forward one, no amount of refining
// could reach the truth. Nothing caught it: the round-trip test covers static
// crystals, where both roots share one matrix and the bug is invisible.
struct Converged {
  bool any = false;
  Vec3 r0;
  double phi = 0.0;
  bool entering = false;
};

//: A root from ewald_intersections, moved into the scan's first turn.
//:
//: ewald_intersections answers in whatever 2 pi interval its arithmetic lands
//: in -- a leaving root of -259.9 degrees, say, for a reflection at 100.1 on a
//: scan from 0 to 180 -- and converge_root looks the crystal up at
//: z_from_phi of the angle it is given. At -259.9 degrees that is image -2600,
//: which setting_at clamps to image zero: the iteration then used the crystal
//: from the START of the scan for a reflection in the middle of it.
//:
//: That was invisible in the first half of a scan, where the crystal at image
//: zero is nearly the crystal at the reflection, and grew through the second.
//: On an 1800 image sweep, leaving reflections were exact to image 900 and then
//: up to 0.96 images late; entering reflections, whose roots happened to come
//: back inside the scan, were never affected. It had been there since the first
//: scan-varying version.
//:
//: Wrapped at the call rather than inside converge_root, which is also given
//: seeds already placed in later turns and must leave those where they are.
//:
//: WHICH 2 pi WINDOW MATTERS AT THE EDGES. The first version of this wrapped
//: into [start, start + 2 pi), and cured 34005 wrong predictions while making
//: 91 right ones wrong: a root a little BEFORE the scan start, -0.05 degrees on
//: a scan from 0, went to 359.95 -- the far end of the window, where the
//: crystal was looked up from the end of the scan. The same mistake as the one
//: being fixed, at the other end.
//:
//: So for a scan of less than a turn the window is centred on the scan's
//: middle, which leaves anything just outside either end on the side it came
//: from. A scan of a turn or more has no such middle to centre on, and
//: emit_turns expects the first turn, so it keeps [start, start + 2 pi): there
//: an angle just before the start lands at the end of the first turn, which is
//: a real position in the scan with a crystal model of its own.
double into_scan(const Experiment &e, double phi) {
  const double lo = std::min(e.scan.phi_start(), e.scan.phi_end());
  const double span = std::abs(e.scan.phi_end() - e.scan.phi_start());
  if (span >= 2.0 * kPi) return wrap_from(phi, lo);
  const double middle = lo + 0.5 * span;
  return wrap_from(phi, middle - kPi);
}

Converged converge_root(const Experiment &e, const Vec3 &h, double phi_seed,
                        bool entering_seed) {
  Converged out;
  out.phi = phi_seed;
  out.entering = entering_seed;
  out.r0 = e.crystal->A * h;
  if (!e.crystal->scan_varying()) {
    out.any = true;
    return out;
  }
  // ewald_intersections answers in a principal 2 pi interval whatever turn it
  // is asked about, and on a scan of several turns that was two bugs at once.
  //
  // The proximity test compared an angle in turn k against two roots in turn
  // zero, both about 2 pi k away, so which one was "nearer" depended on which
  // side of zero each fell -- and it could pick the OTHER crossing, flipping
  // the entering flag. And once one pass had assigned a principal angle, the
  // next evaluated the setting matrix at turn zero's frame rather than turn
  // k's, which for a crystal that moves is a different matrix.
  //
  // Together they put 218640 of 583040 predictions of a four turn sweep within
  // five frames of another observation of the same reflection with the same
  // flag: turns collapsing onto one another. The count was right and the
  // positions were not. A test using the same A at every scan point could not
  // see it, the geometry then being the static one; a crystal that actually
  // moves could.
  //
  // So each root is brought into the turn of the current estimate before it is
  // compared or kept, and the iteration stays in the turn it was seeded in.
  const double two_pi = 2.0 * kPi;
  const auto in_this_turn = [&](double root) {
    return root + two_pi * std::round((out.phi - root) / two_pi);
  };
  for (int pass = 0; pass < 3; ++pass) {
    const Vec3 trial = e.setting_at(e.scan.z_from_phi(out.phi)) * h;
    const Intersections cross = ewald_intersections(e, trial);
    if (!cross.any) return out;
    const double first = in_this_turn(cross.phi[0]);
    const double second = in_this_turn(cross.phi[1]);
    // Stay on the root we started from: the same side of the sphere, which is
    // what distinguishes the two crossings, and then the nearer in angle.
    int which;
    if (cross.entering[0] == entering_seed && cross.entering[1] != entering_seed) {
      which = 0;
    } else if (cross.entering[1] == entering_seed &&
               cross.entering[0] != entering_seed) {
      which = 1;
    } else {
      which = std::abs(first - out.phi) <= std::abs(second - out.phi) ? 0 : 1;
    }
    out.r0 = trial;
    out.phi = which == 0 ? first : second;
    out.entering = cross.entering[which];
  }
  out.any = true;
  return out;
}

// Turn one Ewald intersection into a Prediction, or reject it.
bool build(const Experiment &e, const PredictOptions &options, int h, int k,
           int l, const Vec3 &r0, double phi, bool entering,
           Prediction *out) {
  if (!options.allow_outside_scan) {
    const double start = e.scan.phi_start();
    const double end = e.scan.phi_end();
    const double lo = std::min(start, end);
    const double span = std::abs(end - start);
    // Wrapped into the scan's own interval, so a scan crossing 360 degrees and
    // a reflection predicted at -179 still meet.
    //
    // Only wrap a phi that is OUTSIDE the scan. wrap_from returns a value in
    // [lo, lo + 2 pi), so wrapping unconditionally collapses every turn of a
    // multi-turn sweep into the first one -- which it did: on ten rotations,
    // every shoebox landed in the first 3627 frames of 36000 and the same
    // pixels were integrated ten times over. A phi already in the scan is
    // already in the turn its caller meant it to be in.
    if (phi < lo || phi > lo + span) {
      const double wrapped = wrap_from(phi, lo);
      if (wrapped > lo + span) return false;
      phi = wrapped;
    }
  }

  const Mat3 r = e.goniometer.rotation_at(phi);
  const Vec3 s1 = e.beam.s0() + r * r0;

  auto hit = e.detector.intersect(s1);
  if (!hit) return false;

  out->h = h;
  out->k = k;
  out->l = l;
  out->entering = entering;
  out->phi = phi;
  out->panel = std::get<0>(*hit);
  out->px_fast = std::get<1>(*hit);
  out->px_slow = std::get<2>(*hit);
  out->z = e.scan.z_from_phi(phi);
  out->s1 = s1;
  return true;
}

}  // namespace

//: Emit one prediction per rotation of the scan.
//:
//: A reciprocal lattice point that crosses the Ewald sphere at phi crosses it
//: again at phi + 2 pi, so a sweep of ten full turns records every reflection
//: ten times over -- which is the entire reason for collecting one. Predicting
//: it once left a tenth of the reflections and, on a real ten rotation sweep,
//: 749786 predictions where DIALS made 7516507.
//:
//: Each turn's root is converged from its own seed rather than taken as
//: phi + 2 pi k, because with a scan-varying crystal the setting matrix at
//: turn nine is not the one at turn zero and the reflection is not in quite
//: the same place.
void emit_turns(const Experiment &e, const PredictOptions &options, int h, int k,
                int l, const Vec3 &hkl, const Converged &first, bool entering,
                std::vector<Prediction> *out) {
  const double lo = std::min(e.scan.phi_start(), e.scan.phi_end());
  const double span = std::abs(e.scan.phi_end() - e.scan.phi_start());
  const double two_pi = 2.0 * kPi;
  const int turns =
      options.allow_outside_scan
          ? 1
          : std::max(1, static_cast<int>(std::floor(span / two_pi)) + 1);
  const double base = wrap_from(first.phi, lo);
  for (int turn = 0; turn < turns; ++turn) {
    const double seed = base + two_pi * static_cast<double>(turn);
    if (!options.allow_outside_scan && seed > lo + span) break;
    // The first turn uses the root it was given. Re-converging it from the
    // wrapped seed is not the same calculation and lost nine reflections of a
    // thirty degree sweep, which is a scan of less than one turn and should
    // not have been touched at all by a change about many.
    //
    // Later turns converge from their own seed: converge_root returns the
    // crossing nearest it, so seeding a turn along is what puts the root in
    // that turn. For a scan-static crystal it comes back to the principal
    // value and the shift is put back by hand.
    Converged c = first;
    if (turn > 0) {
      c = converge_root(e, hkl, seed, entering);
      if (!c.any) continue;
      // converge_root answers with a crossing from ewald_intersections, which
      // is in a principal 2 pi interval whatever it was seeded with, so the
      // turn is lost on the way out. Put it back: the nearest turn to the
      // seed, which is the turn this iteration is for.
      //
      // Doing this only for a scan-static crystal is what hid the bug. The
      // synthetic test used one and passed; every real experiment has a
      // scan-varying crystal and every turn collapsed into the first.
      const double two_pi = 2.0 * kPi;
      c.phi += two_pi * std::round((seed - c.phi) / two_pi);
    }
    Prediction p;
    if (build(e, options, h, k, l, c.r0, c.phi, c.entering, &p)) {
      out->push_back(p);
    }
  }
}

std::vector<Prediction> predict_indices(
    const Experiment &e, const std::vector<std::array<int, 3>> &indices,
    const PredictOptions &options) {
  std::vector<Prediction> out;
  if (!e.crystal) return out;


  for (const std::array<int, 3> &hkl : indices) {
    if (hkl[0] == 0 && hkl[1] == 0 && hkl[2] == 0) continue;
    const Vec3 h{static_cast<double>(hkl[0]), static_cast<double>(hkl[1]),
                 static_cast<double>(hkl[2])};
    const Intersections cross = ewald_intersections(e, e.crystal->A * h);
    if (!cross.any) continue;
    for (int i = 0; i < 2; ++i) {
      const Converged c =
          converge_root(e, h, into_scan(e, cross.phi[i]), cross.entering[i]);
      if (!c.any) continue;
      emit_turns(e, options, hkl[0], hkl[1], hkl[2], h, c, cross.entering[i], &out);
    }
  }
  return out;
}

std::vector<Prediction> predict(const Experiment &e,
                                const PredictOptions &options) {
  std::vector<Prediction> out;
  if (!e.crystal) return out;

  double d_min = options.d_min;
  if (d_min <= 0.0) {
    // The absolute limit: a reciprocal lattice point further out than the
    // diameter of the Ewald sphere can never diffract.
    d_min = 0.5 * e.beam.wavelength;
  }
  const double q_max = 1.0 / d_min;
  const std::array<int, 3> bounds =
      index_bounds(*e.crystal, d_min, options.max_index);
  const Mat3 &A = e.crystal->A;

  // One independent unit of work per h, which on a large sweep is 39 per cent
  // of an integration and was all on one thread. Each h collects into its own
  // vector and they are joined in h order, so the output does not depend on
  // the thread count: everything downstream is indexed by position here, and a
  // prediction list that reorders itself would make every comparison between
  // two runs meaningless.
  std::size_t threads = options.threads;
  if (threads == 0) {
    threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;
  }
  const int first_h = -bounds[0];
  const std::size_t rows =
      static_cast<std::size_t>(2 * bounds[0] + 1);
  std::vector<std::vector<Prediction>> by_h(rows);

  const auto one_h = [&](std::size_t row) {
    const int h = first_h + static_cast<int>(row);
    std::vector<Prediction> &out = by_h[row];
    for (int k = -bounds[1]; k <= bounds[1]; ++k) {
      for (int l = -bounds[2]; l <= bounds[2]; ++l) {
        if (h == 0 && k == 0 && l == 0) continue;
        const Vec3 hkl{static_cast<double>(h), static_cast<double>(k),
                       static_cast<double>(l)};
        const Vec3 r0 = A * hkl;
        if (r0.norm() > q_max) continue;

        const Intersections cross = ewald_intersections(e, r0);
        if (!cross.any) continue;
        for (int i = 0; i < 2; ++i) {
          const Converged c =
              converge_root(e, hkl, into_scan(e, cross.phi[i]),
                            cross.entering[i]);
          if (!c.any) continue;
          emit_turns(e, options, h, k, l, hkl, c, cross.entering[i], &out);
        }
      }
    }
  };

  if (threads <= 1 || rows < 2) {
    for (std::size_t row = 0; row < rows; ++row) one_h(row);
  } else {
    std::atomic<std::size_t> next{0};
    std::vector<std::thread> pool;
    const std::size_t n = std::min(threads, rows);
    pool.reserve(n - 1);
    const auto run = [&]() {
      for (;;) {
        const std::size_t row = next.fetch_add(1);
        if (row >= rows) break;
        one_h(row);
      }
    };
    for (std::size_t t = 1; t < n; ++t) pool.emplace_back(run);
    run();
    for (std::thread &t : pool) t.join();
  }

  std::size_t total = 0;
  for (const std::vector<Prediction> &v : by_h) total += v.size();
  out.reserve(total);
  for (std::vector<Prediction> &v : by_h) {
    out.insert(out.end(), v.begin(), v.end());
    v.clear();
    v.shrink_to_fit();
  }
  return out;
}

}  // namespace mxi
