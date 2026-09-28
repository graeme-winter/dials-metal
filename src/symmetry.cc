#include "symmetry.hh"

#include <cctype>
#include <cmath>
#include <stdexcept>

#include <gemmi/symmetry.hpp>

namespace mxi {

struct SpaceGroup::Impl {
  const gemmi::SpaceGroup *group = nullptr;
  gemmi::GroupOps ops;
  gemmi::ReciprocalAsu asu;
  explicit Impl(const gemmi::SpaceGroup *g)
      : group(g), ops(g->operations()), asu(g) {}
};

SpaceGroup SpaceGroup::from_name(const std::string &name) {
  const gemmi::SpaceGroup *g = gemmi::find_spacegroup_by_name(name);
  if (g == nullptr)
    throw std::runtime_error("no space group called '" + name + "'");
  SpaceGroup out;
  out.impl_ = std::make_shared<const Impl>(g);
  return out;
}

SpaceGroup SpaceGroup::from_hall(const std::string &hall) {
  // An .expt's Hall symbol may carry a leading space, " P 1"; gemmi's table
  // may or may not match it verbatim, so the operators are compared instead.
  const gemmi::GroupOps ops = gemmi::symops_from_hall(hall.c_str());
  const gemmi::SpaceGroup *g = gemmi::find_spacegroup_by_ops(ops);
  if (g == nullptr)
    throw std::runtime_error("no tabulated space group has the Hall symbol '" +
                             hall + "'");
  SpaceGroup out;
  out.impl_ = std::make_shared<const Impl>(g);
  return out;
}

std::string SpaceGroup::name() const { return impl_->group->xhm(); }
std::string SpaceGroup::hall() const {
  return std::string(" ") + impl_->group->hall;
}
std::string SpaceGroup::laue() const { return impl_->group->laue_str(); }
int SpaceGroup::number() const { return impl_->group->number; }
std::size_t SpaceGroup::order() const { return impl_->ops.sym_ops.size(); }

Miller SpaceGroup::unique(const Miller &hkl) const {
  const auto a = impl_->asu.to_asu(hkl, impl_->ops);
  return {a.first[0], a.first[1], a.first[2]};
}

bool SpaceGroup::friedel_plus(const Miller &hkl) const {
  if (impl_->ops.is_reflection_centric(hkl))
    return true;
  return impl_->asu.to_asu(hkl, impl_->ops).second % 2 == 1;
}

bool SpaceGroup::absent(const Miller &hkl) const {
  return impl_->ops.is_systematically_absent(hkl);
}

bool SpaceGroup::centric(const Miller &hkl) const {
  return impl_->ops.is_reflection_centric(hkl);
}

namespace {

//: One axis of a change of basis, "b+c" or "-1/2a+1/2b+1/2c", as coefficients
//: of the old a, b and c.
std::array<double, 3> parse_axis(const std::string &text) {
  std::array<double, 3> out{0.0, 0.0, 0.0};
  std::size_t i = 0;
  bool any = false;
  while (i < text.size()) {
    while (i < text.size() && std::isspace(static_cast<unsigned char>(text[i])))
      ++i;
    if (i == text.size())
      break;
    double sign = 1.0;
    if (text[i] == '+' || text[i] == '-') {
      sign = text[i] == '-' ? -1.0 : 1.0;
      ++i;
    } else if (any) {
      throw std::runtime_error("expected + or - in '" + text + "'");
    }
    double coefficient = 1.0;
    const std::size_t start = i;
    while (i < text.size() &&
           (std::isdigit(static_cast<unsigned char>(text[i])) ||
            text[i] == '.' || text[i] == '/' || text[i] == '*'))
      ++i;
    if (i > start) {
      std::string number = text.substr(start, i - start);
      if (!number.empty() && number.back() == '*')
        number.pop_back();
      const std::size_t slash = number.find('/');
      coefficient = slash == std::string::npos
                        ? std::stod(number)
                        : std::stod(number.substr(0, slash)) /
                              std::stod(number.substr(slash + 1));
    }
    if (i == text.size() ||
        (text[i] != 'a' && text[i] != 'b' && text[i] != 'c'))
      throw std::runtime_error("expected a, b or c in '" + text + "'");
    out[static_cast<std::size_t>(text[i] - 'a')] += sign * coefficient;
    ++i;
    any = true;
  }
  if (!any)
    throw std::runtime_error("an empty axis in a change of basis");
  return out;
}

} // namespace

ChangeOfBasis ChangeOfBasis::parse(const std::string &text) {
  ChangeOfBasis cb;
  std::size_t start = 0;
  for (std::size_t j = 0; j < 3; ++j) {
    const std::size_t comma = text.find(',', start);
    if ((j < 2) == (comma == std::string::npos))
      throw std::runtime_error("a change of basis has three axes: '" + text +
                               "'");
    const std::array<double, 3> axis = parse_axis(
        text.substr(start, j < 2 ? comma - start : std::string::npos));
    for (std::size_t i = 0; i < 3; ++i)
      cb.M.m[i * 3 + j] = axis[i];
    start = comma + 1;
  }
  bool ok = false;
  (void)cb.M.inverse(&ok);
  if (!ok || std::abs(cb.M.determinant()) < 1e-9)
    throw std::runtime_error("'" + text +
                             "' is not a change of basis: singular");
  return cb;
}

Miller ChangeOfBasis::apply(const Miller &hkl) const {
  Miller out{};
  for (std::size_t j = 0; j < 3; ++j) {
    double v = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
      v += M.m[i * 3 + j] * hkl[i];
    const double r = std::round(v);
    if (std::abs(v - r) > 1e-6)
      throw std::runtime_error(
          "reindexed to a fractional index: the new cell is "
          "not a supercell of the old");
    out[j] = static_cast<int>(r);
  }
  return out;
}

void ChangeOfBasis::apply(Crystal &crystal) const {
  const Mat3 t = M.inverse().transpose();
  crystal.A = crystal.A * t;
  for (Mat3 &a : crystal.A_points)
    a = a * t;
}

void reindex(ExperimentList &experiments, Table &reflections,
             const ChangeOfBasis &cb, const SpaceGroup &group) {
  for (Experiment &e : experiments) {
    if (!e.crystal)
      continue;
    cb.apply(*e.crystal);
    e.crystal->space_group_hall = group.hall();
  }
  if (!reflections.has("miller_index"))
    return;
  Column miller = reflections.at("miller_index");
  for (std::size_t i = 0; i < reflections.nrows; ++i) {
    const Miller h{static_cast<int>(miller.ints[i * 3]),
                   static_cast<int>(miller.ints[i * 3 + 1]),
                   static_cast<int>(miller.ints[i * 3 + 2])};
    if (h[0] == 0 && h[1] == 0 && h[2] == 0)
      continue;
    const Miller n = cb.apply(h);
    for (std::size_t k = 0; k < 3; ++k)
      miller.ints[i * 3 + k] = n[k];
  }
  reflections.set("miller_index", std::move(miller));
}

} // namespace mxi

// ---- lattice symmetry --------------------------------------------------

#include <algorithm>
#include <set>

#include <gemmi/twin.hpp>

namespace mxi {

namespace {

Rotation from_gemmi(const gemmi::Op::Rot &r) {
  Rotation out{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out[static_cast<std::size_t>(i * 3 + j)] = r[i][j] / gemmi::Op::DEN;
  return out;
}

Rotation multiply(const Rotation &a, const Rotation &b) {
  Rotation out{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      int v = 0;
      for (int k = 0; k < 3; ++k)
        v += a[static_cast<std::size_t>(i * 3 + k)] *
             b[static_cast<std::size_t>(k * 3 + j)];
      out[static_cast<std::size_t>(i * 3 + j)] = v;
    }
  return out;
}

const Rotation kIdentity{1, 0, 0, 0, 1, 0, 0, 0, 1};

std::vector<Rotation> closure(const std::vector<Rotation> &generators) {
  std::set<Rotation> group{kIdentity};
  std::vector<Rotation> frontier{kIdentity};
  while (!frontier.empty()) {
    std::vector<Rotation> next;
    for (const Rotation &a : frontier)
      for (const Rotation &g : generators) {
        const Rotation p = multiply(a, g);
        if (group.insert(p).second)
          next.push_back(p);
      }
    frontier.swap(next);
  }
  std::vector<Rotation> out(group.begin(), group.end());
  std::stable_partition(out.begin(), out.end(),
                        [](const Rotation &r) { return r == kIdentity; });
  return out;
}

gemmi::UnitCell to_gemmi(const UnitCell &c) {
  return gemmi::UnitCell(c.a, c.b, c.c, c.alpha, c.beta, c.gamma);
}

std::string axis_text(const std::array<int, 3> &v) {
  std::string out;
  const char *names = "abc";
  for (int i = 0; i < 3; ++i) {
    const int c = v[static_cast<std::size_t>(i)];
    if (c == 0)
      continue;
    if (c < 0)
      out += "-";
    else if (!out.empty())
      out += "+";
    if (std::abs(c) != 1)
      out += std::to_string(std::abs(c)) + "*";
    out += names[i];
  }
  return out;
}

} // namespace

Miller apply(const Rotation &r, const Miller &h) {
  Miller out{};
  for (int j = 0; j < 3; ++j)
    out[static_cast<std::size_t>(j)] =
        r[static_cast<std::size_t>(0 * 3 + j)] * h[0] +
        r[static_cast<std::size_t>(1 * 3 + j)] * h[1] +
        r[static_cast<std::size_t>(2 * 3 + j)] * h[2];
  return out;
}

std::vector<Rotation> lattice_symmetry(const UnitCell &cell, double max_delta) {
  const gemmi::GroupOps ops =
      gemmi::find_lattice_symmetry(to_gemmi(cell), 'P', max_delta);
  std::vector<Rotation> out;
  for (const gemmi::Op &op : ops.sym_ops)
    out.push_back(from_gemmi(op.rot));
  std::stable_partition(out.begin(), out.end(),
                        [](const Rotation &r) { return r == kIdentity; });
  return out;
}

std::vector<std::vector<Rotation>>
subgroups(const std::vector<Rotation> &group) {
  std::set<std::vector<Rotation>> seen;
  for (std::size_t i = 0; i < group.size(); ++i)
    for (std::size_t j = i; j < group.size(); ++j) {
      std::vector<Rotation> g = closure({group[i], group[j]});
      std::vector<Rotation> key = g;
      std::sort(key.begin(), key.end());
      seen.insert(key);
    }
  std::vector<std::vector<Rotation>> out;
  for (const std::vector<Rotation> &key : seen) {
    std::vector<Rotation> g = key;
    std::stable_partition(g.begin(), g.end(),
                          [](const Rotation &r) { return r == kIdentity; });
    out.push_back(std::move(g));
  }
  std::sort(out.begin(), out.end(),
            [](const auto &a, const auto &b) { return a.size() > b.size(); });
  return out;
}

std::optional<Setting> reference_setting(const std::vector<Rotation> &rotations,
                                         const UnitCell &cell) {
  // The Patterson group in the given (primitive) basis: rotations, inversion.
  gemmi::GroupOps base;
  base.cen_ops.push_back({0, 0, 0});
  for (const Rotation &r : rotations) {
    gemmi::Op op = gemmi::Op::identity();
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        op.rot[i][j] = r[static_cast<std::size_t>(i * 3 + j)] * gemmi::Op::DEN;
    base.sym_ops.push_back(op);
    gemmi::Op inverted = op;
    inverted.rot = op.negated_rot();
    base.sym_ops.push_back(inverted);
  }
  // Direct-space lattice vectors, components -2..2, where a conventional axis
  // can lie: along a rotation axis, which some rotation fixes. A group with
  // a single axis -- cyclic -- has its conventional a and b perpendicular to
  // it, where no rotation of the group fixes them: a two-fold reverses them,
  // and for a three- or six-fold they are found by the cell's metric.
  const gemmi::UnitCell metric = to_gemmi(cell);
  const auto cartesian = [&](const std::array<int, 3> &v) {
    return metric.orth.mat.multiply(gemmi::Vec3(v[0], v[1], v[2]));
  };
  std::vector<std::array<int, 3>> all;
  for (int x = -2; x <= 2; ++x)
    for (int y = -2; y <= 2; ++y)
      for (int z = -2; z <= 2; ++z)
        if (x != 0 || y != 0 || z != 0)
          all.push_back({x, y, z});
  const auto image = [](const Rotation &r, const std::array<int, 3> &v) {
    std::array<int, 3> out{};
    for (int i = 0; i < 3; ++i)
      out[static_cast<std::size_t>(i)] =
          r[static_cast<std::size_t>(i * 3)] * v[0] +
          r[static_cast<std::size_t>(i * 3 + 1)] * v[1] +
          r[static_cast<std::size_t>(i * 3 + 2)] * v[2];
    return out;
  };
  std::vector<std::array<int, 3>> candidates, fixed_axes;
  for (const auto &v : all)
    for (const Rotation &r : rotations)
      if (r != kIdentity && image(r, v) == v) {
        candidates.push_back(v);
        fixed_axes.push_back(v);
        break;
      }
  // Do the fixed vectors span space? If not, the group has one axis.
  bool spans = false;
  for (std::size_t i = 0; i < fixed_axes.size() && !spans; ++i)
    for (std::size_t j = i + 1; j < fixed_axes.size() && !spans; ++j) {
      const gemmi::Vec3 a = cartesian(fixed_axes[i]),
                        b = cartesian(fixed_axes[j]);
      if (a.cross(b).length() > 1e-6 * a.length() * b.length())
        spans = true;
    }
  if (rotations.size() <= 1) {
    candidates = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; // P -1: the cell as it is
  } else if (!spans && !fixed_axes.empty()) {
    const gemmi::Vec3 axis = cartesian(fixed_axes.front());
    for (const auto &v : all) {
      bool reversed = false;
      for (const Rotation &r : rotations) {
        const auto rv = image(r, v);
        if (rv[0] == -v[0] && rv[1] == -v[1] && rv[2] == -v[2])
          reversed = true;
      }
      const gemmi::Vec3 c = cartesian(v);
      const bool perpendicular =
          std::abs(c.dot(axis)) < 0.035 * c.length() * axis.length();
      if (reversed || perpendicular)
        candidates.push_back(v);
    }
  }

  const std::size_t want = base.sym_ops.size();
  std::optional<Setting> best;
  int best_cost = 1 << 30;
  double best_beta = 1e9;
  int best_negatives = 1 << 30;
  for (const auto &u : candidates)
    for (const auto &v : candidates)
      for (const auto &w : candidates) {
        const int det = u[0] * (v[1] * w[2] - v[2] * w[1]) -
                        v[0] * (u[1] * w[2] - u[2] * w[1]) +
                        w[0] * (u[1] * v[2] - u[2] * v[1]);
        if (det <= 0 || det > 4)
          continue;
        int cost = 0;
        for (int i = 0; i < 3; ++i)
          cost += std::abs(u[static_cast<std::size_t>(i)]) +
                  std::abs(v[static_cast<std::size_t>(i)]) +
                  std::abs(w[static_cast<std::size_t>(i)]);
        int negatives = 0;
        for (int i = 0; i < 3; ++i)
          negatives += (u[static_cast<std::size_t>(i)] < 0) +
                       (v[static_cast<std::size_t>(i)] < 0) +
                       (w[static_cast<std::size_t>(i)] < 0);
        if (cost > best_cost)
          continue;
        // M: columns the new axes. The operators move by cob = M^-1.
        gemmi::Op m = gemmi::Op::identity();
        for (int i = 0; i < 3; ++i) {
          m.rot[i][0] = u[static_cast<std::size_t>(i)] * gemmi::Op::DEN;
          m.rot[i][1] = v[static_cast<std::size_t>(i)] * gemmi::Op::DEN;
          m.rot[i][2] = w[static_cast<std::size_t>(i)] * gemmi::Op::DEN;
        }
        gemmi::GroupOps ops = base;
        ops.change_basis_backward(m);
        bool integral = true;
        for (const gemmi::Op &op : ops.sym_ops)
          for (int i = 0; i < 3 && integral; ++i)
            for (int j = 0; j < 3; ++j)
              if (op.rot[i][j] % gemmi::Op::DEN != 0)
                integral = false;
        if (!integral)
          continue;
        const char centring = ops.find_centering();
        const gemmi::SpaceGroup *found = nullptr;
        for (const gemmi::SpaceGroup &sg : gemmi::spacegroup_tables::main) {
          if (!sg.is_reference_setting() || sg.hall[0] != '-' ||
              sg.hall[1] != centring)
            continue;
          const gemmi::GroupOps ref = sg.operations();
          if (ref.sym_ops.size() * ref.cen_ops.size() !=
              want * ops.cen_ops.size())
            continue;
          if (ops.is_same_as(ref)) {
            found = &sg;
            break;
          }
        }
        if (!found)
          continue;
        const std::string text =
            axis_text(u) + "," + axis_text(v) + "," + axis_text(w);
        const ChangeOfBasis cb = ChangeOfBasis::parse(text);
        Crystal c;
        // Only the metric matters here: beta of the new cell, for monoclinic.
        const gemmi::UnitCell g = to_gemmi(cell);
        c.A = Mat3{g.frac.mat[0][0], g.frac.mat[1][0], g.frac.mat[2][0],
                   g.frac.mat[0][1], g.frac.mat[1][1], g.frac.mat[2][1],
                   g.frac.mat[0][2], g.frac.mat[1][2], g.frac.mat[2][2]};
        cb.apply(c);
        // beta nearest 90 decides only for monoclinic, one two-fold; for any
        // other group it chose a permutation of a triclinic cell's axes.
        const double beta =
            rotations.size() == 2 ? std::abs(c.cell().beta - 90.0) : 0.0;
        // Fewest coefficients; then, for monoclinic, beta nearest 90; then the
        // fewest negative signs, so that the operator reads as DIALS' do.
        const bool better =
            cost < best_cost ||
            (cost == best_cost &&
             (beta < best_beta - 1e-6 || (std::abs(beta - best_beta) <= 1e-6 &&
                                          negatives < best_negatives)));
        if (better) {
          best_cost = cost;
          best_beta = beta;
          best_negatives = negatives;
          Setting s{SpaceGroup::from_hall(std::string(" ") + found->hall), cb,
                    text};
          best = std::move(s);
        }
      }
  return best;
}

UnitCell reindexed_cell(const UnitCell &cell, const ChangeOfBasis &cb) {
  const gemmi::UnitCell g = to_gemmi(cell);
  gemmi::Vec3 axes[3];
  for (int j = 0; j < 3; ++j)
    axes[j] = g.orth.mat.multiply(
        gemmi::Vec3(cb.M.m[static_cast<std::size_t>(0 * 3 + j)],
                    cb.M.m[static_cast<std::size_t>(1 * 3 + j)],
                    cb.M.m[static_cast<std::size_t>(2 * 3 + j)]));
  const auto angle = [](const gemmi::Vec3 &a, const gemmi::Vec3 &b) {
    return std::acos(std::fmax(-1.0, std::fmin(1.0, a.cos_angle(b)))) * 180.0 /
           std::acos(-1.0);
  };
  return UnitCell{axes[0].length(),        axes[1].length(),
                  axes[2].length(),        angle(axes[1], axes[2]),
                  angle(axes[0], axes[2]), angle(axes[0], axes[1])};
}

double d_spacing(const UnitCell &cell, const Miller &hkl) {
  return to_gemmi(cell).calculate_d(gemmi::Miller{hkl[0], hkl[1], hkl[2]});
}

std::vector<SpaceGroup>
space_groups_with_patterson(const SpaceGroup &patterson) {
  gemmi::GroupOps target = gemmi::symops_from_hall(patterson.hall().c_str());
  std::vector<SpaceGroup> out;
  for (const gemmi::SpaceGroup &sg : gemmi::spacegroup_tables::main) {
    if (!sg.is_reference_setting())
      continue;
    // Macromolecules are chiral, so their crystals are in groups of proper
    // rotations only (the Sohncke groups), as dials.symmetry assumes: a
    // centrosymmetric group is its own Patterson group, and would otherwise
    // be offered for every crystal of a centric Laue class.
    const gemmi::GroupOps full = sg.operations();
    bool proper = true;
    for (const gemmi::Op &op : full.sym_ops)
      if (op.det_rot() < 0)
        proper = false;
    if (!proper)
      continue;
    gemmi::GroupOps ops = full.derive_symmorphic();
    ops.add_inversion();
    if (ops.is_same_as(target))
      out.push_back(SpaceGroup::from_hall(std::string(" ") + sg.hall));
  }
  std::stable_sort(out.begin(), out.end(),
                   [](const SpaceGroup &a, const SpaceGroup &b) {
                     return a.number() < b.number();
                   });
  return out;
}

} // namespace mxi
