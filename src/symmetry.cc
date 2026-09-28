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
