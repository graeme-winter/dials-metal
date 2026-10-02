#include "nxmx_import.hh"

#include <hdf5.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <functional>
#include <random>
#include <stdexcept>
#include <utility>

#include "linalg.hh"

namespace mxi {

namespace {

constexpr double kPi = 3.14159265358979323846;

// The link information H5Literate gives its callback changed type in 1.12.
#if H5_VERSION_GE(1, 12, 0)
using LinkInfo = H5L_info2_t;
#else
using LinkInfo = H5L_info_t;
#endif
constexpr double kHcKeVA = 12.398419843320026; // h c, keV A

//: An HDF5 identifier closed with its owner.
class H5 {
public:
  H5(hid_t id, herr_t (*close)(hid_t)) : id_(id), close_(close) {}
  H5(const H5 &) = delete;
  H5 &operator=(const H5 &) = delete;
  ~H5() {
    if (id_ >= 0)
      close_(id_);
  }
  hid_t get() const { return id_; }
  bool ok() const { return id_ >= 0; }

private:
  hid_t id_;
  herr_t (*close_)(hid_t);
};

bool exists(hid_t file, const std::string &path) {
  if (path.empty() || path == "/")
    return true;
  // Each component in turn, as H5Lexists wants its parents to exist.
  std::string so_far;
  std::size_t at = path[0] == '/' ? 1 : 0;
  if (path[0] == '/')
    so_far = "";
  while (at <= path.size()) {
    const std::size_t next = path.find('/', at);
    const std::string part = path.substr(
        at, next == std::string::npos ? std::string::npos : next - at);
    so_far += "/" + part;
    if (H5Lexists(file, so_far.c_str(), H5P_DEFAULT) <= 0)
      return false;
    if (next == std::string::npos)
      break;
    at = next + 1;
  }
  return H5Oexists_by_name(file, path.c_str(), H5P_DEFAULT) > 0;
}

std::vector<double> read_doubles(hid_t file, const std::string &path) {
  H5 d(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
  if (!d.ok())
    throw std::runtime_error("cannot open " + path);
  H5 space(H5Dget_space(d.get()), H5Sclose);
  const hssize_t n = H5Sget_simple_extent_npoints(space.get());
  std::vector<double> out(static_cast<std::size_t>(n > 0 ? n : 0));
  if (n > 0 && H5Dread(d.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, out.data()) < 0)
    throw std::runtime_error("cannot read " + path + " as numbers");
  return out;
}

std::string read_string(hid_t file, const std::string &path) {
  H5 d(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
  if (!d.ok())
    throw std::runtime_error("cannot open " + path);
  H5 type(H5Dget_type(d.get()), H5Tclose);
  if (H5Tget_class(type.get()) != H5T_STRING)
    throw std::runtime_error(path + " is not a string");
  std::string out;
  if (H5Tis_variable_str(type.get()) > 0) {
    H5 mem(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(mem.get(), H5T_VARIABLE);
    H5Tset_cset(mem.get(), H5Tget_cset(type.get()));
    H5 space(H5Dget_space(d.get()), H5Sclose);
    if (H5Sget_simple_extent_npoints(space.get()) < 1)
      return out;
    std::vector<char *> values(
        static_cast<std::size_t>(H5Sget_simple_extent_npoints(space.get())),
        nullptr);
    if (H5Dread(d.get(), mem.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                values.data()) >= 0 &&
        values[0])
      out = values[0];
    H5Dvlen_reclaim(mem.get(), space.get(), H5P_DEFAULT, values.data());
  } else {
    const std::size_t size = H5Tget_size(type.get());
    std::vector<char> buffer(size + 1, '\0');
    H5 mem(H5Tcopy(type.get()), H5Tclose);
    if (H5Dread(d.get(), mem.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                buffer.data()) >= 0)
      out.assign(buffer.data(), strnlen(buffer.data(), size));
  }
  return out;
}

std::optional<std::string> attr_string(hid_t object, const char *name) {
  if (H5Aexists(object, name) <= 0)
    return std::nullopt;
  H5 a(H5Aopen(object, name, H5P_DEFAULT), H5Aclose);
  H5 type(H5Aget_type(a.get()), H5Tclose);
  if (H5Tget_class(type.get()) != H5T_STRING)
    return std::nullopt;
  std::string out;
  if (H5Tis_variable_str(type.get()) > 0) {
    H5 mem(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(mem.get(), H5T_VARIABLE);
    H5Tset_cset(mem.get(), H5Tget_cset(type.get()));
    char *value = nullptr;
    if (H5Aread(a.get(), mem.get(), &value) >= 0 && value) {
      out = value;
      H5free_memory(value);
    }
  } else {
    const std::size_t size = H5Tget_size(type.get());
    std::vector<char> buffer(size + 1, '\0');
    if (H5Aread(a.get(), type.get(), buffer.data()) >= 0)
      out.assign(buffer.data(), strnlen(buffer.data(), size));
  }
  return out;
}

std::optional<std::vector<double>> attr_doubles(hid_t object,
                                                const char *name) {
  if (H5Aexists(object, name) <= 0)
    return std::nullopt;
  H5 a(H5Aopen(object, name, H5P_DEFAULT), H5Aclose);
  H5 space(H5Aget_space(a.get()), H5Sclose);
  const hssize_t n = H5Sget_simple_extent_npoints(space.get());
  std::vector<double> out(static_cast<std::size_t>(n > 0 ? n : 0));
  if (n > 0 && H5Aread(a.get(), H5T_NATIVE_DOUBLE, out.data()) < 0)
    return std::nullopt;
  return out;
}

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

//: mm per unit of length; a unit not known is a failure, a missing one mm.
double to_mm(const std::optional<std::string> &units, const std::string &what,
             std::vector<std::string> *notes) {
  if (!units || units->empty()) {
    notes->push_back(what + " has no units: taken as mm");
    return 1.0;
  }
  const std::string u = lower(*units);
  if (u == "m" || u == "metre" || u == "meter" || u == "metres" ||
      u == "meters")
    return 1000.0;
  if (u == "mm" || u == "millimetre" || u == "millimeter" ||
      u == "millimetres" || u == "millimeters")
    return 1.0;
  if (u == "cm")
    return 10.0;
  if (u == "um" || u == "micron" || u == "microns" || u == "micrometre" ||
      u == "micrometer")
    return 1e-3;
  if (u == "nm")
    return 1e-6;
  throw std::runtime_error(what + " has units '" + *units +
                           "', which are not a length known here");
}

//: degrees per unit of angle.
double to_degrees(const std::optional<std::string> &units,
                  const std::string &what, std::vector<std::string> *notes) {
  if (!units || units->empty()) {
    notes->push_back(what + " has no units: taken as degrees");
    return 1.0;
  }
  const std::string u = lower(*units);
  if (u == "deg" || u == "degree" || u == "degrees")
    return 1.0;
  if (u == "rad" || u == "radian" || u == "radians")
    return 180.0 / kPi;
  throw std::runtime_error(what + " has units '" + *units +
                           "', which are not an angle known here");
}

//: One step of a NeXus transformation chain, in McStas mm and degrees.
struct Step {
  std::string path, name, type;
  Vec3 vector, offset;
  Vec3
      raw; //: the vector as the file gives it, which need not be of unit length
  std::vector<double> values; //: mm or degrees, one a frame where it moves
  std::string depends_on;
};

std::string parent_of(const std::string &path) {
  const std::size_t slash = path.rfind('/');
  return slash == 0 || slash == std::string::npos ? "/" : path.substr(0, slash);
}

std::string resolve(const std::string &from_group, const std::string &target) {
  if (target.empty() || target[0] == '/')
    return target;
  std::filesystem::path p = std::filesystem::path(from_group) / target;
  return p.lexically_normal().generic_string();
}

Step read_step(hid_t file, const std::string &path,
               std::vector<std::string> *notes) {
  H5 d(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
  if (!d.ok())
    throw std::runtime_error("the transformation " + path + " does not open");
  Step s;
  s.path = path;
  s.name = path.substr(path.rfind('/') + 1);
  s.type = lower(attr_string(d.get(), "transformation_type").value_or(""));
  const auto vec = attr_doubles(d.get(), "vector");
  if (!vec || vec->size() != 3)
    throw std::runtime_error(path + " has no vector");
  s.vector = Vec3{(*vec)[0], (*vec)[1], (*vec)[2]};
  s.raw = s.vector;
  if (s.vector.norm() > 0.0)
    s.vector = s.vector / s.vector.norm();
  const auto off = attr_doubles(d.get(), "offset");
  if (off && off->size() == 3) {
    const double k =
        to_mm(attr_string(d.get(), "offset_units"), path + "'s offset", notes);
    s.offset = Vec3{(*off)[0] * k, (*off)[1] * k, (*off)[2] * k};
  }
  s.depends_on = attr_string(d.get(), "depends_on").value_or(".");
  const std::vector<double> raw = read_doubles(file, path);
  const auto units = attr_string(d.get(), "units");
  double k = 1.0;
  if (s.type == "translation")
    k = to_mm(units, path, notes);
  else if (s.type == "rotation")
    k = to_degrees(units, path, notes);
  else
    throw std::runtime_error(path + " is neither a translation nor a rotation");
  for (double v : raw)
    s.values.push_back(v * k);
  if (s.values.empty())
    s.values.push_back(0.0);
  return s;
}

//: The chain from `start`, innermost first, to ".".
std::vector<Step> chain(hid_t file, const std::string &start,
                        std::vector<std::string> *notes) {
  std::vector<Step> out;
  std::string at = start;
  while (!at.empty() && at != ".") {
    if (out.size() > 64)
      throw std::runtime_error("the transformation chain from " + start +
                               " does not end");
    if (!exists(file, at))
      throw std::runtime_error("the transformation " + at + " does not exist");
    out.push_back(read_step(file, at, notes));
    const std::string next = out.back().depends_on;
    at = next == "." ? next : resolve(parent_of(at), next);
  }
  return out;
}

Mat3 rotation_deg(const Vec3 &axis, double degrees) {
  return rotation(axis, degrees * kPi / 180.0);
}

//: A point and the rotation applied to it by a chain, at each step's first
//: value: p -> R p + t, outermost last.
struct Pose {
  Mat3 R = Mat3::identity();
  Vec3 t;
};

Pose pose_of(const std::vector<Step> &steps) {
  Pose p;
  for (const Step &s : steps) { // innermost first: each wraps what came before
    Mat3 r = Mat3::identity();
    Vec3 shift = s.offset;
    if (s.type == "rotation")
      r = rotation_deg(s.vector, s.values[0]);
    else
      shift = shift + s.vector * s.values[0];
    p.R = r * p.R;
    p.t = r * p.t + shift;
  }
  return p;
}

//: McStas, which NeXus uses, to imgCIF, which DIALS uses: a half turn about y.
Vec3 to_imgcif(const Vec3 &v) { return Vec3{-v.x, v.y, -v.z}; }

//: A vector as JSON, -0 written as 0: the frame's half turn negates zeros.
json::Value vec(const Vec3 &v) {
  return json::Array{v.x + 0.0, v.y + 0.0, v.z + 0.0};
}

std::string material_symbol(const std::string &raw) {
  const std::string m = lower(raw);
  if (m == "si" || m == "silicon")
    return "Si";
  if (m == "cdte" || m == "cadmium telluride")
    return "CdTe";
  if (m == "gaas" || m == "gallium arsenide")
    return "GaAs";
  if (m == "ge" || m == "germanium")
    return "Ge";
  return raw;
}

std::string uuid4() {
  std::random_device rd;
  std::mt19937_64 gen((static_cast<std::uint64_t>(rd()) << 32) ^ rd());
  std::uniform_int_distribution<int> hex(0, 15);
  const char *digits = "0123456789abcdef";
  std::string s = "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx";
  for (char &c : s) {
    if (c == 'x')
      c = digits[hex(gen)];
    else if (c == 'y')
      c = digits[8 + hex(gen) % 4];
  }
  return s;
}

std::string now_utc() {
  const std::time_t t = std::time(nullptr);
  char buffer[32];
  std::strftime(buffer, sizeof buffer, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
  return buffer;
}

//: Every group under `root` whose NX_class is `cls`, by path.
std::vector<std::string> groups_of_class(hid_t file, const std::string &root,
                                         const std::string &cls) {
  std::vector<std::string> out;
  struct Ctx {
    hid_t file;
    std::string root, cls;
    std::vector<std::string> *out;
  } ctx{file, root, cls, &out};
  H5 g(H5Gopen2(file, root.c_str(), H5P_DEFAULT), H5Gclose);
  if (!g.ok())
    return out;
  H5Literate(
      g.get(), H5_INDEX_NAME, H5_ITER_INC, nullptr,
      [](hid_t group, const char *name, const LinkInfo *,
         void *data) -> herr_t {
        auto *c = static_cast<Ctx *>(data);
        // A group if it opens as one: the same in every HDF5 version, where
        // asking an object's type is not.
        H5 sub(H5Gopen2(group, name, H5P_DEFAULT), H5Gclose);
        if (sub.ok() &&
            attr_string(sub.get(), "NX_class").value_or("") == c->cls)
          c->out->push_back(c->root + "/" + name);
        return 0;
      },
      &ctx);
  return out;
}

std::optional<std::string>
first_existing(hid_t file, std::initializer_list<std::string> paths) {
  for (const std::string &p : paths)
    if (exists(file, p))
      return p;
  return std::nullopt;
}

} // namespace

double attenuation_coefficient(const std::string &material, double wavelength) {
  // Hubbell and Seltzer, NIST mass attenuation coefficients with coherent
  // scattering: energy keV, mu/rho cm^2/g, above the K edge for the MX range.
  struct Table {
    const char *name;
    double density; // g/cm^3
    std::vector<std::pair<double, double>> points;
  };
  static const Table silicon{"Si",
                             2.33,
                             {{2.0, 2777.0},
                              {3.0, 978.4},
                              {4.0, 452.9},
                              {5.0, 245.0},
                              {6.0, 147.0},
                              {8.0, 64.68},
                              {10.0, 33.89},
                              {15.0, 10.34},
                              {20.0, 4.464},
                              {30.0, 1.436},
                              {40.0, 0.7012},
                              {50.0, 0.4385},
                              {60.0, 0.3207},
                              {80.0, 0.2228},
                              {100.0, 0.1835}}};
  if (material != "Si" || !(wavelength > 0.0))
    return 0.0;
  const Table &t = silicon;
  const double energy = kHcKeVA / wavelength;
  if (energy < t.points.front().first || energy > t.points.back().first)
    return 0.0;
  for (std::size_t i = 0; i + 1 < t.points.size(); ++i) {
    const auto [e0, m0] = t.points[i];
    const auto [e1, m1] = t.points[i + 1];
    if (energy >= e0 && energy <= e1) {
      const double f =
          (std::log(energy) - std::log(e0)) / (std::log(e1) - std::log(e0));
      const double mu_rho =
          std::exp(std::log(m0) + f * (std::log(m1) - std::log(m0)));
      return mu_rho * t.density / 10.0; // 1/cm to 1/mm
    }
  }
  return 0.0;
}

json::Value import_nxmx(const std::string &master, const ImportOverrides &o,
                        std::vector<std::string> *notes) {
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  H5 file(H5Fopen(master.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
  if (!file.ok())
    throw std::runtime_error("cannot open " + master + " as HDF5");
  const hid_t f = file.get();

  // The beam.
  double wavelength = 0.0;
  if (o.wavelength) {
    wavelength = *o.wavelength;
    notes->push_back("the wavelength given, " + std::to_string(wavelength) +
                     " A");
  } else if (auto p = first_existing(
                 f, {"/entry/instrument/beam/incident_wavelength",
                     "/entry/sample/beam/incident_wavelength"})) {
    const std::vector<double> w = read_doubles(f, *p);
    if (w.empty())
      throw std::runtime_error(*p + " is empty");
    H5 d(H5Dopen2(f, p->c_str(), H5P_DEFAULT), H5Dclose);
    const std::string u =
        lower(attr_string(d.get(), "units").value_or("angstrom"));
    const double k = (u == "nm") ? 10.0 : (u == "m" ? 1e10 : 1.0);
    wavelength = w[0] * k;
    if (w.size() > 1 &&
        std::fabs(w.back() - w.front()) > 1e-9 * std::fabs(w.front()))
      notes->push_back(*p + " varies; the first value used");
  } else if (auto e =
                 first_existing(f, {"/entry/instrument/beam/incident_energy",
                                    "/entry/sample/beam/incident_energy"})) {
    const std::vector<double> energy = read_doubles(f, *e);
    H5 d(H5Dopen2(f, e->c_str(), H5P_DEFAULT), H5Dclose);
    const std::string u = lower(attr_string(d.get(), "units").value_or("eV"));
    const double kev = energy.at(0) * (u == "kev" ? 1.0 : 1e-3);
    wavelength = kHcKeVA / kev;
    notes->push_back("the wavelength from incident_energy");
  } else {
    throw std::runtime_error(
        "no incident_wavelength or incident_energy: give --wavelength");
  }
  json::Object beam{{"__id__", "monochromatic"},
                    {"direction", json::Array{0.0, 0.0, 1.0}},
                    {"wavelength", wavelength},
                    {"divergence", 0.0},
                    {"sigma_divergence", 0.0},
                    {"polarization_normal", json::Array{0.0, 1.0, 0.0}},
                    {"polarization_fraction", 0.999},
                    {"flux", 0.0},
                    {"transmission", 1.0},
                    {"probe", "x-ray"},
                    {"sample_to_source_distance", 0.0}};

  // The detector: a panel per module, its origin and axes through the
  // transformation chains, then in DIALS's frame.
  const std::string detector = "/entry/instrument/detector";
  if (!exists(f, detector))
    throw std::runtime_error("no " + detector);
  std::vector<std::string> modules =
      groups_of_class(f, detector, "NXdetector_module");
  if (modules.empty())
    throw std::runtime_error("no NXdetector_module under " + detector);
  std::string material = "Si";
  if (exists(f, detector + "/sensor_material"))
    material = material_symbol(read_string(f, detector + "/sensor_material"));
  double thickness = 0.0;
  if (exists(f, detector + "/sensor_thickness")) {
    H5 d(H5Dopen2(f, (detector + "/sensor_thickness").c_str(), H5P_DEFAULT),
         H5Dclose);
    thickness = read_doubles(f, detector + "/sensor_thickness").at(0) *
                to_mm(attr_string(d.get(), "units"), "sensor_thickness", notes);
  }
  double mu = o.mu ? *o.mu : attenuation_coefficient(material, wavelength);
  if (!o.mu && mu == 0.0)
    notes->push_back(
        "no attenuation coefficient tabulated for " + material +
        " at this wavelength: mu 0, so no parallax correction; give --mu");
  double trusted_max = 0.0;
  if (o.trusted_max) {
    trusted_max = *o.trusted_max;
  } else if (auto s = first_existing(
                 f,
                 {detector + "/saturation_value",
                  detector +
                      "/detectorSpecific/countrate_correction_count_cutoff"})) {
    trusted_max = read_doubles(f, *s).at(0);
  } else {
    int bits = 16;
    if (exists(f, detector + "/bit_depth_readout"))
      bits = static_cast<int>(
          read_doubles(f, detector + "/bit_depth_readout").at(0));
    else if (exists(f, detector + "/bit_depth_image"))
      bits = static_cast<int>(
          read_doubles(f, detector + "/bit_depth_image").at(0));
    trusted_max = std::ldexp(1.0, bits) - 2.0;
    notes->push_back("no saturation_value or count cutoff: the trusted range's "
                     "top taken as " +
                     std::to_string(static_cast<long long>(trusted_max)) +
                     ", just below the bad-pixel marker; give --trusted-max");
  }

  json::Array panels;
  for (const std::string &module : modules) {
    const std::vector<Step> fast_chain =
        chain(f, module + "/fast_pixel_direction", notes);
    const std::vector<Step> slow_chain =
        chain(f, module + "/slow_pixel_direction", notes);
    // The fast and slow directions' own steps give the pixel size and the
    // directions; everything they depend on places the module.
    const Step &fast = fast_chain.front();
    const Step &slow = slow_chain.front();
    const std::vector<Step> rest(fast_chain.begin() + 1, fast_chain.end());
    const Pose place = pose_of(rest);
    Vec3 origin = place.t + place.R * fast.offset;
    Vec3 fast_axis = place.R * fast.vector;
    const Pose slow_place =
        pose_of(std::vector<Step>(slow_chain.begin() + 1, slow_chain.end()));
    Vec3 slow_axis = slow_place.R * slow.vector;
    origin = to_imgcif(origin);
    fast_axis = to_imgcif(fast_axis);
    slow_axis = to_imgcif(slow_axis);
    const std::vector<double> size = read_doubles(f, module + "/data_size");
    const std::vector<double> start =
        exists(f, module + "/data_origin")
            ? read_doubles(f, module + "/data_origin")
            : std::vector<double>{0.0, 0.0};
    if (size.size() < 2)
      throw std::runtime_error(module + "/data_size is not two numbers");
    // NeXus gives slow then fast; DIALS fast then slow.
    const long long nx = static_cast<long long>(size[size.size() - 1]);
    const long long ny = static_cast<long long>(size[size.size() - 2]);
    const long long ox = static_cast<long long>(start[start.size() - 1]);
    const long long oy = static_cast<long long>(start[start.size() - 2]);
    panels.push_back(json::Object{
        {"name", module},
        {"type", "SENSOR_PAD"},
        {"fast_axis", vec(fast_axis)},
        {"slow_axis", vec(slow_axis)},
        {"origin", vec(origin)},
        {"raw_image_offset", json::Array{ox, oy}},
        {"image_size", json::Array{nx, ny}},
        {"pixel_size", json::Array{fast.values[0], slow.values[0]}},
        {"trusted_range", json::Array{0.0, trusted_max}},
        {"thickness", thickness},
        {"material", material},
        {"mu", mu},
        {"identifier", ""},
        {"mask", json::Array{}},
        {"gain", 1.0},
        {"pedestal", 0.0},
        {"px_mm_strategy",
         json::Object{{"type", mu > 0.0 && thickness > 0.0
                                   ? "ParallaxCorrectedPxMmStrategy"
                                   : "SimplePxMmStrategy"}}}});
  }
  if (panels.size() > 1)
    notes->push_back(std::to_string(panels.size()) +
                     " modules, so as many panels: a detector of several "
                     "panels is untested here");

  // The overrides of where the detector is, on the first panel's plane.
  auto &p0 = panels[0].as_object();
  const auto get = [](const json::Value &v) {
    const auto &a = v.as_array();
    return Vec3{a[0].as_number(), a[1].as_number(), a[2].as_number()};
  };
  if (o.distance || o.beam_centre) {
    const Vec3 fast_axis = get(p0["fast_axis"]),
               slow_axis = get(p0["slow_axis"]);
    Vec3 normal = fast_axis.cross(slow_axis);
    normal = normal / normal.norm();
    Vec3 origin = get(p0["origin"]);
    if (origin.dot(normal) < 0)
      normal = normal * -1.0;
    Vec3 shift;
    if (o.distance) {
      const double now = origin.dot(normal);
      shift = normal * (*o.distance - now);
      notes->push_back("the distance given, " + std::to_string(*o.distance) +
                       " mm");
    }
    Vec3 moved = origin + shift;
    if (o.beam_centre) {
      // Where the beam meets the panel's plane, and the origin put so that
      // that point is the pixel given.
      const Vec3 s0 = Vec3{0.0, 0.0, -1.0}; // the beam, towards the detector
      const double t = moved.dot(normal) / s0.dot(normal);
      const Vec3 hit = s0 * t;
      const auto px = p0["pixel_size"].as_array();
      moved = hit - fast_axis * ((*o.beam_centre)[0] * px[0].as_number()) -
              slow_axis * ((*o.beam_centre)[1] * px[1].as_number());
      notes->push_back("the beam centre given, pixel " +
                       std::to_string((*o.beam_centre)[0]) + ", " +
                       std::to_string((*o.beam_centre)[1]));
    }
    const Vec3 delta = moved - origin;
    for (json::Value &panel : panels)
      panel.as_object()["origin"] =
          vec(get(panel.as_object()["origin"]) + delta);
  }
  json::Object hierarchy{
      {"name", ""},
      {"type", ""},
      {"fast_axis", json::Array{1.0, 0.0, 0.0}},
      {"slow_axis", json::Array{0.0, 1.0, 0.0}},
      {"origin", json::Array{0.0, 0.0, 0.0}},
      {"raw_image_offset", json::Array{0, 0}},
      {"image_size", json::Array{0, 0}},
      {"pixel_size", json::Array{0.0, 0.0}},
      {"trusted_range", json::Array{0.0, 0.0}},
      {"thickness", 0.0},
      {"material", ""},
      {"mu", 0.0},
      {"identifier", ""},
      {"mask", json::Array{}},
      {"gain", 1.0},
      {"pedestal", 0.0},
      {"px_mm_strategy", json::Object{{"type", "SimplePxMmStrategy"}}}};
  json::Array children;
  for (std::size_t k = 0; k < panels.size(); ++k)
    children.push_back(json::Object{{"panel", static_cast<long long>(k)}});
  hierarchy["children"] = children;

  // The goniometer: the rotations of the sample's chain, innermost first; the
  // scan axis the one that moves.
  std::string sample_start;
  if (exists(f, "/entry/sample/depends_on"))
    sample_start =
        resolve("/entry/sample", read_string(f, "/entry/sample/depends_on"));
  if (sample_start.empty())
    throw std::runtime_error("/entry/sample has no depends_on: no goniometer");
  json::Array axes, angles, names;
  int scan_axis = -1;
  std::vector<double> scan_values;
  for (const Step &s : chain(f, sample_start, notes)) {
    if (s.type != "rotation")
      continue;
    const bool moves = s.values.size() > 1 &&
                       std::fabs(s.values.back() - s.values.front()) > 1e-12;
    if (moves) {
      if (scan_axis >= 0)
        throw std::runtime_error("more than one goniometer axis moves: " +
                                 s.name + " and the scan axis before it");
      scan_axis = static_cast<int>(axes.size());
      scan_values = s.values;
    }
    // As the file gives it, as dials.import writes it: normalised where used.
    axes.push_back(vec(to_imgcif(s.raw)));
    angles.push_back(moves ? 0.0 : s.values[0]);
    names.push_back(s.name);
  }
  if (scan_axis < 0)
    throw std::runtime_error("no goniometer axis moves: not a rotation scan");
  json::Object goniometer{{"axes", axes},
                          {"angles", angles},
                          {"names", names},
                          {"scan_axis", scan_axis}};

  // The scan.
  const long long frames = static_cast<long long>(scan_values.size());
  long long first = 1, last = frames;
  if (o.image_range) {
    first = (*o.image_range)[0];
    last = (*o.image_range)[1];
    if (first < 1 || last > frames || first > last)
      throw std::runtime_error("--image-range outside the " +
                               std::to_string(frames) + " images");
  }
  double exposure = 0.0;
  if (exists(f, detector + "/count_time"))
    exposure = read_doubles(f, detector + "/count_time").at(0);
  json::Array oscillation, epochs, exposures, indices;
  for (long long i = first - 1; i < last; ++i) {
    oscillation.push_back(scan_values[static_cast<std::size_t>(i)]);
    epochs.push_back(0.0);
    exposures.push_back(exposure);
    indices.push_back(i);
  }
  json::Object scan{{"image_range", json::Array{first, last}},
                    {"batch_offset", 0},
                    {"valid_image_ranges", json::Object{}},
                    {"properties", json::Object{{"epochs", epochs},
                                                {"exposure_time", exposures},
                                                {"oscillation", oscillation}}}};

  const std::string template_path =
      std::filesystem::absolute(master).lexically_normal().string();
  json::Object imageset{{"__id__", "ImageSequence"},
                        {"template", template_path},
                        {"single_file_indices", indices},
                        {"mask", json::Value()},
                        {"gain", json::Value()},
                        {"pedestal", json::Value()},
                        {"dx", json::Value()},
                        {"dy", json::Value()},
                        {"params", json::Object{{"dynamic_shadowing", "Auto"},
                                                {"multi_panel", false}}}};
  json::Object experiment{
      {"__id__", "Experiment"}, {"identifier", uuid4()}, {"beam", 0},
      {"detector", 0},          {"goniometer", 0},       {"scan", 0},
      {"imageset", 0}};
  return json::Object{
      {"__id__", "ExperimentList"},
      {"experiment", json::Array{experiment}},
      {"history", json::Array{now_utc() + "|mxi_import|0.1"}},
      {"imageset", json::Array{imageset}},
      {"beam", json::Array{beam}},
      {"detector",
       json::Array{json::Object{{"panels", panels}, {"hierarchy", hierarchy}}}},
      {"goniometer", json::Array{goniometer}},
      {"scan", json::Array{scan}},
      {"crystal", json::Array{}},
      {"profile", json::Array{}},
      {"scaling_model", json::Array{}}};
}

} // namespace mxi
