#include "expt.h"

#include <cmath>

namespace mxi {

namespace {

Vec3 vec3(const json::Value &v, const char *what) {
  double a[3];
  if (!v.numbers(3, a)) throw ExptError(std::string("expected three numbers for ") + what);
  return {a[0], a[1], a[2]};
}

json::Value to_json(const Vec3 &v) {
  return json::Array{json::Value(v.x), json::Value(v.y), json::Value(v.z)};
}

Beam read_beam(const json::Value &v) {
  Beam b;
  b.direction = vec3(v["direction"], "beam direction");
  b.wavelength = v["wavelength"].as_number();
  if (!(b.wavelength > 0.0)) throw ExptError("beam wavelength is not positive");
  return b;
}

Panel read_panel(const json::Value &v) {
  Panel p;
  p.name = v["name"].as_string();
  p.fast = vec3(v["fast_axis"], "panel fast axis");
  p.slow = vec3(v["slow_axis"], "panel slow axis");
  p.origin = vec3(v["origin"], "panel origin");
  if (!v["pixel_size"].numbers(2, p.pixel_size)) {
    throw ExptError("panel pixel_size must be two numbers");
  }
  const std::vector<double> size = v["image_size"].numbers();
  if (size.size() != 2) throw ExptError("panel image_size must be two numbers");
  p.image_size[0] = static_cast<std::int64_t>(size[0]);
  p.image_size[1] = static_cast<std::int64_t>(size[1]);

  const std::vector<double> trusted = v["trusted_range"].numbers();
  if (trusted.size() == 2) {
    p.trusted_min = trusted[0];
    p.trusted_max = trusted[1];
  }

  p.mu = v["mu"].as_number();
  p.thickness = v["thickness"].as_number();
  // The strategy decides it, not the presence of mu and thickness: a panel can
  // carry both and still be told not to apply the correction, and guessing
  // from the numbers would silently overrule the file.
  const std::string strategy = v["px_mm_strategy"]["type"].as_string();
  p.parallax = strategy == "ParallaxCorrectedPxMmStrategy" ||
               strategy == "OffsetParallaxCorrectedPxMmStrategy";
  if (p.parallax && !(p.mu > 0.0 && p.thickness > 0.0)) {
    throw ExptError(
        "panel asks for the parallax correction but has no mu or thickness");
  }
  return p;
}

Goniometer read_goniometer(const json::Value &v) {
  if (v.contains("axes")) {
    std::vector<Vec3> axes;
    for (const json::Value &a : v["axes"].as_array()) {
      axes.push_back(vec3(a, "goniometer axis"));
    }
    const std::vector<double> angles = v["angles"].numbers();
    if (angles.size() != axes.size()) {
      throw ExptError("goniometer has a different number of axes and angles");
    }
    const double scan_axis = v["scan_axis"].as_number(-1.0);
    if (scan_axis < 0.0 || scan_axis >= static_cast<double>(axes.size())) {
      throw ExptError("goniometer scan_axis is out of range");
    }
    return Goniometer::from_axes(axes, angles,
                                 static_cast<std::size_t>(scan_axis));
  }

  // The single-axis form.
  Goniometer g;
  g.axis = vec3(v["rotation_axis"], "rotation axis");
  const std::vector<double> fixed = v["fixed_rotation"].numbers();
  const std::vector<double> setting = v["setting_rotation"].numbers();
  if (fixed.size() == 9) {
    for (std::size_t i = 0; i < 9; ++i) g.fixed.m[i] = fixed[i];
  }
  if (setting.size() == 9) {
    for (std::size_t i = 0; i < 9; ++i) g.setting.m[i] = setting[i];
  }
  return g;
}

Scan read_scan(const json::Value &v) {
  const std::vector<double> range = v["image_range"].numbers();
  if (range.size() != 2) throw ExptError("scan image_range must be two numbers");
  const auto first = static_cast<std::int64_t>(range[0]);
  const auto last = static_cast<std::int64_t>(range[1]);

  if (v["properties"]["oscillation"].is_array()) {
    Scan s = Scan::from_oscillation(v["properties"]["oscillation"].numbers(),
                                    first, last);
    s.batch_offset = static_cast<std::int64_t>(v["batch_offset"].as_number());
    return s;
  }
  const std::vector<double> pair = v["oscillation"].numbers();
  if (pair.size() != 2) {
    throw ExptError(
        "scan has neither properties.oscillation nor a two-element "
        "oscillation; refusing to assume zero, which would compare equal to "
        "anything");
  }
  Scan s;
  s.first_image = first;
  s.last_image = last;
  s.osc_start = pair[0];
  s.osc_width = pair[1];
  s.batch_offset = static_cast<std::int64_t>(v["batch_offset"].as_number());
  return s;
}

Crystal read_crystal(const json::Value &v) {
  Crystal c = Crystal::from_real_space(vec3(v["real_space_a"], "real_space_a"),
                                       vec3(v["real_space_b"], "real_space_b"),
                                       vec3(v["real_space_c"], "real_space_c"));
  if (v.contains("space_group_hall_symbol")) {
    c.space_group_hall = v["space_group_hall_symbol"].as_string();
  }
  return c;
}

int model_index(const json::Value &e, const char *key) {
  const json::Value &v = e[key];
  if (!v.is_number()) return -1;
  const double d = v.as_number();
  return d < 0 ? -1 : static_cast<int>(d);
}

}  // namespace

ExperimentList experiments_from_json(const json::Value &document) {
  if (!document["experiment"].is_array()) {
    throw ExptError("no 'experiment' array at the top level");
  }

  std::vector<Beam> beams;
  for (const json::Value &v : document["beam"].as_array()) beams.push_back(read_beam(v));
  std::vector<Detector> detectors;
  for (const json::Value &v : document["detector"].as_array()) {
    Detector d;
    const json::Value &panels = v.is_array() ? v : v["panels"];
    for (const json::Value &p : panels.as_array()) d.panels.push_back(read_panel(p));
    if (d.panels.empty()) throw ExptError("detector has no panels");
    detectors.push_back(std::move(d));
  }
  std::vector<Goniometer> goniometers;
  for (const json::Value &v : document["goniometer"].as_array()) {
    goniometers.push_back(read_goniometer(v));
  }
  std::vector<Scan> scans;
  for (const json::Value &v : document["scan"].as_array()) scans.push_back(read_scan(v));
  std::vector<Crystal> crystals;
  for (const json::Value &v : document["crystal"].as_array()) {
    crystals.push_back(read_crystal(v));
  }

  ExperimentList list;
  for (const json::Value &e : document["experiment"].as_array()) {
    Experiment x;
    x.identifier = e["identifier"].as_string();

    const int b = model_index(e, "beam");
    const int d = model_index(e, "detector");
    const int g = model_index(e, "goniometer");
    const int s = model_index(e, "scan");
    const int c = model_index(e, "crystal");
    if (b >= 0 && static_cast<std::size_t>(b) < beams.size()) x.beam = beams[b];
    if (d >= 0 && static_cast<std::size_t>(d) < detectors.size()) x.detector = detectors[d];
    if (g >= 0 && static_cast<std::size_t>(g) < goniometers.size()) x.goniometer = goniometers[g];
    if (s >= 0 && static_cast<std::size_t>(s) < scans.size()) x.scan = scans[s];
    if (c >= 0 && static_cast<std::size_t>(c) < crystals.size()) x.crystal = crystals[c];
    list.experiments.push_back(std::move(x));
  }
  return list;
}

ExperimentList read_experiments(const std::string &path) {
  return experiments_from_json(json::parse_file(path));
}

namespace {

json::Value panel_to_json(const Panel &p) {
  json::Object o;
  o["name"] = json::Value(p.name.empty() ? std::string("panel") : p.name);
  o["type"] = json::Value("SENSOR_PAD");
  o["fast_axis"] = to_json(p.fast);
  o["slow_axis"] = to_json(p.slow);
  o["origin"] = to_json(p.origin);
  o["image_size"] = json::Array{json::Value(static_cast<long>(p.image_size[0])),
                                json::Value(static_cast<long>(p.image_size[1]))};
  o["pixel_size"] = json::Array{json::Value(p.pixel_size[0]),
                                json::Value(p.pixel_size[1])};
  o["trusted_range"] = json::Array{json::Value(p.trusted_min),
                                   json::Value(p.trusted_max)};
  o["thickness"] = json::Value(p.thickness);
  o["material"] = json::Value("Si");
  o["mu"] = json::Value(p.mu);
  o["identifier"] = json::Value("");
  o["mask"] = json::Value(json::Array{});
  o["gain"] = json::Value(1.0);
  o["pedestal"] = json::Value(0.0);
  json::Object strategy;
  strategy["type"] = json::Value(p.parallax ? "ParallaxCorrectedPxMmStrategy"
                                            : "SimplePxMmStrategy");
  o["px_mm_strategy"] = json::Value(std::move(strategy));
  return json::Value(std::move(o));
}

// Append `value` to `pool` unless an identical entry is already there, and
// return its index. This is what lets four sweeps of one crystal serialise as
// one crystal, matching how DIALS shares models.
std::size_t intern(std::vector<json::Value> &pool, json::Value value) {
  const std::string text = json::dump(value, 0);
  for (std::size_t i = 0; i < pool.size(); ++i) {
    if (json::dump(pool[i], 0) == text) return i;
  }
  pool.push_back(std::move(value));
  return pool.size() - 1;
}

}  // namespace

json::Value experiments_to_json(const ExperimentList &list) {
  std::vector<json::Value> beams, detectors, goniometers, scans, crystals;
  json::Array experiments;

  for (const Experiment &e : list) {
    json::Object x;
    x["__id__"] = json::Value("Experiment");
    x["identifier"] = json::Value(e.identifier);

    json::Object beam;
    beam["__id__"] = json::Value("monochromatic");
    beam["direction"] = to_json(e.beam.direction);
    beam["wavelength"] = json::Value(e.beam.wavelength);
    beam["polarization_normal"] = to_json({0.0, 1.0, 0.0});
    beam["polarization_fraction"] = json::Value(0.999);
    beam["probe"] = json::Value("x-ray");
    x["beam"] = json::Value(static_cast<long>(intern(beams, json::Value(std::move(beam)))));

    json::Array panels;
    for (const Panel &p : e.detector.panels) panels.push_back(panel_to_json(p));
    json::Object detector;
    detector["panels"] = json::Value(std::move(panels));
    x["detector"] = json::Value(
        static_cast<long>(intern(detectors, json::Value(std::move(detector)))));

    // The goniometer is written in the single-axis form, with the composed
    // fixed and setting rotations. That is lossless for everything downstream
    // -- those two matrices are all any calculation uses -- but it does not
    // round trip the individual axes and their names.
    json::Object gonio;
    gonio["rotation_axis"] = to_json(e.goniometer.axis);
    json::Array fixed, setting;
    for (std::size_t i = 0; i < 9; ++i) {
      fixed.push_back(json::Value(e.goniometer.fixed.m[i]));
      setting.push_back(json::Value(e.goniometer.setting.m[i]));
    }
    gonio["fixed_rotation"] = json::Value(std::move(fixed));
    gonio["setting_rotation"] = json::Value(std::move(setting));
    x["goniometer"] = json::Value(
        static_cast<long>(intern(goniometers, json::Value(std::move(gonio)))));

    json::Object scan;
    scan["image_range"] =
        json::Array{json::Value(static_cast<long>(e.scan.first_image)),
                    json::Value(static_cast<long>(e.scan.last_image))};
    scan["batch_offset"] = json::Value(static_cast<long>(e.scan.batch_offset));
    json::Array oscillation;
    for (std::int64_t i = 0; i < e.scan.num_images(); ++i) {
      oscillation.push_back(
          json::Value(e.scan.osc_start + static_cast<double>(i) * e.scan.osc_width));
    }
    json::Object properties;
    properties["oscillation"] = json::Value(std::move(oscillation));
    scan["properties"] = json::Value(std::move(properties));
    x["scan"] = json::Value(static_cast<long>(intern(scans, json::Value(std::move(scan)))));

    if (e.crystal) {
      json::Object crystal;
      crystal["__id__"] = json::Value("crystal");
      crystal["real_space_a"] = to_json(e.crystal->real_a());
      crystal["real_space_b"] = to_json(e.crystal->real_b());
      crystal["real_space_c"] = to_json(e.crystal->real_c());
      crystal["space_group_hall_symbol"] = json::Value(e.crystal->space_group_hall);
      x["crystal"] = json::Value(
          static_cast<long>(intern(crystals, json::Value(std::move(crystal)))));
    } else {
      x["crystal"] = json::Value(-1L);
    }
    x["imageset"] = json::Value(-1L);
    experiments.push_back(json::Value(std::move(x)));
  }

  json::Object out;
  out["__id__"] = json::Value("ExperimentList");
  out["experiment"] = json::Value(std::move(experiments));
  const auto pool = [](std::vector<json::Value> &v) {
    return json::Value(json::Array(v.begin(), v.end()));
  };
  out["beam"] = pool(beams);
  out["detector"] = pool(detectors);
  out["goniometer"] = pool(goniometers);
  out["scan"] = pool(scans);
  out["crystal"] = pool(crystals);
  out["imageset"] = json::Value(json::Array{});
  return json::Value(std::move(out));
}

void write_experiments(const std::string &path, const ExperimentList &list) {
  json::dump_file(path, experiments_to_json(list));
}

}  // namespace mxi
