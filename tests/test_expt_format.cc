// The textual form of what we write, not just its meaning.
//
// dials.refine refused an indexed.expt from here with
//
//     DXTBX_ASSERT(obj_type == "float") failure   (dxtbx scan.cc:80)
//
// because a scan starting at zero degrees has an oscillation array whose first
// element is 0.0, the JSON writer collapsed any whole-valued double to an
// integer, and dxtbx reads the type of an array from its first element. The
// file was `[0, 0.1, 0.2, ...]` where dxtbx requires `[0.0, 0.1, 0.2, ...]`.
//
// Nothing in this package could have caught it. The document round-tripped
// through this reader perfectly, because this reader does not care: it parses
// 0 and 0.0 into the same double. **A round-trip test cannot detect a wrong
// convention. Only something written by someone else can** -- which is the
// lesson this repository has now learned four times, and the reason these
// assertions are about characters rather than values.

#include <cmath>
#include <string>

#include "../src/expt.h"
#include "../src/json.h"
#include "check.h"

using namespace mxi;

namespace {

// A minimal experiment whose scan starts at exactly zero, which is the case
// that failed and the case a default dials.import produces.
ExperimentList one_experiment(double osc_start) {
  Experiment e;
  e.beam.direction = {0.0, 0.0, 1.0};
  e.beam.wavelength = 1.0;
  Panel p;
  p.fast = {1.0, 0.0, 0.0};
  p.slow = {0.0, -1.0, 0.0};
  p.origin = {-100.0, 100.0, -200.0};
  p.pixel_size[0] = p.pixel_size[1] = 0.075;
  p.image_size[0] = 4148;
  p.image_size[1] = 4362;
  e.detector.panels.push_back(p);
  e.goniometer.axis = {1.0, 0.0, 0.0};
  e.scan.first_image = 1;
  e.scan.last_image = 10;
  e.scan.osc_start = osc_start;
  e.scan.osc_width = 0.1;
  e.crystal = Crystal::from_real_space({50.0, 0.0, 0.0}, {0.0, 60.0, 0.0},
                                       {0.0, 0.0, 70.0});
  ExperimentList list;
  list.experiments.push_back(e);
  return list;
}

// The text between the first '[' after `key` and its closing ']'.
std::string array_after(const std::string &text, const std::string &key) {
  const std::size_t at = text.find("\"" + key + "\"");
  if (at == std::string::npos) return "";
  const std::size_t open = text.find('[', at);
  if (open == std::string::npos) return "";
  const std::size_t close = text.find(']', open);
  if (close == std::string::npos) return "";
  return text.substr(open + 1, close - open - 1);
}

// Does every number in this fragment carry a decimal point or an exponent?
bool all_floating(const std::string &fragment) {
  std::size_t at = 0;
  while (at < fragment.size()) {
    while (at < fragment.size() && (std::isspace(static_cast<unsigned char>(fragment[at])) ||
                                    fragment[at] == ',')) {
      ++at;
    }
    const std::size_t start = at;
    while (at < fragment.size() && fragment[at] != ',') ++at;
    std::string token = fragment.substr(start, at - start);
    while (!token.empty() && std::isspace(static_cast<unsigned char>(token.back()))) {
      token.pop_back();
    }
    if (token.empty()) continue;
    if (token.find_first_of(".eE") == std::string::npos) return false;
  }
  return true;
}

}  // namespace

TEST(a_scan_starting_at_zero_writes_its_oscillation_as_floats) {
  // The exact failure. Before the fix the first element was `0`.
  const std::string text = json::dump(experiments_to_json(one_experiment(0.0)));
  const std::string oscillation = array_after(text, "oscillation");
  check::is_true(!oscillation.empty(), "there is an oscillation array");
  check::is_true(all_floating(oscillation),
                 "every oscillation value must be written as a float");
  // And specifically the first one, which is the element dxtbx types the
  // array from.
  check::is_true(oscillation.find("0.0") != std::string::npos,
                 "the leading zero must be 0.0 and not 0");
}

TEST(whole_numbers_that_are_measurements_still_write_as_floats) {
  // Not only the leading zero. A detector at exactly 200 mm, a wavelength of
  // exactly 1, an axis component of exactly 1 -- all are doubles that happen
  // to be whole, and all would have been written as integers.
  const std::string text = json::dump(experiments_to_json(one_experiment(0.0)));
  for (const char *key : {"origin", "fast_axis", "slow_axis", "pixel_size",
                          "rotation_axis", "real_space_a"}) {
    const std::string fragment = array_after(text, key);
    check::is_true(!fragment.empty(), std::string("found ") + key);
    check::is_true(all_floating(fragment),
                   std::string(key) + " must be written as floats");
  }
}

TEST(counts_and_ranges_still_write_as_integers) {
  // The other half. Making everything a float would be just as wrong: an
  // image_range of [1.0, 10.0] is not what dxtbx expects either, and the fix
  // for one must not break the other.
  const std::string text = json::dump(experiments_to_json(one_experiment(0.0)));
  for (const char *key : {"image_range", "image_size"}) {
    const std::string fragment = array_after(text, key);
    check::is_true(!fragment.empty(), std::string("found ") + key);
    check::is_true(!all_floating(fragment),
                   std::string(key) + " must be written as integers");
  }
}

TEST(a_whole_number_survives_the_json_writer_as_a_float) {
  // The unit beneath the three above. A double of exactly 5 must not come out
  // as `5`, and an int of 5 must not come out as `5.0`.
  json::Object object;
  object["measured"] = json::Value(5.0);
  object["counted"] = json::Value(5);
  const std::string text = json::dump(json::Value(std::move(object)));
  check::is_true(text.find("\"measured\": 5.0") != std::string::npos ||
                     text.find("\"measured\":5.0") != std::string::npos,
                 "a double writes with a decimal point");
  check::is_true(text.find("\"counted\": 5,") != std::string::npos ||
                     text.find("\"counted\":5,") != std::string::npos ||
                     text.find("\"counted\": 5\n") != std::string::npos ||
                     text.find("\"counted\": 5}") != std::string::npos,
                 "an integer writes without one");
}

TEST(the_value_still_reads_back_as_the_same_number) {
  // The formatting change must not cost precision. Seventeen significant
  // digits recovers any double, and appending .0 to a whole one does not move
  // it.
  for (double v : {0.0, 1.0, -1.0, 200.0, 0.1, 1e-17, 1.0 / 3.0, 6.02214076e23}) {
    json::Object object;
    object["x"] = json::Value(v);
    const json::Value back = json::parse(json::dump(json::Value(std::move(object))));
    check::close(back["x"].as_number(), v, 0.0, "exact round trip");
  }
}
