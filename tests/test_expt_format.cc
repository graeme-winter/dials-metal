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

// --------------------------------------------------------------------------
// what we did not model is not ours to discard
// --------------------------------------------------------------------------

namespace {

// An experiment list as dials.import leaves one: the models this package reads,
// plus an imageset saying where the images are and three blocks it has no
// opinion about at all.
const char *kImported = R"({
  "__id__": "ExperimentList",
  "experiment": [{"__id__": "Experiment", "identifier": "abc",
                  "beam": 0, "detector": 0, "goniometer": 0, "scan": 0,
                  "imageset": 0, "crystal": null, "profile": null}],
  "imageset": [{"__id__": "ImageSequence",
                "template": "/data/ins10_1.nxs",
                "single_file_indices": [0, 1, 2],
                "params": {"dynamic_shadowing": "Auto"}}],
  "beam": [{"direction": [0.0, 0.0, 1.0], "wavelength": 1.0}],
  "detector": [{"panels": [{"fast_axis": [1.0, 0.0, 0.0],
                            "slow_axis": [0.0, -1.0, 0.0],
                            "origin": [-100.0, 100.0, -200.0],
                            "pixel_size": [0.075, 0.075],
                            "image_size": [4148, 4362],
                            "px_mm_strategy": {"type": "SimplePxMmStrategy"}}]}],
  "goniometer": [{"rotation_axis": [1.0, 0.0, 0.0],
                  "fixed_rotation": [1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],
                  "setting_rotation": [1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0]}],
  "scan": [{"image_range": [1, 3], "batch_offset": 0,
            "properties": {"oscillation": [0.0, 0.1, 0.2]}}],
  "crystal": [],
  "profile": [],
  "scaling_model": [],
  "history": {"dials.import": "1.2.3"}
})";

}  // namespace

TEST(blocks_we_do_not_model_survive_a_read_and_a_write) {
  // dials.refine failed on an indexed.expt from here with
  //
  //     IndexError: list index out of range   (experiment_list.py:603)
  //
  // because the experiment said "imageset": 0 and the imageset list had been
  // written out empty. The imageset is the only link from the file to the
  // images; profile, scaling_model and history are equally not ours to throw
  // away.
  const ExperimentList list = experiments_from_json(json::parse(kImported));
  check::equal(static_cast<long long>(list.size()), 1, "one experiment");

  const json::Value written = experiments_to_json(list);
  for (const char *block : {"imageset", "profile", "scaling_model", "history"}) {
    check::is_true(written.as_object().count(block) > 0,
                   std::string(block) + " must survive");
  }
  const json::Value &imageset = written["imageset"];
  check::is_true(imageset.is_array() && imageset.as_array().size() == 1,
                 "the imageset itself, not an empty list");
  check::is_true(imageset.as_array()[0]["template"].as_string() ==
                     "/data/ins10_1.nxs",
                 "and it still names the file");
}

TEST(every_model_reference_resolves_in_what_we_write) {
  // The check dxtbx makes and this failed: each index an experiment gives must
  // land inside the list it names. Absent models are null, never -1, because
  // -1 is a valid index into a Python list and reaches for the last element.
  const ExperimentList list = experiments_from_json(json::parse(kImported));
  const json::Value written = experiments_to_json(list);
  const json::Value &experiment = written["experiment"].as_array()[0];

  for (const char *model :
       {"beam", "detector", "goniometer", "scan", "crystal", "imageset"}) {
    const json::Value &reference = experiment[model];
    if (reference.is_null()) continue;
    check::is_true(reference.is_number(),
                   std::string(model) + " reference is a number or null");
    const double index = reference.as_number();
    check::is_true(index >= 0.0, std::string(model) + " index is not negative");
    check::is_true(written[model].is_array(),
                   std::string(model) + " list is present");
    check::is_true(index < static_cast<double>(written[model].as_array().size()),
                   std::string(model) + " index is inside its list");
  }
}

TEST(a_model_reference_is_written_as_an_integer) {
  // It survived as 0.0 once the block was carried through, because the reader
  // did not record that the document had written it without a decimal point.
  // dxtbx indexes a list with it.
  const ExperimentList list = experiments_from_json(json::parse(kImported));
  const std::string text = json::dump(experiments_to_json(list));
  check::is_true(text.find("\"imageset\": 0.0") == std::string::npos,
                 "no float model reference");
  check::is_true(text.find("\"imageset\": 0") != std::string::npos,
                 "an integer one");
}

TEST(an_experiment_list_built_in_memory_still_writes) {
  // No source document, so nothing to carry through: the models must still be
  // written and the absent ones must be null rather than missing or -1.
  const json::Value written = experiments_to_json(one_experiment(0.0));
  check::is_true(written["imageset"].is_array(), "an imageset list exists");
  check::is_true(written["imageset"].as_array().empty(), "and is empty");
  const json::Value &experiment = written["experiment"].as_array()[0];
  check::is_true(experiment["imageset"].is_null(),
                 "with the reference null, not zero and not -1");
}
