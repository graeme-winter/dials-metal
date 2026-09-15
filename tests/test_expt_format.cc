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

// --------------------------------------------------------------------------
// the columns dials.index produces
// --------------------------------------------------------------------------

#include "../src/derivatives.h"
#include "../src/refine.h"

namespace {

Table observations_of(const ExperimentList &list, std::size_t n) {
  Table t;
  t.nrows = n;
  Column &xyz = t.real_column("xyzobs.px.value", "vec3<double>", 3);
  Column &var = t.real_column("xyzobs.px.variance", "vec3<double>", 3);
  Column &panel = t.int_column("panel", "std::size_t", 1);
  Column &id = t.int_column("id", "int", 1);
  for (std::size_t i = 0; i < n; ++i) {
    xyz.reals[i * 3 + 0] = 300.0 + 220.0 * static_cast<double>(i);
    xyz.reals[i * 3 + 1] = 250.0 + 180.0 * static_cast<double>(i);
    xyz.reals[i * 3 + 2] = 0.5 + static_cast<double>(i);
    var.reals[i * 3 + 0] = var.reals[i * 3 + 1] = var.reals[i * 3 + 2] = 1.0 / 12.0;
    panel.ints[i] = 0;
    id.ints[i] = 0;
  }
  (void)list;
  return t;
}

}  // namespace

TEST(millimetre_centroids_are_the_corrected_conversion_and_an_angle) {
  // dials.refine will not read a table without xyzobs.mm.value: DIALS measures
  // its residual in millimetres and radians, so this is the observation it
  // minimises against.
  const ExperimentList list = one_experiment(0.0);
  Table t = observations_of(list, 5);
  add_observed_columns(list, t);

  const Panel &p = list[0].detector[0];
  const Column &px = t.at("xyzobs.px.value");
  const Column &mm = t.at("xyzobs.mm.value");
  for (std::size_t i = 0; i < t.nrows; ++i) {
    const auto expected = p.px_to_mm(px.real(i, 0), px.real(i, 1));
    check::close(mm.real(i, 0), expected.first, 1e-12, "fast in millimetres");
    check::close(mm.real(i, 1), expected.second, 1e-12, "slow in millimetres");
    // Radians, not degrees and not images.
    check::close(mm.real(i, 2), list[0].scan.phi_from_z(px.real(i, 2)), 1e-12,
                 "the rotation angle in radians");
  }

  // The variance carries the same factors, squared.
  const Column &pv = t.at("xyzobs.px.variance");
  const Column &mv = t.at("xyzobs.mm.variance");
  const double width = Scan::radians(list[0].scan.osc_width);
  check::close(mv.real(0, 0), pv.real(0, 0) * p.pixel_size[0] * p.pixel_size[0],
               1e-15, "fast variance in mm squared");
  check::close(mv.real(0, 2), pv.real(0, 2) * width * width, 1e-18,
               "angular variance in radians squared");
}

TEST(the_millimetre_conversion_is_not_a_division_by_the_pixel_size) {
  // If it were, the parallax correction would be missing from the observation
  // that DIALS refines against, and everything would still look plausible.
  ExperimentList list = one_experiment(0.0);
  list[0].detector.panels[0].parallax = true;
  list[0].detector.panels[0].mu = 3.663;
  list[0].detector.panels[0].thickness = 0.45;
  Table t = observations_of(list, 3);
  add_observed_columns(list, t);
  const Column &px = t.at("xyzobs.px.value");
  const Column &mm = t.at("xyzobs.mm.value");
  const double plain = px.real(2, 0) * list[0].detector[0].pixel_size[0];
  check::is_true(std::abs(mm.real(2, 0) - plain) > 1e-4,
                 "the correction must be in there");
}

TEST(refinement_does_not_move_the_millimetre_centroids) {
  // They are what the spot finder measured through the model as imported. DIALS
  // never recomputes them, which is exactly why they are stale after refinement
  // and must not be used as a join key. Recomputing them here would silently
  // change the observations that refinement was just fitted to.
  const ExperimentList list = one_experiment(0.0);
  Table t = observations_of(list, 5);
  add_observed_columns(list, t);
  std::vector<double> before = t.at("xyzobs.mm.value").reals;

  // Move the detector by a millimetre and recompute what follows the model.
  ExperimentList moved = list;
  const double shift[6] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  moved[0].detector.panels[0] = perturb_panel(list[0].detector[0], shift);
  add_reciprocal_columns(moved, t);

  const std::vector<double> &after = t.at("xyzobs.mm.value").reals;
  for (std::size_t i = 0; i < before.size(); ++i) {
    check::close(after[i], before[i], 0.0, "millimetre centroids are frozen");
  }
  // And s1 did follow it, so the test is not passing because nothing happened.
  check::is_true(t.has("s1"), "s1 was written");
}

TEST(entering_is_a_property_of_the_geometry_not_a_placeholder) {
  // The flag is part of a reflection's identity: the same Miller index can be
  // recorded both entering and exiting in one scan, and a join on the index
  // alone would merge them. A column that is all one value would satisfy any
  // test that only checked it existed.
  const ExperimentList list = one_experiment(0.0);
  Table t = observations_of(list, 12);
  // Spread the observations across the panel so both cases occur.
  Column &xyz = t.real_column("xyzobs.px.value", "vec3<double>", 3);
  for (std::size_t i = 0; i < t.nrows; ++i) {
    xyz.reals[i * 3 + 0] = 100.0 + 350.0 * static_cast<double>(i);
    xyz.reals[i * 3 + 1] = 100.0 + 330.0 * static_cast<double>(i % 6);
    xyz.reals[i * 3 + 2] = 0.5 + static_cast<double>(i % 10);
  }
  add_reciprocal_columns(list, t);
  const Column &entering = t.at("entering");
  std::size_t yes = 0;
  for (std::size_t i = 0; i < t.nrows; ++i) {
    if (entering.integer(i)) ++yes;
  }
  check::is_true(yes > 0 && yes < t.nrows,
                 "both entering and exiting must occur");
}

// --------------------------------------------------------------------------
// the flags column, which is how DIALS asks what a reflection is
// --------------------------------------------------------------------------

TEST(indexing_sets_the_indexed_bit_and_keeps_the_strong_one) {
  // dials.* filters on the flags, not on the Miller indices. A table with
  // correct indices and empty flags processes perfectly and is then invisible
  // to every selection downstream, which is how this went unnoticed.
  Table t;
  t.nrows = 4;
  Column &miller = t.int_column("miller_index", "cctbx::miller::index<>", 3);
  Column &flags = t.int_column("flags", "std::size_t", 1);
  for (std::size_t i = 0; i < t.nrows; ++i) flags.ints[i] = flag::kStrong;
  const int indices[4][3] = {{1, 2, 3}, {0, 0, 0}, {-4, 5, 0}, {0, 0, 0}};
  for (std::size_t i = 0; i < t.nrows; ++i) {
    for (std::size_t k = 0; k < 3; ++k) miller.ints[i * 3 + k] = indices[i][k];
  }

  set_indexed_flags(t);
  const Column &after = t.at("flags");
  // 36 is strong | indexed, which is what a real DIALS indexed.refl carries.
  check::equal(after.integer(0), flag::kStrong | flag::kIndexed, "indexed row");
  check::equal(after.integer(1), flag::kStrong, "unindexed row keeps strong");
  check::equal(after.integer(2), flag::kStrong | flag::kIndexed, "negative indices count");
  check::equal(after.integer(3), flag::kStrong, "and the last one too");
}

TEST(setting_a_flag_does_not_destroy_the_column) {
  // int_column REPLACES a column with a zeroed one, which is right for a
  // derived column recomputed in full and wrong for one that must be read
  // before it is written. Setting the indexed bit that way silently threw away
  // the strong bit dials.find_spots had set, and nothing failed until
  // something downstream filtered on it.
  Table t;
  t.nrows = 2;
  Column &flags = t.int_column("flags", "std::size_t", 1);
  flags.ints[0] = flag::kStrong;
  flags.ints[1] = flag::kStrong | flag::kObserved;

  Column &again = t.modify_int_column("flags", "std::size_t", 1);
  check::equal(again.ints[0], flag::kStrong, "the value survives");
  check::equal(again.ints[1], flag::kStrong | flag::kObserved, "all bits do");

  // And the destructive one really is destructive, so the distinction is real
  // rather than two names for the same thing.
  Column &replaced = t.int_column("flags", "std::size_t", 1);
  check::equal(replaced.ints[0], 0, "int_column clears");
}

TEST(an_index_removed_on_a_later_pass_clears_its_flag) {
  // The bit is cleared as well as set. A reflection indexed on one macrocycle
  // and dropped on the next must stop claiming to be indexed, or the flag and
  // the Miller index disagree and the two ways of asking give different
  // answers.
  Table t;
  t.nrows = 1;
  Column &miller = t.int_column("miller_index", "cctbx::miller::index<>", 3);
  Column &flags = t.int_column("flags", "std::size_t", 1);
  miller.ints[0] = 1;
  flags.ints[0] = flag::kStrong;
  set_indexed_flags(t);
  check::equal(t.at("flags").integer(0), flag::kStrong | flag::kIndexed, "set");

  for (std::size_t k = 0; k < 3; ++k) {
    t.int_column("miller_index", "cctbx::miller::index<>", 3);
  }
  set_indexed_flags(t);
  check::equal(t.at("flags").integer(0), flag::kStrong, "and cleared again");
}

TEST(refinement_records_which_reflections_it_used_and_which_it_threw_out) {
  // The four values a real DIALS indexed.refl carries are 32, 36, 44 and
  // 131108: strong; strong | indexed; strong | indexed | used_in_refinement;
  // and strong | indexed | centroid_outlier. Nothing else appears, and the
  // outlier bit never coincides with the used bit.
  Table table;
  table.nrows = 5;
  Column &flags = table.int_column("flags", "std::size_t", 1);
  for (std::size_t i = 0; i < table.nrows; ++i) flags.ints[i] = flag::kStrong | flag::kIndexed;

  RefineResult result;
  result.rows_used = {0, 2};
  result.rows_rejected = {1, 4};
  set_refinement_flags(result, table);

  const Column &after = table.at("flags");
  check::equal(after.integer(0), flag::kStrong | flag::kIndexed | flag::kUsedInRefinement, "used");
  check::equal(after.integer(1), flag::kStrong | flag::kIndexed | flag::kCentroidOutlier, "rejected");
  check::equal(after.integer(2), flag::kStrong | flag::kIndexed | flag::kUsedInRefinement, "used");
  // Neither used nor rejected: never a candidate, because its rotation angle
  // was not determined. It keeps the indexed bit and gains nothing.
  check::equal(after.integer(3), flag::kStrong | flag::kIndexed, "never a candidate");
  check::equal(after.integer(4), flag::kStrong | flag::kIndexed | flag::kCentroidOutlier, "rejected");

  // And a second pass must not leave the first one's verdict behind.
  RefineResult again;
  again.rows_used = {1};
  again.rows_rejected = {};
  set_refinement_flags(again, table);
  const Column &twice = table.at("flags");
  check::equal(twice.integer(0), flag::kStrong | flag::kIndexed, "used bit cleared");
  check::equal(twice.integer(1), flag::kStrong | flag::kIndexed | flag::kUsedInRefinement,
               "and the outlier bit too, before the new verdict");
}
