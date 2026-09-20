// Reading an experiment list, and the shapes of one that matter.
//
// The documents below are cut down from what dials.import writes -- the fields
// that are read, plus enough of the ones that are not to check that they are
// stepped over rather than tripped on. The awkward cases are the ones that
// arrive from real use: a still with no scan at all, a sliced import whose
// image_range does not start at one, a segmented detector, and a file that is
// not an experiment list because someone passed the reflections by mistake.

#include <cstdio>
#include <cstring>
#include <string>

#include "expt.hh"

namespace {

int failures = 0;

void check(const std::string &what, unsigned long long actual,
           unsigned long long expected) {
  if (actual == expected)
    return;
  std::printf("FAIL %s: %llu, expected %llu\n", what.c_str(), actual, expected);
  failures++;
}

void check_text(const std::string &what, const std::string &actual,
                const std::string &expected) {
  if (actual == expected)
    return;
  std::printf("FAIL %s: \"%s\", expected \"%s\"\n", what.c_str(),
              actual.c_str(), expected.c_str());
  failures++;
}

const char *kPath = "test_expt.expt";

void put(const std::string &text) {
  std::FILE *file = std::fopen(kPath, "wb");
  if (file == nullptr) {
    std::printf("FAIL cannot write %s\n", kPath);
    failures++;
    return;
  }
  std::fwrite(text.data(), 1, text.size(), file);
  std::fclose(file);
}

// A sweep, as dials.import leaves one, with the fields this reads and a
// scattering of the ones it does not.
const char *kSweep = R"({
  "__id__": "ExperimentList",
  "experiment": [
    {
      "__id__": "Experiment",
      "identifier": "f2d0a1c4-6b3e-4a11-9e77-0c5f1b2d3e4f",
      "beam": 0,
      "detector": 0,
      "goniometer": 0,
      "scan": 0,
      "imageset": 0,
      "crystal": null,
      "profile": null
    }
  ],
  "imageset": [
    {
      "__id__": "ImageSequence",
      "template": "/data/ins10_1_master.h5",
      "mask": "",
      "gain": "",
      "params": { "dynamic_shadowing": "Auto", "multi_panel": false }
    }
  ],
  "beam": [
    {
      "__id__": "MonochromaticBeam",
      "direction": [0.0, 0.0, 1.0],
      "wavelength": 0.9795,
      "divergence": 0.0,
      "polarization_normal": [0.0, 1.0, 0.0],
      "polarization_fraction": 0.999,
      "flux": 0.0,
      "transmission": 1.0
    }
  ],
  "detector": [
    {
      "panels": [
        {
          "name": "Panel",
          "type": "SENSOR_PAD",
          "fast_axis": [1.0, 0.0, 0.0],
          "slow_axis": [0.0, -1.0, 0.0],
          "origin": [-155.55, 166.4, -265.27],
          "raw_image_offset": [0, 0],
          "image_size": [4148, 4362],
          "pixel_size": [0.075, 0.075],
          "trusted_range": [0.0, 65535.0],
          "thickness": 0.45,
          "material": "Si",
          "mu": 3.9206,
          "identifier": "",
          "mask": [[1029, 0, 1041, 4362], [2069, 0, 2081, 4362]],
          "gain": 1.0,
          "pedestal": 0.0,
          "px_mm_strategy": { "type": "ParallaxCorrectedPxMmStrategy" }
        }
      ],
      "hierarchy": { "__id__": "ExtendedDetector" }
    }
  ],
  "goniometer": [
    {
      "__id__": "Goniometer",
      "rotation_axis": [1.0, 0.0, 0.0],
      "fixed_rotation": [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
      "setting_rotation": [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    }
  ],
  "scan": [
    {
      "__id__": "Scan",
      "image_range": [1, 600],
      "batch_offset": 0,
      "oscillation": [0.0, 0.1],
      "exposure_time": [0.0665, 0.0665],
      "epochs": [0.0, 0.0665],
      "valid_image_ranges": {}
    }
  ],
  "crystal": [],
  "profile": [],
  "scaling_model": []
})";

// The imageset block, which says where the images actually are. Without it the
// tool has to be told the master file separately, and the two can disagree --
// spots from one file indexed against the geometry of another, which nothing
// downstream catches.
const char *kWithImageset = R"({
  "__id__": "ExperimentList",
  "experiment": [{"__id__": "Experiment", "imageset": 0, "scan": 0,
                  "identifier": "abc"}],
  "imageset": [{"__id__": "ImageSequence",
                "template": "/data/ins10_1.nxs",
                "single_file_indices": [0, 1, 2, 3, 4],
                "dx": null, "dy": null,
                "params": {"dynamic_shadowing": "Auto", "multi_panel": false}}],
  "scan": [{"image_range": [1, 5]}]
})";

// A numbered sequence of CBFs, which this tool cannot read. What matters is
// that it says so rather than handing a path with hashes in it to HDF5.
const char *kTemplatedImageset = R"({
  "__id__": "ExperimentList",
  "experiment": [{"__id__": "Experiment", "imageset": 0}],
  "imageset": [{"__id__": "ImageSequence",
                "template": "/data/lys_1_####.cbf"}]
})";

// Indices with a gap. A sliced or filtered import can produce this, and
// reading first..last as a range would quietly read the frames left out.
const char *kGappyImageset = R"({
  "__id__": "ExperimentList",
  "experiment": [{"__id__": "Experiment", "imageset": 0}],
  "imageset": [{"__id__": "ImageSequence", "template": "/data/x.nxs",
                "single_file_indices": [0, 1, 4, 5]}]
})";

// No imageset at all, which is not an error: an .expt can be assembled without
// one, and everything that worked before this field existed must still work.
const char *kNoImageset = R"({
  "__id__": "ExperimentList",
  "experiment": [{"__id__": "Experiment", "scan": 0}],
  "scan": [{"image_range": [1, 5]}]
})";

// Ten frames in the imageset, five in the scan: one has been sliced and the
// other has not.
const char *kMismatchedImageset = R"({
  "__id__": "ExperimentList",
  "experiment": [{"__id__": "Experiment", "imageset": 0, "scan": 0}],
  "imageset": [{"__id__": "ImageSequence", "template": "/data/x.nxs",
                "single_file_indices": [0,1,2,3,4,5,6,7,8,9]}],
  "scan": [{"image_range": [1, 5]}]
})";

} // namespace

int main() {
  {
    put(kSweep);
    const expt::Info info = expt::read(kPath);
    check("one experiment", info.experiments, 1);
    check_text("the identifier", info.identifier,
               "f2d0a1c4-6b3e-4a11-9e77-0c5f1b2d3e4f");
    check("it has a scan", info.has_scan ? 1 : 0, 1);
    check("the first image", static_cast<unsigned long long>(info.first_image),
          1);
    check("the last image", static_cast<unsigned long long>(info.last_image),
          600);
    check("so 600 images", static_cast<unsigned long long>(info.images()), 600);
    check("it has a detector", info.has_detector ? 1 : 0, 1);
    check("of one panel", info.panels, 1);
    // dxtbx writes image_size fast then slow, which is the other way round
    // from every frame dimension here. Getting this the wrong way round is the
    // whole reason it is checked.
    check("4148 pixels across a row", info.image_fast, 4148);
    check("and 4362 rows", info.image_slow, 4362);
  }

  // A slice, as dials.import ... image_range=501,600 or a sliced experiment
  // gives. The array index of image 501 is 500, which is what a reflection's z
  // has to be measured in.
  {
    std::string text = kSweep;
    const std::size_t at = text.find("\"image_range\": [1, 600]");
    if (at == std::string::npos) {
      std::printf("FAIL the test document has changed shape\n");
      failures++;
    } else {
      text.replace(at, std::strlen("\"image_range\": [1, 600]"),
                   "\"image_range\": [501, 600]");
      put(text);
      const expt::Info info = expt::read(kPath);
      check("the first image of a slice",
            static_cast<unsigned long long>(info.first_image), 501);
      check("and its length", static_cast<unsigned long long>(info.images()),
            100);
    }
  }

  // A still: no scan, which is not an error. dxtbx writes the model index as
  // null rather than leaving the field out.
  {
    std::string text = kSweep;
    const std::size_t at = text.find("\"scan\": 0");
    text.replace(at, std::strlen("\"scan\": 0"), "\"scan\": null");
    put(text);
    const expt::Info info = expt::read(kPath);
    check("a still has no scan", info.has_scan ? 1 : 0, 0);
    check("and no images to speak of",
          static_cast<unsigned long long>(info.images()), 0);
    check("but still a detector", info.panels, 1);
  }

  // Two panels, which this reports and the caller refuses: the finder treats a
  // frame as one panel, and a segmented detector would need a panel number per
  // spot.
  {
    std::string text = kSweep;
    const std::size_t panels = text.find("\"panels\": [");
    const std::size_t open = text.find('{', panels);
    const std::size_t close = text.find("\n        }", open);
    const std::string panel =
        text.substr(open, close + std::strlen("\n        }") - open);
    text.insert(close + std::strlen("\n        }"), ", " + panel);
    put(text);
    const expt::Info info = expt::read(kPath);
    check("two panels are reported", info.panels, 2);
    check("with the first one's size", info.image_fast, 4148);
  }

  // No experiments at all, which is a valid empty list.
  {
    put("{\"__id__\": \"ExperimentList\", \"experiment\": []}");
    const expt::Info info = expt::read(kPath);
    check("an empty list has no experiments", info.experiments, 0);
  }

  // Not an experiment list. A .refl handed in by mistake is msgpack, so it is
  // not even JSON; a .json of something else parses and has no experiments in
  // it. Both have to say so rather than come back empty.
  {
    put("{\"nrows\": 12, \"data\": {}}");
    bool threw = false;
    try {
      expt::read(kPath);
    } catch (const std::exception &) {
      threw = true;
    }
    check("a document with no experiment list throws", threw ? 1 : 0, 1);
  }
  {
    put(std::string("\x93\xbb") + "dials::af::reflection_table");
    bool threw = false;
    try {
      expt::read(kPath);
    } catch (const std::exception &) {
      threw = true;
    }
    check("a msgpack file throws", threw ? 1 : 0, 1);
  }
  {
    put(R"({"experiment": [{"scan": 0}], "scan": [{"image_range": [1, 2)");
    bool threw = false;
    try {
      expt::read(kPath);
    } catch (const std::exception &) {
      threw = true;
    }
    check("a truncated document throws", threw ? 1 : 0, 1);
  }
  {
    bool threw = false;
    try {
      expt::read("no_such_file_at_all.expt");
    } catch (const std::exception &) {
      threw = true;
    }
    check("a missing file throws", threw ? 1 : 0, 1);
  }

  // An identifier with an escape in it, since the reader has to decode one.
  {
    put(R"({"experiment": [{"identifier": "a\/b\u00e9c", "scan": null}]})");
    const expt::Info info = expt::read(kPath);
    check_text("escapes are decoded", info.identifier,
               "a/b\xc3\xa9"
               "c");
  }

  {
    put(kWithImageset);
    const expt::Info info = expt::read(kPath);
    check("it has an imageset", info.has_imageset ? 1 : 0, 1);
    check_text("the image file", info.image_file, "/data/ins10_1.nxs");
    check("one imageset", info.imagesets, 1);
    check("five frames", info.frames, 5);
    check("not templated", info.templated ? 1 : 0, 0);
    check("the indices are contiguous", info.contiguous_indices ? 1 : 0, 1);
    // The indices count from zero and image_range counts from one. Conflating
    // those is the off-by-one this field exists to make visible rather than
    // to hide.
    check("the first image is one",
          static_cast<unsigned long long>(info.first_image), 1);
    check("the first index is zero",
          static_cast<unsigned long long>(info.first_index), 0);
    check("the last index is four",
          static_cast<unsigned long long>(info.last_index), 4);
    check("imageset and scan agree", info.imageset_matches_scan() ? 1 : 0, 1);
  }

  {
    put(kTemplatedImageset);
    const expt::Info info = expt::read(kPath);
    check("templated is seen", info.templated ? 1 : 0, 1);
    // Kept verbatim. Working out what the hashes stand for is dials.import's
    // job and it has already done it; a second guess here could differ.
    check_text("the path is kept as written", info.image_file,
               "/data/lys_1_####.cbf");
    check("no frames listed", info.frames, 0);
  }

  {
    put(kGappyImageset);
    const expt::Info info = expt::read(kPath);
    check("four frames", info.frames, 4);
    check("first index", static_cast<unsigned long long>(info.first_index), 0);
    check("last index", static_cast<unsigned long long>(info.last_index), 5);
    check("the gap is noticed", info.contiguous_indices ? 1 : 0, 0);
  }

  {
    put(kMismatchedImageset);
    const expt::Info info = expt::read(kPath);
    check("ten frames", info.frames, 10);
    check("five images in the scan",
          static_cast<unsigned long long>(info.images()), 5);
    check("the disagreement is noticed",
          info.imageset_matches_scan() ? 1 : 0, 0);
  }

  {
    // The fixture that was already here carries an imageset, so it doubles as a
    // check that this is read out of a document written the way dials.import
    // writes one rather than only out of the cut-down examples above.
    put(kSweep);
    const expt::Info info = expt::read(kPath);
    check("the existing fixture has an imageset", info.has_imageset ? 1 : 0, 1);
    check_text("and names its file", info.image_file,
               "/data/ins10_1_master.h5");
  }

  {
    // An .expt without one is not an error, and everything that worked before
    // this field existed must still work.
    put(kNoImageset);
    const expt::Info info = expt::read(kPath);
    check("no imageset here", info.has_imageset ? 1 : 0, 0);
    check("no frames", info.frames, 0);
    check("contiguous by default", info.contiguous_indices ? 1 : 0, 1);
    check("matches vacuously", info.imageset_matches_scan() ? 1 : 0, 1);
    check("the scan is still read",
          static_cast<unsigned long long>(info.images()), 5);
  }


  std::remove(kPath);

  std::printf("%s: the experiment list reader, %d failures\n",
              failures == 0 ? "PASS" : "FAILED", failures);
  return failures == 0 ? 0 : 1;
}
