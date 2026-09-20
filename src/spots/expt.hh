// Reading the few facts a spot list needs out of a dxtbx experiment list.
//
// An .expt is JSON, and dials.index wants the .expt and the .refl to agree
// about three things: which images the scan covers, so that a spot's z is the
// array index DIALS expects rather than the detector's own numbering; how big a
// panel is, so that a mismatched frame is caught here rather than as a
// nonsensical lattice later; and the experiment identifier, which is the string
// DIALS uses to tie a table to an experiment.
//
// Nothing else is taken from it, and nothing is written back. This deliberately
// does not try to be dxtbx: the beam, the goniometer and the detector's
// hierarchy are dials.import's business, and reimplementing them here would be
// a second model to keep in step with the first.
//
// The JSON reader is a few hundred lines because pulling three numbers out of a
// document with a regular expression is the kind of thing that works until an
// experiment list gains a field.

#ifndef SPOTFINDER_EXPT_HH
#define SPOTFINDER_EXPT_HH

#include <cstddef>
#include <cstdint>
#include <string>

namespace expt {

struct Info {
  std::size_t experiments = 0;

  // The identifier of the first experiment. Empty if dials.import did not set
  // one, which is not an error.
  std::string identifier;

  // The scan, as image numbers, inclusive and counting from one, exactly as
  // dxtbx stores image_range. The array index of image n is n - 1, and that
  // index is what a reflection's z is measured in.
  bool has_scan = false;
  std::int64_t first_image = 0;
  std::int64_t last_image = 0;

  // Where the images are. dials.import records this in the imageset block:
  // `template` names the file and `single_file_indices` lists the array
  // indices within it, so an .expt on its own says everything needed to find
  // the data. Without this the tool has to be told the master file separately
  // and the two can disagree, which is a mistake nothing downstream catches --
  // spots from one file indexed against the geometry of another.
  bool has_imageset = false;
  std::string image_file;
  std::size_t imagesets = 0;

  // A template with '#' placeholders names a numbered sequence of files, not
  // one file. This tool reads NXmx and so cannot follow it; reported rather
  // than resolved, because guessing which file the placeholders stand for is
  // dials.import's job and it has already done it.
  bool templated = false;

  // The entries of single_file_indices. They are zero-based array indices,
  // where image_range counts from one.
  std::size_t frames = 0;
  std::int64_t first_index = 0;
  std::int64_t last_index = 0;
  //: False when the indices skip, which a sliced or filtered import can
  //: produce. Nothing here handles that yet, so it has to be visible.
  bool contiguous_indices = true;

  // Whether the imageset and the scan describe the same number of images. They
  // can disagree -- an .expt assembled by hand, or a slice taken of one but not
  // the other -- and the difference would otherwise appear as spots on the
  // wrong frames.
  bool imageset_matches_scan() const {
    return !has_scan || !has_imageset || frames == 0 ||
           static_cast<std::int64_t>(frames) == images();
  }

  // The first detector's panels. A segmented detector is reported rather than
  // refused here; the caller decides.
  bool has_detector = false;
  std::size_t panels = 0;
  std::size_t image_fast = 0; // pixels across a row
  std::size_t image_slow = 0; // rows

  std::int64_t images() const {
    return has_scan ? last_image - first_image + 1 : 0;
  }
};

// Throws std::runtime_error if the file cannot be read or is not an experiment
// list. A missing scan or detector is not an error: a still has no scan.
Info read(const std::string &path);

// One line for the startup banner.
std::string describe(const Info &info);

} // namespace expt

#endif // SPOTFINDER_EXPT_HH
