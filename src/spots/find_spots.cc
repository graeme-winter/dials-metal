// mxi_find -- the extended dispersion spot finder on an Apple
// GPU, writing a DIALS reflection table.
//
// The threshold in src/dext.{hh,cc} is a transcription of DIALS'
// DispersionExtendedThreshold and the kernels in src/dext_metal.{cc,metal} are
// a transcription of that. What this adds is everything DIALS does after the
// per-frame pixel list: three-dimensional grouping, centroids, and the file
// format -- src/dials_spots.{hh,cc} and src/refl.{hh,cc}.
//
//   dials.import /data/ins10_1_master.h5
//   mxi_find -gpu -e imported.expt
//   dials.index imported.expt strong.refl
//
// The .expt is dials.import's business and is read rather than written: the
// beam, the goniometer and the detector belong to dxtbx, and a second model of
// them here would be a second thing to keep in step. What is taken from it is
// the scan's image range, the panel size and the experiment identifier.
//
// Frames are read and thresholded in parallel and grouped in one thread, since
// the grouping is inherently sequential -- a reflection spans frames and a
// component is only finished when a frame arrives without it. Frames therefore
// have to reach the grouping in order, and they are collected a chunk at a
// time: dispatch a chunk, wait for it, sort it, group it. A chunk is a few
// frames per thread, so the barrier costs a fraction of one frame's latency
// per chunk and the alternative -- resequencing a stream whose frame numbers
// may have gaps in it, since a chunk the writer never received is a frame that
// does not exist -- is a great deal more machinery for that fraction.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "../log_mirror.hh"
#include "decompress.hh"
#include "dext.hh"
#include "dials_spots.hh"
#include "expt.hh"
#include "histogram.hh"
#include "queue.hh"
#include "refl.hh"
#include "series.hh"
#include "signal_pixel.hh"

#ifdef SPOTFINDER_GPU
#include "dext_gpu.hh"
#endif

namespace {

volatile std::sig_atomic_t interrupted = 0;

void on_signal(int /*signal*/) { interrupted = 1; }

using Clock = std::chrono::steady_clock;

struct Options {
  int threads = 4;
  int timeout_seconds = 60;
  int poll_milliseconds = 200;
  std::string master;                 // the NXmx HDF5 master file
  std::string experiments;            // -e: what dials.import wrote
  std::string output = "strong.refl"; // -o
  bool gpu = false;
  bool shoeboxes = true;
  bool two_d = false;
  bool z_offset_given = false;
  std::int64_t z_offset = 0;
  dials_spots::Options grouping;
  std::size_t chunk = 0; // frames collected before grouping; 0 picks one
};

// Which backend, and whether its arithmetic is the CPU's. A reflection table
// from a fast-math build is not one to compare against DIALS without knowing
// that, so it is in the one line that identifies the binary.
void report_version(const char *program) {
  std::printf("%s %s (%s%s)\n", program, SPOTFINDER_VERSION,
#ifdef SPOTFINDER_GPU
              gpu::backend(),
#else
              "no GPU",
#endif
#ifdef SPOTFINDER_FAST_MATH
              ", fast math"
#else
              ""
#endif
  );
}

// To standard output when it was asked for, and to standard error when it is
// the answer to a mistake: `mxi_find --help | less` showed nothing while it
// all went to standard error.
void usage(const char *program, std::FILE *to = stderr) {
  std::fprintf(
      to,
      "usage: %s [-j threads] [-gpu] [-e imported.expt] [-o strong.refl]\n"
      "       [options] [master.nxs]\n"
      "\n"
      "  master.nxs         an NXmx HDF5 master file, or -x master.nxs.\n"
      "                     Optional when -e names an .expt: dials.import\n"
      "                     already recorded the file in its imageset block,\n"
      "                     and -x overrides it if the data has moved.\n"
      "  -e imported.expt   what dials.import wrote, for the scan range, the\n"
      "                     panel size and the experiment identifier\n"
      "  -o file            where to write the reflection table (strong.refl)\n"
      "  -j threads         frames read and thresholded at once (default 4)\n"
      "  -gpu               run the threshold on the GPU; 16-bit only under\n"
      "                     Metal, which has no double precision\n"
      "  --no-shoeboxes     leave out the pixel data, which is most of the\n"
      "                     file and is not needed for indexing\n"
      "  --min-spot-size N  contiguous pixels a spot needs (3)\n"
      "  --max-spot-size N  and the most it may have (1000)\n"
      "  --max-separation D peak to centroid, in pixels; 0 turns it off (2)\n"
      "  --2d               group each frame on its own, as for stills\n"
      "  --z-offset N       array index of image number 0; taken from -e when\n"
      "                     that is given, and 0 otherwise\n"
      "  --chunk N          frames held between groupings (a few per thread)\n"
      "  -t timeout         seconds to wait for the file to appear (default "
      "60)\n"
      "  -p poll-ms         interval between checks while waiting (default "
      "200)\n"
      "  --version          what this binary is, and what it was built with\n",
      program);
}

bool parse_options(int argc, char **argv, Options *options) {
  for (int i = 1; i < argc; i++) {
    const std::string flag = argv[i];
    const bool has_value = i + 1 < argc;
    if (flag == "-j" && has_value) {
      options->threads = std::atoi(argv[++i]);
    } else if (flag == "-t" && has_value) {
      options->timeout_seconds = std::atoi(argv[++i]);
    } else if (flag == "-p" && has_value) {
      options->poll_milliseconds = std::atoi(argv[++i]);
    } else if (flag == "-x" && has_value) {
      options->master = argv[++i];
    } else if (flag == "-e" && has_value) {
      options->experiments = argv[++i];
    } else if (flag == "-o" && has_value) {
      options->output = argv[++i];
    } else if (flag == "-gpu") {
      options->gpu = true;
    } else if (flag == "--no-shoeboxes") {
      options->shoeboxes = false;
    } else if (flag == "--2d") {
      options->two_d = true;
    } else if (flag == "--min-spot-size" && has_value) {
      options->grouping.min_spot_size =
          static_cast<std::size_t>(std::atoll(argv[++i]));
    } else if (flag == "--max-spot-size" && has_value) {
      options->grouping.max_spot_size =
          static_cast<std::size_t>(std::atoll(argv[++i]));
    } else if (flag == "--max-separation" && has_value) {
      options->grouping.max_separation = std::atof(argv[++i]);
    } else if (flag == "--z-offset" && has_value) {
      options->z_offset = std::atoll(argv[++i]);
      options->z_offset_given = true;
    } else if (flag == "--chunk" && has_value) {
      options->chunk = static_cast<std::size_t>(std::atoll(argv[++i]));
    } else if (flag == "--version") {
      report_version("mxi_find");
      std::exit(0);
    } else if (flag == "-h" || flag == "--help") {
      // Asked for, so not an error: it exited 2, like an unknown flag, which
      // made a script checking whether this is installed think it was broken.
      usage(argv[0], stdout);
      std::exit(0);
    } else if (!flag.empty() && flag[0] != '-' && options->master.empty()) {
      // The master file may be named without -x, since it is the only thing
      // this reads and having to flag it would be ceremony.
      options->master = flag;
    } else {
      usage(argv[0]);
      return false;
    }
  }
  if (options->threads < 1) {
    usage(argv[0]);
    return false;
  }
  // The master file may come from the .expt instead of the command line.
  // dials.import already recorded where the images are, in the imageset block,
  // so making the operator repeat it only creates an opportunity for the two to
  // disagree -- and spots found in one file and indexed against the geometry of
  // another is a mistake nothing downstream catches.
  //
  // An explicit -x still wins, because the .expt records an absolute path and a
  // dataset that has moved since import would otherwise be unusable.
  if (options->master.empty() && !options->experiments.empty()) {
    expt::Info info;
    try {
      info = expt::read(options->experiments);
    } catch (const std::exception &error) {
      std::fprintf(stderr, "%s\n", error.what());
      return false;
    }
    if (!info.has_imageset) {
      std::fprintf(stderr,
                   "%s has no imageset, so it does not say where the images "
                   "are; name the master file, or pass -x\n",
                   options->experiments.c_str());
      return false;
    }
    if (info.templated) {
      // Hashes stand for a numbered sequence of files, which this reads none
      // of. Saying so beats handing the path to HDF5 and reporting whatever it
      // makes of it.
      std::fprintf(stderr,
                   "%s names a file template, '%s', which is a numbered "
                   "sequence rather than one NXmx file; this reads NXmx only\n",
                   options->experiments.c_str(), info.image_file.c_str());
      return false;
    }
    options->master = info.image_file;
    std::fprintf(stdout, "Images: %s, from %s\n", options->master.c_str(),
                 options->experiments.c_str());
  }
  if (options->master.empty()) {
    usage(argv[0]);
    return false;
  }
  options->grouping.two_d = options->two_d;
  if (options->chunk == 0)
    options->chunk = static_cast<std::size_t>(options->threads) * 4 + 4;
  return true;
}

void sleep_for(int milliseconds) {
  std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds));
}

// ---------------------------------------------------------------------------
// The frame buffer, and the threshold.
//
// The frame buffer is allocated through gpu::host_alloc when a device is in
// use, which on Metal is a shared MTLBuffer: the frame is decompressed straight
// into memory the kernels read and there is no copy in either direction.
// ---------------------------------------------------------------------------

#ifdef SPOTFINDER_GPU
// Only meaningful in a build that has a device to be on; without one the
// threshold has a single path and this would be an unused variable.
bool on_device = false;
#endif

class FrameBuffer {
public:
  ~FrameBuffer() { release(); }
  FrameBuffer() = default;
  FrameBuffer(const FrameBuffer &) = delete;
  FrameBuffer &operator=(const FrameBuffer &) = delete;

  std::uint8_t *get(std::size_t bytes) {
    if (bytes > capacity_) {
      release();
#ifdef SPOTFINDER_GPU
      // Shared with the device on Metal, so the frame is decompressed straight
      // into memory the kernels read and there is no copy at all.
      if (on_device)
        data_ = static_cast<std::uint8_t *>(gpu::host_alloc(bytes));
      if (data_ != nullptr)
        pinned_ = true;
#endif
      if (data_ == nullptr)
        data_ = new std::uint8_t[bytes];
      capacity_ = bytes;
    }
    return data_;
  }

private:
  void release() {
    if (data_ == nullptr)
      return;
      // Spelled out for both builds. With the GPU compiled out, the dangling
      // 'else' used here left a bare delete indented as though it belonged to
      // the early return above it -- correct, and exactly the shape that gets
      // edited wrongly by the next person to read it.
#ifdef SPOTFINDER_GPU
    if (pinned_) {
      gpu::host_free(data_);
    } else {
      delete[] data_;
    }
#else
    delete[] data_;
#endif
    data_ = nullptr;
    pinned_ = false;
    capacity_ = 0;
  }

  std::uint8_t *data_ = nullptr;
  std::size_t capacity_ = 0;
  bool pinned_ = false;
};

template <typename T>
void threshold(const std::uint8_t *frame, std::vector<SignalPixel> *signal,
               std::size_t height, std::size_t width) {
  const T *const pixels = reinterpret_cast<const T *>(frame);
  int status = 0;
#ifdef SPOTFINDER_GPU
  if (on_device) {
    status = gpu::find<T>(pixels, *signal, height, width);
  } else
#endif
  {
    static thread_local dext_scratch<T> scratch;
    status = dext<T>(pixels, *signal, height, width, scratch);
  }
  if (status == -2) {
    throw std::runtime_error(
        "more signal pixels than the device buffer holds; a frame that dense "
        "is not a frame of spots");
  }
  if (status != 0)
    throw std::runtime_error("the threshold rejected the frame");
}

// ---------------------------------------------------------------------------
// One chunk of frames in flight: the workers fill it and the main thread
// groups it once every frame in it is accounted for.
// ---------------------------------------------------------------------------

struct Found {
  std::int64_t number = 0;
  std::vector<SignalPixel> pixels;
};

class Chunk {
public:
  void expect(std::size_t frames) {
    std::lock_guard<std::mutex> lock(mutex_);
    outstanding_ += frames;
  }

  // Called once per dispatched key, whatever became of it.
  void done(Found *found) {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (found != nullptr)
        frames_.push_back(std::move(*found));
      outstanding_--;
    }
    empty_.notify_all();
  }

  // Set when the last worker has gone. Without it, a chunk dispatched to
  // threads that all failed to open the series would be waited on forever.
  void abandon() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      abandoned_ = true;
    }
    empty_.notify_all();
  }

  bool abandoned() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return abandoned_ && outstanding_ > 0;
  }

  // Frames of this chunk, in ascending order, once they are all in.
  std::vector<Found> collect() {
    std::unique_lock<std::mutex> lock(mutex_);
    empty_.wait(lock, [this] { return outstanding_ == 0 || abandoned_; });
    std::vector<Found> frames = std::move(frames_);
    frames_.clear();
    std::sort(frames.begin(), frames.end(), [](const Found &a, const Found &b) {
      return a.number < b.number;
    });
    return frames;
  }

private:
  mutable std::mutex mutex_;
  std::condition_variable empty_;
  std::size_t outstanding_ = 0;
  bool abandoned_ = false;
  std::vector<Found> frames_;
};

series::Info open_series(series::Series *source, const Options &options) {
  const Clock::time_point deadline =
      Clock::now() + std::chrono::seconds(options.timeout_seconds);
  bool announced = false;
  while (interrupted == 0) {
    series::Info info;
    if (source->try_open(&info))
      return info;
    if (Clock::now() > deadline)
      throw std::runtime_error("timed out waiting for " + source->describe());
    if (!announced) {
      std::fprintf(stdout, "Waiting for %s\n", source->describe().c_str());
      announced = true;
    }
    sleep_for(options.poll_milliseconds);
  }
  throw std::runtime_error("interrupted while waiting for a series");
}

// What the experiment list says, against what the series says. A mismatch here
// is worth stopping for: spot finding would succeed and index nothing, and the
// reason would not be visible in either file.
void reconcile(const expt::Info &experiments, const series::Info &series,
               Options *options) {
  std::fprintf(stdout, "Experiments: %s\n", describe(experiments).c_str());

  if (experiments.experiments != 1) {
    throw std::runtime_error("this writes one experiment's spots, and " +
                             options->experiments + " holds " +
                             std::to_string(experiments.experiments) +
                             "; slice it with dials.split_experiments first");
  }
  if (experiments.panels > 1) {
    throw std::runtime_error(
        "the detector in " + options->experiments + " has " +
        std::to_string(experiments.panels) +
        " panels, and a frame is treated here as one panel; a segmented "
        "detector needs a panel number per spot, which this does not do");
  }
  if (experiments.image_slow > 0 && (experiments.image_slow != series.height ||
                                     experiments.image_fast != series.width)) {
    throw std::runtime_error(
        "the panel in " + options->experiments + " is " +
        std::to_string(experiments.image_fast) + " x " +
        std::to_string(experiments.image_slow) + " and the frames are " +
        std::to_string(series.width) + " x " + std::to_string(series.height) +
        " (fast x slow); these are not the same images");
  }
  // A scan covering fewer images than the file is not a mismatch to warn
  // about: the frames outside it are simply not read. It used to warn that z
  // would be wrong outside the scan and then read them all anyway, which was
  // the warning describing the bug rather than preventing it. A scan claiming
  // MORE images than the file has is still an error, since there is nothing
  // to read.
  if (experiments.has_scan && series.images > 0 &&
      static_cast<std::uint64_t>(experiments.images()) > series.images) {
    throw std::runtime_error(
        options->experiments + " covers " +
        std::to_string(experiments.images()) + " images and the file has " +
        std::to_string(series.images) + "; the scan reaches past the data");
  }

  // The imageset and the scan can disagree about how many images there are --
  // one sliced and the other not, or an .expt assembled by hand. The scan sets
  // z, so a mismatch means z is measured against a range the images do not
  // cover.
  if (!experiments.imageset_matches_scan()) {
    std::fprintf(stderr,
                 "warning: %s lists %llu frames in its imageset and %lld in "
                 "its scan; z follows the scan\n",
                 options->experiments.c_str(),
                 static_cast<unsigned long long>(experiments.frames),
                 static_cast<long long>(experiments.images()));
  }
  // Gaps in single_file_indices mean the imageset is not the whole file in
  // order, which nothing here allows for: frames are read as a contiguous run,
  // so a gap would silently shift every spot after it onto the wrong image.
  if (experiments.has_imageset && !experiments.contiguous_indices) {
    throw std::runtime_error(
        "the imageset in " + options->experiments +
        " skips frames, and this reads a contiguous run; z would be wrong "
        "from the first gap onwards");
  }

  // z is measured in scan array indices, where image n is n - 1. A frame
  // arrives numbered by its index IN THE FILE, and the .expt says which file
  // index each scan image is: single_file_indices. So the offset is the
  // difference between the two, which is nothing when they line up.
  //
  // It was first_image - 1 unconditionally, which assumed a sliced import's
  // frames are numbered from zero -- true of a file holding only the slice,
  // false of a master file for the whole run, where image 6 already arrives
  // as frame 5. Reading only the scan's frames made that visible: a scan of
  // images 6 to 10 put its spots at z of 10.5 to 14.5 instead of 5.5 to 9.5.
  //
  // Without single_file_indices the file index of image n is taken to be
  // n - 1, the same assumption the frame restriction makes, so that the frames
  // read and the z they are given cannot disagree.
  if (!options->z_offset_given && experiments.has_scan) {
    const std::int64_t first_index = experiments.frames > 0
                                         ? experiments.first_index
                                         : experiments.first_image - 1;
    options->z_offset = (experiments.first_image - 1) - first_index;
    if (options->z_offset != 0) {
      std::fprintf(stdout,
                   "The scan starts at image %lld, so z starts at %lld\n",
                   static_cast<long long>(experiments.first_image),
                   static_cast<long long>(options->z_offset));
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  // Mirrored to mxi_find.log in the working directory, as DIALS writes
  // dials.find_spots.log; not for a run that only asks for help.
  if (!mxi::only_asks_for_help(argc, argv))
    mxi::mirror_to_log("mxi_find.log");
  Options options;
  if (!parse_options(argc, argv, &options))
    return 2;

  std::signal(SIGINT, on_signal);
  std::signal(SIGTERM, on_signal);

  std::unique_ptr<series::Series> source;
  series::Info info;
  // How many frames will be offered, when fewer than the file holds.
  std::uint64_t restricted = 0;
  expt::Info experiments;
  try {
    if (!options.experiments.empty())
      experiments = expt::read(options.experiments);

    source = series::nxmx(options.master);

    // Before the series is opened, so that asking for a device that is not
    // there fails at once rather than after a wait.
    if (options.gpu) {
#ifdef SPOTFINDER_GPU
      if (!gpu::available()) {
        throw std::runtime_error(std::string("-gpu was given but no usable ") +
                                 gpu::backend() + " device was found");
      }
      on_device = true;
#else
      throw std::runtime_error(
          "-gpu was given but this build has no GPU support; configure with "
          "-DSPOTFINDER_METAL=ON or -DSPOTFINDER_CUDA=ON");
#endif
    }

    info = open_series(source.get(), options);
    if (!options.experiments.empty())
      reconcile(experiments, info, &options);

    // Only the frames the scan covers. The array indices come from the
    // imageset's single_file_indices when it has them, which is what they are;
    // otherwise from image_range, whose image n is array index n - 1.
    //
    // Found by spot finding a scan of 1800 images of a 36000 image file, which
    // warned that z would be wrong outside the scan and then processed all
    // 36000 -- twenty times the work, and 1115233 spots of which a
    // twentieth belonged to the experiment it was asked about.
    if (experiments.has_scan) {
      std::uint64_t first = 0, last = 0;
      if (experiments.frames > 0) {
        first = static_cast<std::uint64_t>(experiments.first_index);
        last = static_cast<std::uint64_t>(experiments.last_index);
      } else {
        first = static_cast<std::uint64_t>(
            std::max<std::int64_t>(experiments.first_image - 1, 0));
        last = static_cast<std::uint64_t>(
            std::max<std::int64_t>(experiments.last_image - 1, 0));
      }
      if (!source->restrict_frames(first, last)) {
        throw std::runtime_error(
            "this series cannot be limited to the frames the .expt's scan "
            "covers, and finding spots on the rest would put them at z "
            "values the experiment does not have");
      }
      // What "done" means now: the scan's frames, not the file's. Left at the
      // file's count, the loop would wait for frames it has been told not to
      // offer until the timeout ended it.
      restricted = last - first + 1;
    }
  } catch (const std::exception &error) {
    std::fprintf(stderr, "%s\n", error.what());
    return 1;
  }

  if (restricted > 0 && restricted != info.images) {
    std::fprintf(stdout,
                 "Reading %llu of the %llu images, the ones the scan covers\n",
                 static_cast<unsigned long long>(restricted),
                 static_cast<unsigned long long>(info.images));
  }
  std::fprintf(
      stdout,
      "Series %s: %llu images of %llu x %llu pixels (fast x slow), from %s, "
      "%d thread%s, "
      "on the %s\n",
      info.name.c_str(), static_cast<unsigned long long>(info.images),
      static_cast<unsigned long long>(info.width),
      static_cast<unsigned long long>(info.height), source->describe().c_str(),
      options.threads, options.threads == 1 ? "" : "s",
      options.gpu ? "GPU" : "CPU");
  std::fprintf(
      stdout,
      "Grouping in %s, %zu to %zu pixels a spot, peak within %.1f of the "
      "centroid\n",
      options.two_d ? "two dimensions" : "three dimensions",
      options.grouping.min_spot_size, options.grouping.max_spot_size,
      options.grouping.max_separation);

  dials_spots::Labeller labeller(static_cast<std::size_t>(info.height),
                                 static_cast<std::size_t>(info.width),
                                 options.grouping);

  work::Queue<std::string> queue(static_cast<std::size_t>(options.threads) * 4);
  Chunk chunk;
  std::atomic<std::uint64_t> read{0};
  std::atomic<std::uint64_t> missing{0};
  std::atomic<std::uint64_t> failed{0};
  std::atomic<int> live{options.threads};
  const Clock::time_point began = Clock::now();

  std::vector<std::thread> workers;
  for (int i = 0; i < options.threads; i++) {
    workers.emplace_back([&, i] {
      // However this thread leaves, the last one out says so: a chunk in
      // flight is otherwise waited on forever.
      struct Leaving {
        std::atomic<int> &live;
        Chunk &chunk;
        ~Leaving() {
          if (live.fetch_sub(1) == 1)
            chunk.abandon();
        }
      } leaving{live, chunk};

      std::unique_ptr<series::Reader> reader;
      try {
        reader = source->reader();
      } catch (const std::exception &error) {
        // The keys this thread would have taken are still in the queue for the
        // others; only if every thread fails does the run stop.
        std::fprintf(stderr, "thread %d: %s\n", i, error.what());
        return;
      }

      FrameBuffer buffer;
      series::Frame frame;
      std::string key;
      while (queue.pop(&key)) {
        Found found;
        bool have = false;
        try {
          if (!reader->read(key, &frame)) {
            missing++;
          } else {
            const std::size_t height = static_cast<std::size_t>(frame.height);
            const std::size_t width = static_cast<std::size_t>(frame.width);
            const std::size_t bytes =
                decompress::frame_bytes(height, width, frame.bit_depth);
            std::uint8_t *const pixels = buffer.get(bytes);
            decompress::image(frame.data, frame.algorithm, frame.bit_depth,
                              height, width, {pixels, bytes});
            found.number = frame.number;
            switch (frame.bit_depth) {
            case 16:
              threshold<std::uint16_t>(pixels, &found.pixels, height, width);
              break;
            case 32:
              threshold<std::uint32_t>(pixels, &found.pixels, height, width);
              break;
            default:
              throw std::runtime_error("the threshold does not support " +
                                       std::to_string(frame.bit_depth) +
                                       "-bit data");
            }
            have = true;
            read++;
          }
        } catch (const std::exception &error) {
          failed++;
          have = false;
          std::fprintf(stderr, "%s: %s\n", key.c_str(), error.what());
        }
        chunk.done(have ? &found : nullptr);
      }
    });
  }

  // Dispatch in chunks, group each chunk once it is complete. Keys are zero
  // padded decimal, so sorting them as text sorts them as numbers.
  std::set<std::string> dispatched;
  std::vector<std::string> waiting;
  Clock::time_point progressed = Clock::now();
  std::int64_t highest = 0;
  bool any = false;
  std::uint64_t grouped = 0;
  std::uint64_t out_of_order = 0;
  bool stop = false;

  const auto run = [&](std::size_t frames) {
    chunk.expect(frames);
    for (std::size_t i = 0; i < frames; i++) {
      if (!queue.push(waiting[i])) {
        // The queue only closes at the end, but if it ever did, the chunk
        // still has to be balanced.
        chunk.done(nullptr);
      }
    }
    waiting.erase(waiting.begin(),
                  waiting.begin() + static_cast<std::ptrdiff_t>(frames));

    const std::vector<Found> collected = chunk.collect();
    if (chunk.abandoned()) {
      std::fprintf(stderr,
                   "every thread has gone, so the rest of the series cannot be "
                   "read; stopping with what was grouped so far\n");
      failed++;
      stop = true;
    }
    for (const Found &found : collected) {
      const std::int64_t z = found.number + options.z_offset;
      if (any && z <= highest) {
        // Only reachable against a live source that published an image out of
        // order. Dropping it is better than throwing away the whole run, and
        // saying so is better than a quiet gap.
        out_of_order++;
        continue;
      }
      try {
        labeller.add(z, found.pixels);
      } catch (const std::exception &error) {
        std::fprintf(stderr, "grouping frame %lld: %s\n",
                     static_cast<long long>(found.number), error.what());
        failed++;
        continue;
      }
      highest = z;
      any = true;
      grouped++;
    }
  };

  while (!stop) {
    std::size_t added = 0;
    try {
      for (const std::string &key : source->ready()) {
        if (!dispatched.insert(key).second)
          continue;
        waiting.push_back(key);
        added++;
      }
    } catch (const std::exception &error) {
      std::fprintf(stderr, "%s\n", error.what());
    }
    if (added > 0) {
      progressed = Clock::now();
      std::sort(waiting.begin(), waiting.end());
    }

    const bool complete = (restricted > 0 && dispatched.size() >= restricted) ||
                          (restricted == 0 && info.images > 0 &&
                           dispatched.size() >= info.images) ||
                          (source->finished() && added == 0);
    const bool timed_out = Clock::now() - progressed >
                           std::chrono::seconds(options.timeout_seconds);
    if (interrupted != 0 || complete || timed_out)
      stop = true;
    if (timed_out && !complete) {
      std::fprintf(stderr, "No new images for %d seconds, giving up\n",
                   options.timeout_seconds);
    }

    // Whole chunks while more may arrive; whatever is left once nothing can.
    while (waiting.size() >= options.chunk)
      run(options.chunk);
    if (stop && !waiting.empty())
      run(waiting.size());

    if (!stop)
      sleep_for(options.poll_milliseconds);
  }

  queue.close();
  for (std::thread &worker : workers)
    worker.join();

  labeller.finish();

  const double seconds =
      std::chrono::duration<double>(Clock::now() - began).count();
  const dials_spots::Counts &counts = labeller.counts();

  std::fprintf(stdout,
               "Thresholded %llu of %llu images in %.1f s (%.1f images/s), "
               "%llu never written, %llu failures\n",
               static_cast<unsigned long long>(read.load()),
               static_cast<unsigned long long>(dispatched.size()), seconds,
               seconds > 0 ? read.load() / seconds : 0.0,
               static_cast<unsigned long long>(missing.load()),
               static_cast<unsigned long long>(failed.load()));
  if (out_of_order > 0) {
    std::fprintf(stderr,
                 "warning: %llu frames arrived after a later one and were "
                 "dropped\n",
                 static_cast<unsigned long long>(out_of_order));
  }
  // dials.find_spots' own summary, wording included, so that a log from this
  // and a log from that can be read the same way and compared line for line.
  // The two "Calculated" lines are one pass here rather than two, so their
  // counts are equal by construction; they are both printed anyway, because a
  // block that reads differently from the familiar one is the thing this
  // replaced.
  const std::uint64_t sized =
      counts.groups - counts.too_small - counts.too_large;

  std::fprintf(stdout, "Found %llu signal pixels on %llu frames\n",
               static_cast<unsigned long long>(counts.signal_pixels),
               static_cast<unsigned long long>(grouped));
  std::fprintf(stdout, "Extracted %llu spots\n",
               static_cast<unsigned long long>(counts.groups));
  std::fprintf(stdout, "Removed %llu spots with size < %zu pixels\n",
               static_cast<unsigned long long>(counts.too_small),
               options.grouping.min_spot_size);
  std::fprintf(stdout, "Removed %llu spots with size > %zu pixels\n",
               static_cast<unsigned long long>(counts.too_large),
               options.grouping.max_spot_size);
  std::fprintf(stdout, "Calculated %llu spot centroids\n",
               static_cast<unsigned long long>(sized));
  std::fprintf(stdout, "Calculated %llu spot intensities\n",
               static_cast<unsigned long long>(sized));
  // Only when the filter ran: "Filtered 48 of 48" would otherwise read as a
  // filter that passed everything rather than one that was switched off.
  if (options.grouping.max_separation > 0.0) {
    std::fprintf(stdout,
                 "Filtered %llu of %llu spots by peak-centroid distance\n",
                 static_cast<unsigned long long>(counts.accepted),
                 static_cast<unsigned long long>(sized));
  }

  // Spots per image, as dials.find_spots draws it: the part of this report a
  // user reads at a glance. A spot's image is floor(z) + 1, image n starting at
  // z = n - 1, over the images the scan covers or, without a scan, the series.
  if (!labeller.spots().empty()) {
    long long first = 1;
    long long last = static_cast<long long>(info.images);
    if (experiments.has_scan) {
      first = experiments.first_image;
      last = experiments.last_image;
    }
    if (last >= first) {
      std::vector<std::size_t> per_image(
          static_cast<std::size_t>(last - first + 1), 0);
      for (const auto &spot : labeller.spots()) {
        const long long image =
            static_cast<long long>(std::floor(spot.position[2])) + 1;
        if (image >= first && image <= last)
          ++per_image[static_cast<std::size_t>(image - first)];
      }
      std::fprintf(stdout, "\nHistogram of spots per image:\n");
      for (const std::string &line : spots::spot_histogram(per_image, first))
        std::fprintf(stdout, "%s\n", line.c_str());
      std::fprintf(stdout, "\n");
    }
  }

  refl::Options writing;
  writing.identifier = experiments.identifier;
  writing.shoeboxes = options.shoeboxes;
  try {
    refl::write(options.output, labeller.spots(), labeller.pixels(),
                static_cast<std::size_t>(info.width), writing);
  } catch (const std::exception &error) {
    std::fprintf(stderr, "%s\n", error.what());
    return 1;
  }

  if (labeller.spots().empty()) {
    std::fprintf(stderr,
                 "No reflections found; %s is a well formed table with no rows "
                 "in it, which dials.index will refuse\n",
                 options.output.c_str());
  } else {
    std::error_code ignored;
    const std::uintmax_t bytes =
        std::filesystem::file_size(options.output, ignored);
    std::fprintf(stdout, "Wrote %zu reflections to %s (%.1f MB%s)\n",
                 labeller.spots().size(), options.output.c_str(),
                 static_cast<double>(bytes) / 1e6,
                 options.shoeboxes ? ", most of it shoeboxes" : "");
  }

  return failed.load() == 0 ? 0 : 1;
}
