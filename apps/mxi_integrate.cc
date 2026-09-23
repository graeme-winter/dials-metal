// mxi_integrate: summation integration on real images.
//
//   mxi_integrate integrated.expt master.nxs -o mine.refl --d-min 2.0
//   mxi_integrate ... --save-shoeboxes            # keep the pixels and the mask
//
// Predicts, builds each reflection's measurement box, fills it from the
// images, fits the background and sums the foreground. The output is a DIALS
// reflection table, so `dials.show`, `dials.image_viewer` and the scaler all
// read it, and `intensity.sum.value` can be compared with DIALS' own directly.
//
// --save-shoeboxes keeps the pixels and the mask in the table, as mxi_mask
// does. They are large -- a megabyte per hundred reflections -- and they are
// the difference between "the intensity is wrong" and "the intensity is wrong
// because the mask is here and the spot is there".

#include <algorithm>
#include <array>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "args.h"
#include "background.h"
#include "expt.h"
#include "integrate.h"
#include "mask.h"
#include "predict.h"
#include "profile_model.h"
#include "refl.h"
#include "reference.h"
#include "shoebox.h"

#include "decompress.hh"
#include "series.hh"

#include <memory>
#include <stdexcept>

using namespace mxi;

namespace {

//: Wall clock, in seconds. Wall rather than CPU: the question is how long
//: someone waits, and on a seven minute job most of the answer may be the
//: disk.
double now_wall() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}


//: Frames decompressed ahead of the loop that consumes them.
//:
//: Decompressing is half of a large integration -- 158 of 310 seconds on ten
//: rotations -- and it is per frame with nothing shared, so it threads. What
//: does NOT thread is the bookkeeping: a shoebox spans several frames and is
//: opened, filled and closed in frame order, so the consumer must see frames
//: in order and one at a time.
//:
//: So the threads only fetch and decompress. They put finished frames into a
//: buffer keyed by frame number and the consumer takes them in sequence,
//: which keeps every shoebox decision on one thread and needs no locking
//: around the boxes at all. HDF5 cannot be entered from two threads at once,
//: so each worker holds its own Reader and the library's own mutex serialises
//: the fetch; the decompression, which is the expensive part, runs outside it.
//:
//: The buffer is bounded because a decompressed frame of a 16M detector is 72
//: MB and running ahead without limit would be a way to exhaust memory
//: instead of time.
class FrameQueue {
 public:
  struct Decoded {
    std::int64_t number = -1;
    unsigned bit_depth = 0;
    std::size_t height = 0;
    std::size_t width = 0;
    std::vector<std::uint8_t> pixels;
  };

  FrameQueue(std::size_t depth) : depth_(std::max<std::size_t>(depth, 1)) {}

  //: Called by a worker when a frame is ready. Never blocks.
  //:
  //: Bounding the buffer here deadlocks: a worker holding frame zero waits for
  //: room while the buffer is full of frames one to eight, and the consumer
  //: waits for frame zero. The bound belongs on how far ahead work is handed
  //: out, not on how much has come back -- see `wait_for_room`, which cannot
  //: deadlock because keys are handed out in order and the frame the consumer
  //: is waiting for was therefore claimed before any frame ahead of it.
  void put(Decoded frame) {
    std::lock_guard<std::mutex> held(lock_);
    ready_.emplace(frame.number, std::move(frame));
    arrived_.notify_all();
  }

  //: Called by a worker before it starts on `which`: holds it back until the
  //: consumer is close enough behind.
  void wait_for_room(std::size_t which) {
    std::unique_lock<std::mutex> held(lock_);
    room_.wait(held, [&] { return which < consumed_ + depth_ || stop_; });
  }

  //: Called by the consumer when it has finished with a frame.
  void consumed(std::size_t which) {
    std::lock_guard<std::mutex> held(lock_);
    consumed_ = which + 1;
    room_.notify_all();
  }

  //: Called by a worker that read nothing for a key, so the consumer does not
  //: wait for a frame that will never come.
  void skip(std::int64_t number) {
    std::unique_lock<std::mutex> held(lock_);
    missing_.insert(number);
    arrived_.notify_all();
  }

  //: The frame with this number, once some worker has produced it. False when
  //: no worker ever will.
  bool take(std::int64_t number, Decoded *out) {
    std::unique_lock<std::mutex> held(lock_);
    arrived_.wait(held, [&] {
      return ready_.count(number) > 0 || missing_.count(number) > 0 || stop_;
    });
    if (stop_) return false;
    const auto found = ready_.find(number);
    if (found == ready_.end()) {
      missing_.erase(number);
      return false;
    }
    *out = std::move(found->second);
    ready_.erase(found);
    room_.notify_all();
    return true;
  }

  void halt() {
    std::lock_guard<std::mutex> held(lock_);
    stop_ = true;
    arrived_.notify_all();
    room_.notify_all();
  }

 private:
  std::mutex lock_;
  std::condition_variable arrived_;
  std::condition_variable room_;
  std::map<std::int64_t, Decoded> ready_;
  std::set<std::int64_t> missing_;
  std::size_t depth_;
  std::size_t consumed_ = 0;
  bool stop_ = false;
};

void usage(const char *program) {
  std::printf(
      "usage: %s [options] EXPT [STRONG_REFL]\n"
      "\n"
      "  -o FILE           where to write (integrated.refl)\n"
      "  --images PATH     the image file; by default the .expt's own\n"
      "                    imageset template is used\n"
      "  --sigma-b B --sigma-m M   profile model. Estimated from STRONG_REFL\n"
      "                    if given, else taken from EXPT's profile block\n"
      "  --n-sigma N       foreground spans plus and minus N sigma (3)\n"
      "  --box-scale S     box is S times wider than the foreground (1.9)\n"
      "  --min-zeta Z      skip reflections nearer the rotation axis than this\n"
      "                    (0.05). A reflection's extent in phi goes as 1/|zeta|,\n"
      "                    so the smallest are forty images deep\n"
      "  --d-min D         resolution limit\n"
      "  --first-image N --last-image N   restrict to part of the scan\n"
      "  --gain G          detector gain, counts per photon (1)\n"

      "  --save-shoeboxes  keep the pixels and the mask in the output\n"
      "  --threads N       threads fetching and decompressing frames; 0 is one\n"
      "                    per core, 1 is none (0)\n"
      "  --window N        frames per window (1000)\n"
      "  --max-boxes N     shoeboxes per window (4000); about 31 kB each\n"
      "  --summation-only  skip profile fitting and its second pass\n"
      "  --grid-points N   the profile grid is 2N+1 a side (4)\n"
      "  --subdivisions N  pixel subdivisions per axis (5, as Kabsch uses)\n"
      "  --regions N       detector divided N by N for reference profiles (3)\n"
      "  --scan-blocks N   the scan divided N ways as well (5)\n"
      "  --reference-signal S   learn from reflections above S sigma (10)\n"
      "  --save-profiles F the learned reference profiles, as text\n"
      "  --timing          where the time went, by phase\n",
      program);
}

//: One frame, decompressed into counts.
struct Frame {
  std::vector<std::int32_t> pixels;
  std::size_t fast = 0;
  std::size_t slow = 0;
};

}  // namespace

int main(int argc, char **argv) {
  const std::set<std::string> known = {
      "-o",          "--sigma-b",    "--sigma-m",   "--n-sigma",
      "--box-scale", "--d-min",      "--first-image", "--last-image",
      "--gain",      "--save-shoeboxes", "--images", "--timing", "--threads", "--window", "--max-boxes", "--grid-points", "--subdivisions", "--regions", "--scan-blocks", "--reference-signal", "--summation-only", "--save-profiles", "--min-zeta"};
  std::set<std::string> takes_value = known;
  takes_value.erase("--save-shoeboxes");
  takes_value.erase("--timing");
  takes_value.erase("--summation-only");
  const Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (args.help) {
    usage(argv[0]);
    return 0;
  }
  if (!args.ok) {
    std::fprintf(stderr, "mxi_integrate: %s\n", args.error.c_str());
    return 2;
  }
  if (args.positional.empty() || args.positional.size() > 2) {
    std::fprintf(stderr,
                 "mxi_integrate: expected an .expt and optionally a .refl of "
                 "strong spots\n");
    usage(argv[0]);
    return 2;
  }
  const std::string strong_path =
      args.positional.size() == 2 ? args.positional[1] : std::string();

  try {
    const double t_start = now_wall();
    double t_profile = 0.0, t_predict = 0.0, t_boxes = 0.0;
    double t_fetch = 0.0, t_decompress = 0.0, t_fill = 0.0, t_open = 0.0;
    double t_region = 0.0;
    double t_integrate = 0.0, t_write = 0.0;

    const ExperimentList experiments = read_experiments(args.positional[0]);
    if (experiments.size() == 0) {
      std::fprintf(stderr, "mxi_integrate: no experiments\n");
      return 1;
    }
    const Experiment &e = experiments[0];
    const Panel &panel = e.detector[0];

    const double t_profile_start = now_wall();
    MaskOptions mask_options;
    mask_options.box_scale = args.number("--box-scale", 1.9);
    // How near the rotation axis a reflection may be and still be integrated.
    //
    // This is the lever on the near-axis reflections, and they are worth a
    // word. A reflection's extent in phi is n_sigma sigma_m / |zeta|, so at the
    // default cut of 0.05 a box is forty images deep and its foreground is two
    // thousand voxels of which a handful hold the reflection. Those are the
    // whole of the summation disagreement with DIALS: setting them aside takes
    // the correlation from 0.939 to 1.000.
    //
    // They also cost time twice over, since a window has to read every frame
    // its longest reflection reaches, which is what drives the re-read rate.
    mask_options.min_zeta = args.number("--min-zeta", 0.05);

    // The profile model, in order of preference: what the command line says,
    // then what this package estimates from the strong spots, then what the
    // .expt was carrying. Estimating is the default because the numbers should
    // not have to be carried by hand between two programs, and because these
    // estimates are the ones the rest of this package was measured with.
    double sigma_b = args.number("--sigma-b", 0.0);
    double sigma_m = args.number("--sigma-m", 0.0);
    double n_sigma = args.number("--n-sigma", 0.0);
    const char *source = "the command line";
    if ((!(sigma_b > 0.0) || !(sigma_m > 0.0)) && !strong_path.empty()) {
      const Table strong = read_reflections(strong_path);
      const std::vector<Shoebox> strong_boxes = decode_shoeboxes(strong);
      if (strong_boxes.empty()) {
        std::fprintf(stderr,
                     "mxi_integrate: %s has no shoeboxes, so the profile model "
                     "cannot be estimated from it\n", strong_path.c_str());
        return 1;
      }
      const Column &s1_in = strong.at("s1");
      const Column &cal_in = strong.at("xyzcal.mm");
      std::vector<Shoebox> selected;
      std::vector<Vec3> beams;
      std::vector<RangeSample> samples;
      for (std::size_t i = 0; i < strong.nrows && i < strong_boxes.size(); ++i) {
        if (strong.has("flags") &&
            (strong.at("flags").integer(i) & flag::kUsedInRefinement) == 0) {
          continue;
        }
        const Vec3 beam{s1_in.real(i, 0), s1_in.real(i, 1), s1_in.real(i, 2)};
        if (std::fabs(compute_zeta(e, beam)) < 0.05) continue;
        if (!has_prediction(strong, i)) continue;
        selected.push_back(strong_boxes[i]);
        beams.push_back(beam);
        for (const RangeSample &sample :
             range_samples(e, strong_boxes[i], cal_in.real(i, 2),
                           compute_zeta(e, beam))) {
          samples.push_back(sample);
        }
      }
      std::size_t used = 0;
      if (!(sigma_b > 0.0)) sigma_b = beam_divergence(e, selected, beams, &used);
      if (!(sigma_m > 0.0)) {
        sigma_m = reflecting_range(samples, Scan::radians(e.scan.osc_width), 0.0);
      }
      source = "estimated from the strong spots";
    }
    if ((!(sigma_b > 0.0) || !(sigma_m > 0.0)) && experiments.profile.present) {
      if (!(sigma_b > 0.0)) sigma_b = experiments.profile.sigma_b;
      if (!(sigma_m > 0.0)) sigma_m = experiments.profile.sigma_m;
      if (!(n_sigma > 0.0)) n_sigma = experiments.profile.n_sigma;
      source = "the profile block of the .expt";
    }
    if (!(sigma_b > 0.0) || !(sigma_m > 0.0)) {
      std::fprintf(stderr,
                   "mxi_integrate: no profile model. Give a .refl of strong "
                   "spots with shoeboxes to estimate it from, or --sigma-b and "
                   "--sigma-m, or an .expt that carries a profile block.\n");
      return 1;
    }
    mask_options.sigma_d = sigma_b;
    mask_options.sigma_m = sigma_m;
    mask_options.n_sigma = n_sigma > 0.0 ? n_sigma : 3.0;
    t_profile = now_wall() - t_profile_start;
    std::printf("profile model: sigma_b %.6f sigma_m %.6f n_sigma %.1f (%s)\n",
                sigma_b, sigma_m, mask_options.n_sigma, source);
    IntegrateOptions integrate_options;
    integrate_options.gain = args.number("--gain", 1.0);

    // The images. series::nxmx unpacks the NXmx virtual dataset so each
    // frame's compressed bytes can be read without HDF5's filter pipeline;
    // decompress::image then turns them into pixels.
    // Where the images are: the .expt says, and --images overrides it. The
    // experiment list already carries the path and being told it twice is how
    // the two come to disagree -- but a file moved since it was written needs
    // the override, so both exist and the program says which it used.
    std::string image_path = args.value("--images", "");
    const char *image_source = "--images";
    if (image_path.empty()) {
      image_path = experiments.image_template;
      image_source = "the .expt's imageset template";
      if (image_path.empty()) {
        std::fprintf(stderr,
                     "mxi_integrate: %s has no imageset template, so there is "
                     "nothing to read images from; give --images\n",
                     args.positional[0].c_str());
        return 1;
      }
      if (image_path.find('#') != std::string::npos) {
        std::fprintf(stderr,
                     "mxi_integrate: the imageset template is %s, which names "
                     "a numbered series rather than one file; give --images\n",
                     image_path.c_str());
        return 1;
      }
    }
    std::printf("images: %s (from %s)\n", image_path.c_str(), image_source);

    std::unique_ptr<series::Series> images = series::nxmx(image_path);
    series::Info info;
    if (!images || !images->try_open(&info)) {
      std::fprintf(stderr,
                   "mxi_integrate: cannot open %s. If the data have moved "
                   "since the .expt was written, give --images\n",
                   image_path.c_str());
      return 1;
    }
    std::printf("%s\n", images->describe().c_str());

    std::size_t workers = static_cast<std::size_t>(args.number("--threads", 0.0));
    if (workers == 0) {
      workers = std::thread::hardware_concurrency();
      if (workers == 0) workers = 1;
    }

    PredictOptions predict_options;
    predict_options.d_min = args.number("--d-min", 0.0);
    // Prediction was 39.7 per cent of a large integration and all of it on one
    // thread. One unit of work per h, joined in h order, so the list is the
    // same whatever the thread count.
    predict_options.threads = workers;
    const double t_predict_start = now_wall();
    const std::vector<Prediction> predictions = predict(e, predict_options);
    t_predict = now_wall() - t_predict_start;
    const double first_image = args.number("--first-image", 0.0);
    const double last_image =
        args.number("--last-image", static_cast<double>(e.scan.num_images()));
    std::printf("%zu reflections predicted\n", predictions.size());

    // TWO PASSES, AND WHY
    //
    // Every shoebox cannot exist at once. Ten rotations of insulin is 7.4
    // million reflections and, at about 3500 voxels each, 233 GB of pixels. So
    // the first pass works out every bounding box and keeps nothing but the
    // box, and the second walks the scan in blocks of images, holding only the
    // shoeboxes of the block it is on.
    //
    // A reflection belongs to the block its FIRST image falls in, and a block
    // reads whatever frames its reflections span -- which is a few images past
    // its own end. Splitting a reflection across two blocks would integrate
    // half of it twice, so the blocks overlap in what they read and never in
    // what they own.
    struct Planned {
      const Prediction *prediction;
      std::int32_t bbox[6];
    };
    const double t_boxes_start = now_wall();
    std::vector<Planned> planned;
    std::map<std::string, std::size_t> refused;
    std::size_t outside_range = 0;
    for (const Prediction &p : predictions) {
      if (p.z < first_image || p.z > last_image) {
        ++outside_range;
        continue;
      }
      Planned item;
      item.prediction = &p;
      BoxRejection why = BoxRejection::kNone;
      if (!integration_bbox(e, p, mask_options, item.bbox, &why)) {
        ++refused[describe(why)];
        continue;
      }
      planned.push_back(item);
    }
    t_boxes = now_wall() - t_boxes_start;
    std::printf("%zu shoeboxes to fill\n", planned.size());
    if (outside_range > 0) {
      std::printf("  %zu outside the image range\n", outside_range);
    }
    for (const auto &entry : refused) {
      std::printf("  %zu %s\n", entry.second, entry.first.c_str());
    }
    if (planned.empty()) return 1;

    // In order of first image, so a block is a contiguous run of this vector.
    std::sort(planned.begin(), planned.end(),
              [](const Planned &a, const Planned &b) {
                return a.bbox[4] < b.bbox[4];
              });

    // ONE PASS OVER THE FRAMES, WITH SHOEBOXES OPEN ACROSS THEM
    //
    // Blocks of images do not work here. A reflection near the rotation axis
    // spans n_sigma sigma_m / |zeta| in phi, which at the zeta cut is about a
    // hundred and fifty frames, and one of those in a block forces the whole
    // block to read that far -- and the next block to read it again. On three
    // hundred frames with blocks of ten that was 2951 reads of 300 frames.
    //
    // So the frames are walked once in order instead. A shoebox opens when the
    // frame it starts on comes round, takes a slice from every frame it spans,
    // and is integrated and released on the frame it ends. Every frame is read
    // exactly once and memory is bounded by what is open at the time, which is
    // what the block size was trying to bound anyway.
    Table out;
    out.nrows = planned.size();
    Column &miller = out.int_column("miller_index", "cctbx::miller::index<>", 3);
    Column &panel_column = out.int_column("panel", "std::size_t", 1);
    Column &id = out.int_column("id", "int", 1);
    Column &imageset = out.int_column("imageset_id", "int", 1);
    Column &flags = out.int_column("flags", "std::size_t", 1);
    Column &entering = out.int_column("entering", "bool", 1);
    Column &bbox = out.int_column("bbox", "int6", 6);
    Column &n_fg = out.int_column("num_pixels.foreground", "int", 1);
    Column &n_bg = out.int_column("num_pixels.background", "int", 1);
    Column &n_val = out.int_column("num_pixels.valid", "int", 1);
    Column &cal_px = out.real_column("xyzcal.px", "vec3<double>", 3);
    Column &cal_mm = out.real_column("xyzcal.mm", "vec3<double>", 3);
    Column &s1_column = out.real_column("s1", "vec3<double>", 3);
    Column &isum = out.real_column("intensity.sum.value", "double", 1);
    Column &ivar = out.real_column("intensity.sum.variance", "double", 1);
    Column &bmean = out.real_column("background.mean", "double", 1);
    Column &bsum = out.real_column("background.sum.value", "double", 1);
    Column &bsumvar = out.real_column("background.sum.variance", "double", 1);
    Column &qe_column = out.real_column("qe", "double", 1);
    Column &lp_column = out.real_column("lp", "double", 1);
    Column &obs_px = out.real_column("xyzobs.px.value", "vec3<double>", 3);
    Column &obs_px_var =
        out.real_column("xyzobs.px.variance", "vec3<double>", 3);
    Column &obs_mm = out.real_column("xyzobs.mm.value", "vec3<double>", 3);
    Column &d_column = out.real_column("d", "double", 1);
    Column &iprf = out.real_column("intensity.prf.value", "double", 1);
    Column &iprf_var = out.real_column("intensity.prf.variance", "double", 1);
    Column &prf_cc = out.real_column("profile.correlation", "double", 1);
    Column &zeta_column = out.real_column("zeta", "double", 1);
    Column &part_column = out.real_column("partiality", "double", 1);
    Column &partial_id = out.int_column("partial_id", "std::size_t", 1);
    Column &n_bg_used = out.int_column("num_pixels.background_used", "int", 1);
    Column &obs_mm_var =
        out.real_column("xyzobs.mm.variance", "vec3<double>", 3);

    const bool save = args.has("--save-shoeboxes");
    std::string shoebox_bytes;
    // Saved shoeboxes come out in the order they close, which is not the order
    // of the table, so they are held by row and encoded at the end.
    std::vector<Shoebox> saved;
    if (save) saved.resize(planned.size());

    std::unique_ptr<series::Reader> reader = images->reader();
    // The keys, in the order the series gives them. Reading each one to build
    // a frame-number-to-key map first would read every chunk twice, which on
    // the thirty degree sweep was most of the run: 18 of 21 seconds spent
    // fetching compressed bytes that were then fetched again.
    //
    // The loop below is driven by the frames as they arrive instead, and
    // checks they arrive in order rather than assuming it -- the Series
    // interface promises keys that can be read, not keys in sequence.
    const std::vector<std::string> keys = images->ready();

    series::Frame frame;
    std::vector<std::uint8_t> pixels;
    std::size_t frames_read = 0;
    std::size_t frames_missing = 0;
    // Which frames any shoebox wanted, so "frames read" can be compared with
    // something rather than left to be compared with the image count.
    std::set<std::int32_t> wanted;
    std::size_t bad_pixels = 0;
    std::size_t integrated = 0;
    std::size_t most_open = 0;
    const Vec3 s0 = e.beam.s0();
    const Vec3 axis = e.goniometer.lab_axis();

    // The shoeboxes currently open, by their row in the table.

    const auto close = [&](std::size_t row, Shoebox *box) -> bool {
      const Prediction &p = *planned[row].prediction;
      const IntegratedReflection r = integrate_shoebox(box, integrate_options);
      miller.ints[row * 3 + 0] = p.h;
      miller.ints[row * 3 + 1] = p.k;
      miller.ints[row * 3 + 2] = p.l;
      panel_column.ints[row] = static_cast<std::int64_t>(p.panel);
      id.ints[row] = 0;
      imageset.ints[row] = 0;
      flags.ints[row] = flag::kPredicted | (r.valid ? flag::kIntegratedSum : 0);
      entering.ints[row] = p.s1.dot(axis.cross(s0)) > 0.0 ? 1 : 0;
      for (int k = 0; k < 6; ++k) bbox.ints[row * 6 + k] = box->bbox[k];
      n_fg.ints[row] = static_cast<std::int64_t>(r.n_foreground);
      n_bg.ints[row] = static_cast<std::int64_t>(r.n_background);
      n_val.ints[row] = static_cast<std::int64_t>(r.n_valid);
      cal_px.reals[row * 3 + 0] = p.px_fast;
      cal_px.reals[row * 3 + 1] = p.px_slow;
      cal_px.reals[row * 3 + 2] = p.z;
      const auto mm = panel.px_to_mm(p.px_fast, p.px_slow);
      cal_mm.reals[row * 3 + 0] = mm.first;
      cal_mm.reals[row * 3 + 1] = mm.second;
      cal_mm.reals[row * 3 + 2] = p.phi;
      for (int k = 0; k < 3; ++k) s1_column.reals[row * 3 + k] = p.s1[k];
      isum.reals[row] = r.intensity;
      ivar.reals[row] = r.variance;
      bmean.reals[row] = r.background_mean;
      bsum.reals[row] = r.background_sum;
      bsumvar.reals[row] = r.background_sum_variance;
      qe_column.reals[row] = quantum_efficiency(panel, p.s1);
      lp_column.reals[row] = lorentz_polarization(e.beam, e.goniometer, p.s1);
      d_column.reals[row] = e.crystal ? resolution(*e.crystal, p.h, p.k, p.l) : 0.0;
      const double z_of = compute_zeta(e, p.s1);
      zeta_column.reals[row] = z_of;
      part_column.reals[row] = partiality(e.scan, p.phi, z_of,
                                         mask_options.sigma_m, box->bbox[4],
                                         box->bbox[5]);
      // Every reflection here is its own, since nothing splits one across two
      // rows; DIALS uses this to tie the pieces of a split reflection together.
      partial_id.ints[row] = static_cast<std::int64_t>(row);
      // Everything not foreground was used: this integrator has no second
      // round of rejection on top of the GLM's own weighting.
      n_bg_used.ints[row] = static_cast<std::int64_t>(r.n_background);
      // The centroid variance in millimetres and radians. The px values it
      // comes from do not reproduce DIALS' and neither will these.
      obs_mm_var.reals[row * 3 + 0] =
          r.centroid_variance_fast * panel.pixel_size[0] * panel.pixel_size[0];
      obs_mm_var.reals[row * 3 + 1] =
          r.centroid_variance_slow * panel.pixel_size[1] * panel.pixel_size[1];
      obs_mm_var.reals[row * 3 + 2] =
          r.centroid_variance_z * Scan::radians(e.scan.osc_width) *
          Scan::radians(e.scan.osc_width);
      const double of = r.centroid_valid ? r.centroid_fast : p.px_fast;
      const double os = r.centroid_valid ? r.centroid_slow : p.px_slow;
      const double oz = r.centroid_valid ? r.centroid_z : p.z;
      obs_px.reals[row * 3 + 0] = of;
      obs_px.reals[row * 3 + 1] = os;
      obs_px.reals[row * 3 + 2] = oz;
      obs_px_var.reals[row * 3 + 0] = r.centroid_variance_fast;
      obs_px_var.reals[row * 3 + 1] = r.centroid_variance_slow;
      obs_px_var.reals[row * 3 + 2] = r.centroid_variance_z;
      const auto obs_mm_pair = panel.px_to_mm(of, os);
      obs_mm.reals[row * 3 + 0] = obs_mm_pair.first;
      obs_mm.reals[row * 3 + 1] = obs_mm_pair.second;
      obs_mm.reals[row * 3 + 2] = e.scan.phi_from_z(oz);
      return r.valid;
    };

    // WINDOWS, AND WHY NOT A PIPELINE
    //
    // Threading only the reads settles at about one core in use. The thread
    // that owns the shoeboxes has 111 of the 310 seconds of work on it --
    // building masks, filling, background and summation -- so the readers
    // finish their lookahead and wait. No amount of reader threads fixes that.
    //
    // The in-order constraint that forced the pipeline is not actually
    // needed. A frame writes only its own z plane of a shoebox, so two frames
    // of the same box can be filled by two threads at once without touching
    // the same double. What needs ordering is nothing; what needs a box to
    // exist is everything.
    //
    // So the scan is cut into windows of frames. Within a window: the boxes
    // are built in parallel, the frames are fetched, decompressed and filled
    // in parallel, and the boxes are integrated in parallel. Nothing is
    // serialised except the window boundary.
    //
    // A box belongs to the window its first frame falls in and the window
    // reads past its own end to finish it, so a window wants to be much
    // longer than a shoebox: at a thousand frames the overhead is under one
    // per cent, and holds about ten thousand boxes, which is 320 MB.
    const std::size_t window =
        static_cast<std::size_t>(std::max(1.0, args.number("--window", 1000.0)));
    // A safety net, not the primary bound. The window has to be long compared
    // with a shoebox or the frames it reads past its own end dominate: at four
    // thousand boxes a three hundred frame sweep becomes six windows of fifty
    // frames and reads 748 of them instead of 321, and decompression goes from
    // 1.13 seconds to 2.66. At twenty thousand -- about 620 MB in flight --
    // the frame count is back to 321.
    const std::size_t max_boxes =
        static_cast<std::size_t>(std::max(1.0, args.number("--max-boxes", 20000.0)));
    std::printf("%zu threads, windows of %zu frames or %zu shoeboxes\n", workers,
                window, max_boxes);

    //: Run `count` units of work over the pool, by index.
    const auto in_parallel = [&](std::size_t count, auto &&body) {
      if (count == 0) return;
      const std::size_t n = std::min(workers, count);
      if (n <= 1) {
        for (std::size_t i = 0; i < count; ++i) body(i);
        return;
      }
      std::atomic<std::size_t> next_unit{0};
      std::vector<std::thread> pool;
      pool.reserve(n - 1);
      const auto run = [&]() {
        for (;;) {
          const std::size_t i = next_unit.fetch_add(1);
          if (i >= count) break;
          body(i);
        }
      };
      for (std::size_t t = 1; t < n; ++t) pool.emplace_back(run);
      run();
      for (std::thread &t : pool) t.join();
    };

    // TWO PASSES OVER THE IMAGES
    //
    // The reference profiles are learned from the reflections themselves, so
    // they cannot exist until something has been integrated. The images are
    // therefore read twice: once to sum and to learn, once to fit. That is the
    // slow way round and it is the right one -- a profile learned from part of
    // a scan and applied to the rest is a different algorithm, and one whose
    // errors would be hard to attribute.
    GridSpec grid_spec;
    grid_spec.n = static_cast<int>(args.number("--grid-points", 4.0));
    grid_spec.sigma_d = sigma_b;
    grid_spec.sigma_m = sigma_m;
    grid_spec.half_width = mask_options.n_sigma;
    grid_spec.subdivisions = static_cast<int>(args.number("--subdivisions", 5.0));
    ReferenceProfiles reference = make_reference(
        grid_spec, static_cast<int>(args.number("--regions", 3.0)),
        static_cast<int>(args.number("--scan-blocks", 5.0)), e.detector.size(),
        first_image, last_image);
    const bool fitting = !args.has("--summation-only");
    // Which reflections are worth learning from: strong, nearly whole, and
    // mostly inside the grid. DIALS marks these `reference_spot`.
    const double least_signal = args.number("--reference-signal", 10.0);

    std::vector<Shoebox> boxes;
    std::size_t at = 0;
    int pass = 0;
    std::size_t references_used = 0;
    std::size_t fitted = 0;
    double t_transform = 0.0, t_fit = 0.0;
    while (at < planned.size()) {
      const std::int32_t window_start = planned[at].bbox[4];
      const std::int32_t window_end =
          window_start + static_cast<std::int32_t>(window);
      // Bounded by frames AND by boxes. Frames alone is not enough: a window
      // of a thousand frames swallows a three hundred frame sweep whole, which
      // is 21032 boxes and 650 MB, and the allocation costs more than the
      // parallelism saves -- measured at 16.7 seconds against 4.5 for the
      // version this replaced.
      std::size_t stop = at;
      std::int32_t highest = window_start;
      while (stop < planned.size() && planned[stop].bbox[4] < window_end &&
             stop - at < max_boxes) {
        highest = std::max(highest, planned[stop].bbox[5]);
        ++stop;
      }
      const std::size_t count = stop - at;

      // The boxes, built in parallel: this is the mask, and it was 46 seconds.
      const double t_open_start = now_wall();
      boxes.assign(count, Shoebox{});
      in_parallel(count, [&](std::size_t i) {
        BoxRejection ignored = BoxRejection::kNone;
        if (build_shoebox(e, *planned[at + i].prediction, mask_options,
                          &boxes[i], &ignored)) {
          boxes[i].data.assign(boxes[i].size(), 0.0f);
        }
      });
      t_open += now_wall() - t_open_start;
      most_open = std::max(most_open, count);

      // Which boxes each frame touches, so a worker filling a frame knows
      // what to write without searching every box in the window.
      std::map<std::int32_t, std::vector<std::size_t>> touching;
      for (std::size_t i = 0; i < count; ++i) {
        for (std::int32_t z = boxes[i].bbox[4]; z < boxes[i].bbox[5]; ++z) {
          touching[z].push_back(i);
        }
      }
      std::vector<std::int32_t> frame_numbers;
      frame_numbers.reserve(touching.size());
      for (const auto &entry : touching) frame_numbers.push_back(entry.first);

      // Fetch, decompress and fill, in parallel over frames. Each frame writes
      // only its own z plane of each box it touches, so two frames of one box
      // never touch the same voxel.
      std::atomic<std::size_t> frames_done{0};
      std::atomic<std::size_t> bad_here{0};
      std::atomic<std::size_t> unread{0};
      std::vector<double> fetch_by_thread(std::max<std::size_t>(workers, 1), 0.0);
      std::vector<double> decompress_by_thread(fetch_by_thread.size(), 0.0);
      std::vector<double> fill_by_thread(fetch_by_thread.size(), 0.0);
      std::atomic<std::size_t> thread_slot{0};
      const double t_region_start = now_wall();
      in_parallel(frame_numbers.size(), [&](std::size_t which) {
        thread_local std::unique_ptr<series::Reader> mine;
        thread_local std::size_t slot = 0;
        thread_local bool first = true;
        if (first) {
          mine = images->reader();
          slot = thread_slot.fetch_add(1) % fetch_by_thread.size();
          first = false;
        }
        const std::int32_t z = frame_numbers[which];
        if (z < 0 || static_cast<std::size_t>(z) >= keys.size()) return;
        series::Frame raw;
        const double t0 = now_wall();
        if (!mine->read(keys[static_cast<std::size_t>(z)], &raw)) {
          // A frame the writer never received: an unallocated chunk. Counted,
          // because a shoebox that spans it is missing a slice and will
          // integrate low, and silently dropping it leaves "frames read" less
          // than the number of images with no explanation.
          unread.fetch_add(1);
          return;
        }
        const double t1 = now_wall();
        const std::size_t height = static_cast<std::size_t>(raw.height);
        const std::size_t width = static_cast<std::size_t>(raw.width);
        const std::size_t bytes =
            decompress::frame_bytes(height, width, raw.bit_depth);
        std::vector<std::uint8_t> pixels(bytes);
        decompress::image(raw.data, raw.algorithm, raw.bit_depth, height, width,
                          {pixels.data(), bytes});
        const double t2 = now_wall();
        frames_done.fetch_add(1);

        const auto fill = [&](auto typed, std::uint32_t bad) {
          using Pixel = decltype(typed);
          const Pixel *const raw_pixels =
              reinterpret_cast<const Pixel *>(pixels.data());
          std::size_t bad_count = 0;
          for (std::size_t i : touching[z]) {
            Shoebox &box = boxes[i];
            if (box.data.empty()) continue;
            const std::int32_t zi = z - box.bbox[4];
            for (std::int32_t y = 0; y < box.ny(); ++y) {
              const std::size_t row =
                  static_cast<std::size_t>(box.bbox[2] + y) * width;
              for (std::int32_t x = 0; x < box.nx(); ++x) {
                const std::size_t index =
                    row + static_cast<std::size_t>(box.bbox[0] + x);
                const std::size_t into = box.at(x, y, zi);
                const std::uint32_t v =
                    static_cast<std::uint32_t>(raw_pixels[index]);
                // The largest representable value is the bad-pixel marker, not
                // a count: excluded from both sums rather than counted as
                // zero, which would drag the background down wherever a module
                // gap crosses a shoebox.
                if (v == bad) {
                  box.mask[into] = 0;
                  ++bad_count;
                } else {
                  box.data[into] = static_cast<float>(v);
                }
              }
            }
          }
          bad_here.fetch_add(bad_count);
        };
        if (raw.bit_depth == 16) {
          fill(std::uint16_t{}, 0xFFFFu);
        } else if (raw.bit_depth == 32) {
          fill(std::uint32_t{}, 0xFFFFFFFFu);
        } else {
          throw std::runtime_error("unsupported bit depth " +
                                   std::to_string(raw.bit_depth));
        }
        const double t3 = now_wall();
        fetch_by_thread[slot] += t1 - t0;
        decompress_by_thread[slot] += t2 - t1;
        fill_by_thread[slot] += t3 - t2;
      });
      t_region += now_wall() - t_region_start;
      frames_read += frames_done.load();
      bad_pixels += bad_here.load();
      frames_missing += unread.load();
      for (const auto &entry : touching) wanted.insert(entry.first);
      // Thread-seconds, summed over threads, NOT the busiest thread's share.
      //
      // Reporting the busiest thread per phase gave fetching 77.1 per cent and
      // decompressing 47.4 per cent of the same run: both are wall clock
      // inside one parallel region, they overlap, and the percentages summed
      // to 124. Thread-seconds are additive and comparable with each other,
      // and the region's own wall clock is reported beside them so the
      // difference between work and waiting is visible.
      const auto summed = [](const std::vector<double> &v) {
        double all = 0.0;
        for (double x : v) all += x;
        return all;
      };
      t_fetch += summed(fetch_by_thread);
      t_decompress += summed(decompress_by_thread);
      t_fill += summed(fill_by_thread);

      if (pass == 0) {
        // Integrate, in parallel over boxes.
        const double t_close_start = now_wall();
        std::atomic<std::size_t> done{0};
        in_parallel(count, [&](std::size_t i) {
          if (boxes[i].data.empty()) return;
          if (close(at + i, &boxes[i])) done.fetch_add(1);
        });
        integrated += done.load();
        t_integrate += now_wall() - t_close_start;

        // Learn from the ones worth learning from. Serial, because it
        // accumulates into shared profiles and is a small part of the cost;
        // threading it would need a profile per thread and a reduction, which
        // is work to do when it shows up in the timing and not before.
        if (fitting) {
          const double t0 = now_wall();
          for (std::size_t i = 0; i < count; ++i) {
            if (boxes[i].data.empty()) continue;
            const std::size_t row = at + i;
            const double signal = isum.reals[row];
            const double sigma = std::sqrt(std::max(ivar.reals[row], 1e-12));
            if (!(signal > least_signal * sigma)) continue;
            if (part_column.reals[row] < 0.99) continue;
            const Prediction &q = *planned[row].prediction;
            const Transformed t =
                transform_shoebox(e, boxes[i], q.s1, q.phi, grid_spec);
            if (!t.valid || t.outside > 0.05) continue;
            const std::size_t region =
                reference.region_of(panel, static_cast<std::size_t>(q.panel),
                                    q.px_fast, q.px_slow, q.z);
            if (add_reference(&reference, region, t)) ++references_used;
          }
          t_transform += now_wall() - t0;
        }
      } else {
        // Fit, in parallel over boxes: each writes only its own row.
        const double t0 = now_wall();
        std::atomic<std::size_t> done{0};
        in_parallel(count, [&](std::size_t i) {
          if (boxes[i].data.empty()) return;
          const std::size_t row = at + i;
          const Prediction &q = *planned[row].prediction;
          // The background the GLM found in the first pass, put back so the
          // fit subtracts the same thing the sum did.
          boxes[i].background.assign(boxes[i].size(),
                                     static_cast<float>(bmean.reals[row]));
          const Transformed t =
              transform_shoebox(e, boxes[i], q.s1, q.phi, grid_spec);
          if (!t.valid) return;
          // A weighted average of the nearby profiles, not the nearest one:
          // taking the nearest makes the model jump at a cell boundary, so two
          // reflections either side of one are fitted with different profiles.
          const std::vector<double> local =
              profile_at(reference, panel, static_cast<std::size_t>(q.panel),
                         q.px_fast, q.px_slow, q.z);
          // Fitted against the PIXELS, with the profile carried onto them,
          // rather than against the grid with the pixels carried onto it. The
          // two give the same intensity and very different variances: the grid
          // has more points than the shoebox has pixels and one pixel's counts
          // reach several of them, so treating its points as independent
          // overcounts the information. Measured against DIALS, the grid fit
          // claimed a variance 0.35 of the summed one at high resolution where
          // DIALS has 0.84 -- errors too small by 1.6, and everything weighted
          // by them wrong.
          const std::vector<double> on_pixels = profile_on_pixels(
              e, boxes[i], q.s1, q.phi, grid_spec, local);
          const ProfileFit fit =
              on_pixels.empty()
                  ? fit_profile(local, t, integrate_options.gain)
                  : fit_on_pixels(boxes[i], on_pixels, integrate_options.gain);
          if (!fit.valid) return;
          iprf.reals[row] = fit.intensity;
          iprf_var.reals[row] = fit.variance;
          prf_cc.reals[row] = fit.correlation;
          flags.ints[row] |= flag::kIntegratedPrf;
          done.fetch_add(1);
        });
        fitted += done.load();
        t_fit += now_wall() - t0;
      }

      if (save) {
        for (std::size_t i = 0; i < count; ++i) saved[at + i] = std::move(boxes[i]);
      }
      boxes.clear();
      at = stop;
      if (at >= planned.size() && pass == 0 && fitting) {
        // Between the passes: the profiles are what they are going to be.
        finalise_reference(&reference);
        const std::string profile_path = args.value("--save-profiles", "");
        if (!profile_path.empty()) {
          std::FILE *f = std::fopen(profile_path.c_str(), "w");
          if (f == nullptr) {
            std::fprintf(stderr, "mxi_integrate: cannot write %s\n",
                         profile_path.c_str());
          } else {
            // Plain text, one profile a block, because the thing that reads
            // this is a person with a plotting script and not a program that
            // needs a format.
            std::fprintf(f, "# reference profiles\n");
            std::fprintf(f, "side %d\n", grid_spec.side());
            std::fprintf(f, "sigma_b %.9g\n", grid_spec.sigma_d);
            std::fprintf(f, "sigma_m %.9g\n", grid_spec.sigma_m);
            std::fprintf(f, "half_width %.9g\n", grid_spec.half_width);
            std::fprintf(f, "divisions %d\n", reference.divisions);
            std::fprintf(f, "blocks %d\n", reference.blocks);
            std::fprintf(f, "panels %zu\n", reference.panels);
            for (std::size_t r = 0; r < reference.profile.size(); ++r) {
              std::fprintf(f, "profile %zu spots %zu\n", r, reference.spots[r]);
              for (double v : reference.profile[r]) std::fprintf(f, "%.9g\n", v);
            }
            std::fclose(f);
            std::printf("wrote %s\n", profile_path.c_str());
          }
        }
        std::size_t empty = 0;
        for (std::size_t r = 0; r < reference.spots.size(); ++r) {
          if (reference.spots[r] < 10) ++empty;
        }
        std::printf(
            "learned %zu reference profiles from %zu reflections (%zu of %zu "
            "regions borrowed the detector average)\n",
            reference.region_count() - empty, references_used, empty,
            reference.region_count());
        if (references_used == 0) {
          std::printf(
              "  nothing to learn from, so no profile fitting: raise "
              "--reference-signal or check the summation\n");
        } else {
          pass = 1;
          at = 0;
          continue;
        }
      }
    }
    if (save) shoebox_bytes = encode_shoeboxes(saved);
    std::printf("at most %zu shoeboxes in a window\n", most_open);
    if (fitting) std::printf("%zu of %zu profile fitted\n", fitted, planned.size());

    std::printf("%zu frames read (%zu wanted by a shoebox, each read %.2f "
                "times), %zu bad pixels masked\n",
                frames_read, wanted.size(),
                wanted.empty() ? 0.0
                               : static_cast<double>(frames_read) /
                                     static_cast<double>(wanted.size()),
                bad_pixels);
    if (frames_missing > 0) {
      std::printf(
          "  %zu reads found no frame: an unallocated chunk is a frame the "
          "writer never received, and a shoebox spanning one is missing a "
          "slice\n",
          frames_missing);
    }
    std::printf("%zu of %zu integrated\n", integrated, planned.size());

    if (save) {
      Table::Opaque column;
      column.type = "Shoebox<>";
      column.rows = planned.size();
      std::printf("shoeboxes kept: %.1f MB\n", shoebox_bytes.size() / 1e6);
      column.bytes = std::move(shoebox_bytes);
      out.set_opaque("shoebox", std::move(column));
    }
    if (!e.identifier.empty()) out.identifiers[0] = e.identifier;

    const std::string path = args.value("-o", "integrated.refl");
    const double t_write_start = now_wall();
    write_reflections(path, out);
    t_write = now_wall() - t_write_start;
    std::printf("wrote %s\n", path.c_str());

    if (args.has("--timing")) {
      const double total = now_wall() - t_start;
      const auto line = [&](const char *name, double seconds) {
        std::printf("  %-26s %8.3f s  %5.1f%%\n", name, seconds,
                    total > 0.0 ? 100.0 * seconds / total : 0.0);
      };
      std::printf("\ntiming\n");
      line("the profile model", t_profile);
      line("prediction", t_predict);
      line("bounding boxes", t_boxes);
      line("opening shoeboxes", t_open);
      line("reading frames (wall)", t_region);
      std::printf("    of which, in thread-seconds over %zu threads:\n",
                  workers);
      const auto thread_line = [&](const char *name, double seconds) {
        std::printf("      %-22s %8.3f s  %5.2f x wall\n", name, seconds,
                    t_region > 0.0 ? seconds / t_region : 0.0);
      };
      thread_line("fetching", t_fetch);
      thread_line("decompressing", t_decompress);
      thread_line("filling shoeboxes", t_fill);
      line("background and summation", t_integrate);
      line("learning profiles", t_transform);
      line("profile fitting", t_fit);
      line("writing", t_write);
      std::printf("  %-26s %8.3f s\n", "total", total);
    }
    return 0;
  } catch (const std::exception &error) {
    std::fprintf(stderr, "mxi_integrate: %s\n", error.what());
    return 1;
  }
}
