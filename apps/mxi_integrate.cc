// mxi_integrate: summation integration on real images.
//
//   mxi_integrate integrated.expt master.nxs -o mine.refl --d-min 2.0
//   mxi_integrate ... --save-shoeboxes            # keep the pixels and the
//   mask
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
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include "../src/expt.hh"
#include "../src/refl.hh"
#include "args.hh"
#include "background.hh"
#include "integrate.hh"
#include "log_mirror.hh"
#include "mask.hh"
#include "postrefine.hh"
#include "predict.hh"
#include "profile_model.hh"
#include "reference.hh"
#include "shoebox.hh"
#include "summary.hh"

#include "decompress.hh"
#include "series.hh"

#include <memory>
#include <stdexcept>

namespace mxi {

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
    if (stop_)
      return false;
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
      "usage: %s [options] EXPT [INDEXED_REFL]\n"
      "\n"
      "  -o FILE           where to write (integrated.refl)\n"
      "  --images PATH     the image file; by default the .expt's own\n"
      "                    imageset template is used\n"
      "  --postrefine      integrate, refine against the centres integration\n"
      "                    measured, and integrate again with the refined\n"
      "                    models, which are written to --output-expt\n"
      "                    (integrated.expt). Removes the z offset the spot\n"
      "                    finder's centres leave in a refined model\n"
      "  --postrefine-points N   control points for its scan-varying pass;\n"
      "                    default one per 36 degrees and two more, at least "
      "5\n"
      "  --sigma-b B --sigma-m M   profile model. Estimated from INDEXED_REFL\n"
      "                    if given -- indexed or refined reflections, not a\n"
      "                    spot finder's -- else taken from EXPT's profile "
      "block\n"
      "  --n-sigma N       foreground spans plus and minus N sigma (3)\n"
      "  --box-scale S     box is S times wider than the foreground (1.9)\n"
      "  --min-zeta Z      skip reflections nearer the rotation axis than "
      "this\n"
      "                    (0.05). A reflection's extent in phi goes as "
      "1/|zeta|,\n"
      "                    so the smallest are forty images deep\n"
      "  --d-min D         resolution limit\n"
      "  --first-image N --last-image N   restrict to part of the scan\n"
      "  --gain G          detector gain, counts per photon (1)\n"

      "  --save-shoeboxes  keep the pixels and the mask in the output\n"
      "  --threads N       threads fetching and decompressing frames; 0 is "
      "one\n"
      "                    per core, 1 is none (0)\n"
      "  --window N        frames read per chunk (64); boxes stay open across\n"
      "                    chunks, so this bounds memory, not re-reading\n"
      "  --max-boxes N     boxes opened per chunk (20000); shortens the chunk\n"
      "                    when it bites, about 31 kB a box\n"
      "  --summation-only  skip profile fitting and its second pass\n"
      "  --grid-points N   the profile grid is 2N+1 a side (4)\n"
      "  --subdivisions N  pixel subdivisions per axis (5, as Kabsch uses)\n"
      "  --regions N       detector divided N by N for reference profiles (3)\n"
      "  --scan-blocks N   the scan divided N ways as well (5)\n"
      "  --reference-signal S   learn from reflections above S sigma (10)\n"
      "  --least-measured F  a fit needs this fraction of the reflection to\n"
      "                    have been recorded (0.6); below it the intensity "
      "is\n"
      "                    written but not flagged as fitted\n"
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

} // namespace

int run_program(int argc, char **argv);

// Integrate; refine against the centres that integration measured; integrate
// again with the refined models. Composed of whole runs of this program rather
// than of a function, because the integration is still one body in
// run_program: the models are written to --output-expt between the two, and
// the second run reads them as any later program would.
int integrate_with_postrefinement(const Arguments &args, const char *program,
                                  const std::set<std::string> &takes_value) {
  const std::string out_refl = args.value("-o", "integrated.refl");
  const std::string out_expt = args.value("--output-expt", "integrated.expt");
  const std::string first_refl = out_refl + ".before-postrefinement.refl";
  std::size_t points = 0;
  if (args.has("--postrefine-points")) {
    const double n = args.number("--postrefine-points", 0.0);
    if (!(n >= 2.0)) {
      std::fprintf(stderr,
                   "mxi_integrate: --postrefine-points needs at least 2\n");
      return 2;
    }
    points = static_cast<std::size_t>(n);
  }

  // A command line for one run: this one's options, less those that belong to
  // the post-refinement, with its own models and output. The first run keeps no
  // shoeboxes or profiles, since only the second run's are the result.
  const auto command = [&](const std::string &expt, const std::string &refl,
                           bool first) {
    std::vector<std::string> words = {program, expt};
    if (args.positional.size() == 2)
      words.push_back(args.positional[1]);
    for (const auto &[flag, value] : args.options) {
      if (flag == "--postrefine" || flag == "--postrefine-points" ||
          flag == "--output-expt" || flag == "-o")
        continue;
      if (first && (flag == "--save-shoeboxes" || flag == "--save-profiles"))
        continue;
      words.push_back(flag);
      if (takes_value.count(flag))
        words.push_back(value);
    }
    words.push_back("-o");
    words.push_back(refl);
    return words;
  };
  const auto run = [](std::vector<std::string> words) {
    std::vector<char *> argv;
    for (std::string &w : words)
      argv.push_back(w.data());
    argv.push_back(nullptr);
    return run_program(static_cast<int>(words.size()), argv.data());
  };

  std::printf("Integrating with the models as given, then refining against the "
              "centres that measures, then integrating again\n\n");
  std::printf("=== Integration with the models as given ===\n");
  int status = run(command(args.positional[0], first_refl, true));
  if (status != 0)
    return status;

  try {
    ExperimentList experiments = read_experiments(args.positional[0]);
    const Table first = read_reflections(first_refl);
    const PostrefineResult result = postrefine(experiments, first, points);
    std::printf("\n=== Post-refinement ===\n");
    std::printf(
        "Refining against %zu of the %zu integrated reflections: summed, "
        "and with a centre of mass\n",
        result.selected, result.candidates);
    const RefineResult &st = result.refinement.static_pass;
    const RefineResult &sv = result.refinement.varying_pass;
    if (st.n_used == 0) {
      std::fprintf(stderr,
                   "mxi_integrate: post-refinement fitted nothing, so there "
                   "is nothing to integrate again with\n");
      std::remove(first_refl.c_str());
      return 1;
    }
    std::printf(
        "  scan-static:  %3zu parameters, %zu reflections, RMSD %.4f px, "
        "%.4f px, %.4f images\n",
        st.n_parameters, st.n_used, st.rmsd_x, st.rmsd_y, st.rmsd_z);
    if (result.refinement.varied) {
      std::printf(
          "  scan-varying: %3zu parameters, %zu reflections, RMSD %.4f px, "
          "%.4f px, %.4f images; %zu control points, the detector held\n",
          sv.n_parameters, sv.n_used, sv.rmsd_x, sv.rmsd_y, sv.rmsd_z,
          experiments[0].crystal ? experiments[0].crystal->A_points.size() : 0);
    }
    write_experiments(out_expt, experiments);
    std::printf("Wrote the post-refined models to %s\n\n", out_expt.c_str());
  } catch (const std::exception &error) {
    std::fprintf(stderr, "mxi_integrate: post-refinement: %s\n", error.what());
    std::remove(first_refl.c_str());
    return 1;
  }
  std::remove(first_refl.c_str());

  std::printf("=== Integration with the post-refined models ===\n");
  return run(command(out_expt, out_refl, false));
}

int run_program(int argc, char **argv) {
  const std::set<std::string> known = {"-o",
                                       "--sigma-b",
                                       "--sigma-m",
                                       "--n-sigma",
                                       "--box-scale",
                                       "--d-min",
                                       "--first-image",
                                       "--last-image",
                                       "--gain",
                                       "--save-shoeboxes",
                                       "--images",
                                       "--timing",
                                       "--threads",
                                       "--window",
                                       "--max-boxes",
                                       "--grid-points",
                                       "--subdivisions",
                                       "--regions",
                                       "--scan-blocks",
                                       "--reference-signal",
                                       "--summation-only",
                                       "--save-profiles",
                                       "--min-zeta",
                                       "--least-measured",
                                       "--postrefine",
                                       "--postrefine-points",
                                       "--output-expt"};
  std::set<std::string> takes_value = known;
  takes_value.erase("--save-shoeboxes");
  takes_value.erase("--timing");
  takes_value.erase("--summation-only");
  takes_value.erase("--postrefine");
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
                 "indexed reflections\n");
    usage(argv[0]);
    return 2;
  }
  const std::string strong_path =
      args.positional.size() == 2 ? args.positional[1] : std::string();

  if (args.has("--postrefine"))
    return integrate_with_postrefinement(args, argv[0], takes_value);
  // Refused rather than ignored: a run that silently drops an option has used
  // settings nobody chose.
  for (const char *needs : {"--postrefine-points", "--output-expt"}) {
    if (args.has(needs)) {
      std::fprintf(stderr, "mxi_integrate: %s is for --postrefine\n", needs);
      return 2;
    }
  }

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
                     "cannot be estimated from it\n",
                     strong_path.c_str());
        return 1;
      }
      // The profile model needs INDEXED spots -- each spot's predicted
      // diffracted beam -- and a spot finder's strong.refl has none. Said in
      // those terms, because "no column 's1'" is true and names nothing a
      // user can do.
      if (!strong.has("s1") || !strong.has("xyzcal.mm")) {
        std::fprintf(
            stderr,
            "mxi_integrate: %s has no predictions for its spots, so the "
            "profile model cannot be estimated from it. Give the "
            "refined or indexed reflections (refined.refl), not the "
            "spot finder's, or set --sigma-b and --sigma-m\n",
            strong_path.c_str());
        return 1;
      }
      const Column &s1_in = strong.at("s1");
      const Column &cal_in = strong.at("xyzcal.mm");
      std::vector<Shoebox> selected;
      std::vector<Vec3> beams;
      std::vector<RangeSample> samples;
      for (std::size_t i = 0; i < strong.nrows && i < strong_boxes.size();
           ++i) {
        if (strong.has("flags") &&
            (strong.at("flags").integer(i) & flag::kUsedInRefinement) == 0) {
          continue;
        }
        const Vec3 beam{s1_in.real(i, 0), s1_in.real(i, 1), s1_in.real(i, 2)};
        if (std::fabs(compute_zeta(e, beam)) < 0.05)
          continue;
        if (!has_prediction(strong, i))
          continue;
        selected.push_back(strong_boxes[i]);
        beams.push_back(beam);
        for (const RangeSample &sample :
             range_samples(e, strong_boxes[i], cal_in.real(i, 2),
                           compute_zeta(e, beam))) {
          samples.push_back(sample);
        }
      }
      std::size_t used = 0;
      if (!(sigma_b > 0.0))
        sigma_b = beam_divergence(e, selected, beams, &used);
      if (!(sigma_m > 0.0)) {
        sigma_m =
            reflecting_range(samples, Scan::radians(e.scan.osc_width), 0.0);
      }
      source = "estimated from the indexed spots";
    }
    if ((!(sigma_b > 0.0) || !(sigma_m > 0.0)) && experiments.profile.present) {
      if (!(sigma_b > 0.0))
        sigma_b = experiments.profile.sigma_b;
      if (!(sigma_m > 0.0))
        sigma_m = experiments.profile.sigma_m;
      if (!(n_sigma > 0.0))
        n_sigma = experiments.profile.n_sigma;
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
    std::printf("Profile model: sigma_b %.4f deg, sigma_m %.4f deg (%s); "
                "foreground to %.1f sigma\n",
                sigma_b, sigma_m, source, mask_options.n_sigma);
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
    std::printf("Images: %s (from %s)\n", image_path.c_str(), image_source);

    std::unique_ptr<series::Series> images = series::nxmx(image_path);
    series::Info info;
    if (!images || !images->try_open(&info)) {
      std::fprintf(stderr,
                   "mxi_integrate: cannot open %s. If the data have moved "
                   "since the .expt was written, give --images\n",
                   image_path.c_str());
      return 1;
    }

    std::size_t workers =
        static_cast<std::size_t>(args.number("--threads", 0.0));
    if (workers == 0) {
      workers = std::thread::hardware_concurrency();
      if (workers == 0)
        workers = 1;
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
    std::printf("Predicted %zu reflections\n", predictions.size());

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
    {
      std::size_t left_out = outside_range;
      for (const auto &entry : refused)
        left_out += entry.second;
      std::printf("Integrating %zu of them", planned.size());
      if (left_out > 0) {
        std::printf("; not integrated:\n");
        if (outside_range > 0)
          std::printf("  %8zu  outside the image range\n", outside_range);
        for (const auto &entry : refused)
          std::printf("  %8zu  %s\n", entry.second, entry.first.c_str());
      } else {
        std::printf("\n");
      }
    }
    if (planned.empty())
      return 1;

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
    Column &miller =
        out.int_column("miller_index", "cctbx::miller::index<>", 3);
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
    // Where the reflection WAS against where it was predicted, in the frame
    // the images are in: pixels on the detector and images along the scan.
    // Both ends are in the same frame -- the prediction is the pixel that
    // fires, parallax included, and the centroid is of the counts that pixel
    // recorded -- so this is a straight difference and needs no correction.
    //
    // Not a DIALS column. For post-analysis of how well positions were
    // predicted, which is otherwise a join and a subtraction every time.
    Column &res_px = out.real_column("xyzres.px.value", "vec3<double>", 3);
    Column &res_px_var =
        out.real_column("xyzres.px.variance", "vec3<double>", 3);
    Column &iprf = out.real_column("intensity.prf.value", "double", 1);
    Column &iprf_var = out.real_column("intensity.prf.variance", "double", 1);
    Column &prf_cc = out.real_column("profile.correlation", "double", 1);
    // How much of each reflection the detector actually recorded. Not a DIALS
    // column; written because 11 per cent of reflections here touch a module
    // gap and 5 per cent lose half their box to one, and nothing else in the
    // table says so.
    Column &measured = out.real_column("profile.measured", "double", 1);
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
    if (save)
      saved.resize(planned.size());

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
      // DIALS' convention for a reflection that reaches a masked pixel, read
      // from its source and confirmed on a real integration. A foreground
      // with pixels missing makes the summed intensity not this reflection's
      // intensity -- summation cannot put back what the gap took -- so it is
      // not flagged as integrated by summation, and it says why. The profile
      // fitted intensity, which does put it back, keeps its own flag.
      //
      // Setting kIntegratedSum on these regardless let 1415 gap-crossing
      // reflections into dials.scale where DIALS lets 62, and the combined
      // intensity leans on the sum for strong reflections -- so a truncated
      // sum flagged as good is what scaling rejected, clustered along every
      // module edge.
      flags.ints[row] = flag::kPredicted | summation_flags(r);
      entering.ints[row] = p.s1.dot(axis.cross(s0)) > 0.0 ? 1 : 0;
      for (int k = 0; k < 6; ++k)
        bbox.ints[row * 6 + k] = box->bbox[k];
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
      for (int k = 0; k < 3; ++k)
        s1_column.reals[row * 3 + k] = p.s1[k];
      isum.reals[row] = r.intensity;
      ivar.reals[row] = r.variance;
      bmean.reals[row] = r.background_mean;
      bsum.reals[row] = r.background_sum;
      bsumvar.reals[row] = r.background_sum_variance;
      qe_column.reals[row] = quantum_efficiency(panel, p.s1);
      lp_column.reals[row] = lorentz_polarization(e.beam, e.goniometer, p.s1);
      d_column.reals[row] =
          e.crystal ? resolution(*e.crystal, p.h, p.k, p.l) : 0.0;
      const double z_of = compute_zeta(e, p.s1);
      zeta_column.reals[row] = z_of;
      part_column.reals[row] =
          partiality(e.scan, p.phi, z_of, mask_options.sigma_m, box->bbox[4],
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
      obs_mm_var.reals[row * 3 + 2] = r.centroid_variance_z *
                                      Scan::radians(e.scan.osc_width) *
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
      // The residual, NaN where there was no signal to find a centre in. The
      // observed column falls back to the prediction there because dials.scale
      // stops on a NaN in it; this column is ours, so it can say "nothing
      // measured" honestly instead of claiming a residual of exactly zero.
      //
      // The variance is the centroid's only. The prediction has an
      // uncertainty too, from the refined model's covariance, and it is not
      // propagated here -- so a pull of residual over sigma larger than one
      // is prediction error PLUS whatever that leaves out.
      {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        // From the UNCLIPPED centre of mass, whose uncertainty was checked;
        // not from xyzobs.px.value, whose is overconfident. So this is not
        // xyzobs - xyzcal, and is not meant to be.
        const bool has = r.unbiased_valid;
        res_px.reals[row * 3 + 0] = has ? r.unbiased_fast - p.px_fast : nan;
        res_px.reals[row * 3 + 1] = has ? r.unbiased_slow - p.px_slow : nan;
        res_px.reals[row * 3 + 2] = has ? r.unbiased_z - p.z : nan;
        res_px_var.reals[row * 3 + 0] = has ? r.unbiased_variance_fast : nan;
        res_px_var.reals[row * 3 + 1] = has ? r.unbiased_variance_slow : nan;
        res_px_var.reals[row * 3 + 2] = has ? r.unbiased_variance_z : nan;
      }
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
    // The scan is read in chunks of frames, and a box stays open from the
    // chunk its first frame is in to the chunk its last frame is in -- so a
    // frame is read once a pass however deep the boxes crossing it are.
    // Within a chunk the boxes are built, the frames fetched, decompressed
    // and filled, and the finished boxes integrated, each in parallel.
    //
    // The chunk is short because the boxes now outlive it: what is open at
    // once is the chunk's worth of boxes plus those still running on from
    // before, and a thousand frames of those is tens of gigabytes. Sixty-four
    // frames keeps sixteen threads decompressing four frames each.
    const std::size_t window =
        static_cast<std::size_t>(std::max(1.0, args.number("--window", 64.0)));
    // A cap on the boxes opened in one chunk, which shortens the chunk rather
    // than dropping boxes when it bites.
    const std::size_t max_boxes = static_cast<std::size_t>(
        std::max(1.0, args.number("--max-boxes", 20000.0)));

    //: Run `count` units of work over the pool, by index.
    // Run body(i, worker) for i in [0, count), where worker is this call's
    // number for the thread running it: 0 for the calling thread, 1 to n - 1
    // for the threads started here. It is stable within one call and means
    // nothing across calls, so anything kept per thread must be indexed by it
    // and not by thread_local state.
    //
    // thread_local lane numbers were a data race. The threads are new on every
    // call but the calling thread takes part in all of them, so its lane,
    // chosen once, outlived the call that chose it -- while each later call's
    // new threads were numbered from zero again. A probe caught two threads
    // sharing lane 0 in 32 calls of a thirty image run: the calling thread did
    // the whole of three small early calls on lane 0 and kept it. In profile
    // learning that meant two threads adding into the same partial profiles at
    // once. ThreadSanitizer saw nothing on that run only because the two rarely
    // both learned a reflection in the same call.
    const auto in_parallel_by_worker = [&](std::size_t count, auto &&body) {
      if (count == 0)
        return;
      const std::size_t n = std::min(workers, count);
      if (n <= 1) {
        for (std::size_t i = 0; i < count; ++i)
          body(i, std::size_t{0});
        return;
      }
      std::atomic<std::size_t> next_unit{0};
      std::vector<std::thread> pool;
      pool.reserve(n - 1);
      const auto run = [&](std::size_t worker) {
        for (;;) {
          const std::size_t i = next_unit.fetch_add(1);
          if (i >= count)
            break;
          body(i, worker);
        }
      };
      for (std::size_t t = 1; t < n; ++t)
        pool.emplace_back(run, t);
      run(0);
      for (std::thread &t : pool)
        t.join();
    };
    const auto in_parallel = [&](std::size_t count, auto &&body) {
      in_parallel_by_worker(count,
                            [&](std::size_t i, std::size_t) { body(i); });
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
    grid_spec.subdivisions =
        static_cast<int>(args.number("--subdivisions", 5.0));
    ReferenceProfiles reference = make_reference(
        grid_spec, static_cast<int>(args.number("--regions", 3.0)),
        static_cast<int>(args.number("--scan-blocks", 5.0)), e.detector.size(),
        first_image, last_image);
    const bool fitting = !args.has("--summation-only");
    // Which reflections are worth learning from: strong, nearly whole, and
    // mostly inside the grid. DIALS marks these `reference_spot`.
    const double least_signal = args.number("--reference-signal", 10.0);
    // How much of a reflection must have been measured for its fit to count.
    const double least_measured = args.number("--least-measured", 0.6);

    std::vector<Shoebox> boxes;
    std::size_t at = 0;
    int pass = 0;
    std::size_t references_used = 0;
    double t_transform = 0.0, t_fit = 0.0;
    // Boxes OUTLIVE a chunk of frames. Each chunk opens the boxes that start
    // in it, reads its frames once into every open box, and closes the boxes
    // whose last frame it held.
    //
    // It used to be windows of boxes that read every frame they touched and
    // then discarded the boxes. A box starting in one window and ending past
    // it made that window read on into the next one's frames, and the next
    // window read them again for its own -- and near-axis reflections run forty
    // to eighty frames deep, so every boundary repeated them. On a 3600 frame
    // Eiger 16M run that was 20224 frames read for 3600 wanted, 5.62 reads a
    // frame over two passes where 2.0 is the floor. Decompression was the
    // same 10 to 12 ms a frame as the spot finder; there were five times as
    // many of them.
    std::vector<Shoebox> active;
    std::vector<std::size_t> active_rows;
    std::int32_t chunk_start = planned.empty() ? 0 : planned[0].bbox[4];
    while (at < planned.size() || !active.empty()) {
      if (active.empty() && at < planned.size()) {
        chunk_start = std::max(chunk_start, planned[at].bbox[4]);
      }
      const std::int32_t chunk_limit =
          chunk_start + static_cast<std::int32_t>(window);
      // Open everything starting before the chunk's end -- capped by
      // --max-boxes, but never part way through the boxes of one frame, and
      // when the cap bites the chunk ENDS where opening stopped. Otherwise a
      // box not yet opened would miss the frames the chunk reads.
      std::size_t stop = at;
      while (stop < planned.size() && planned[stop].bbox[4] < chunk_limit) {
        if (stop - at >= max_boxes && planned[stop].bbox[4] > chunk_start)
          break;
        ++stop;
      }
      const std::int32_t chunk_end =
          (stop < planned.size() && planned[stop].bbox[4] < chunk_limit)
              ? planned[stop].bbox[4]
              : chunk_limit;
      const std::size_t opening = stop - at;

      // The new boxes, built in parallel: this is the mask, and it was 46
      // seconds before it was.
      const double t_open_start = now_wall();
      const std::size_t first_new = active.size();
      active.resize(first_new + opening);
      active_rows.resize(first_new + opening);
      in_parallel(opening, [&](std::size_t i) {
        active_rows[first_new + i] = at + i;
        BoxRejection ignored = BoxRejection::kNone;
        Shoebox &box = active[first_new + i];
        if (build_shoebox(e, *planned[at + i].prediction, mask_options, &box,
                          &ignored)) {
          box.data.assign(box.size(), 0.0f);
        }
      });
      at = stop;
      t_open += now_wall() - t_open_start;
      most_open = std::max(most_open, active.size());

      // Which open boxes each of THIS chunk's frames touches.
      std::map<std::int32_t, std::vector<std::size_t>> touching;
      for (std::size_t i = 0; i < active.size(); ++i) {
        const std::int32_t lo = std::max(active[i].bbox[4], chunk_start);
        const std::int32_t hi = std::min(active[i].bbox[5], chunk_end);
        for (std::int32_t z = lo; z < hi; ++z)
          touching[z].push_back(i);
      }
      std::vector<std::int32_t> frame_numbers;
      frame_numbers.reserve(touching.size());
      for (const auto &entry : touching)
        frame_numbers.push_back(entry.first);

      // Fetch, decompress and fill, in parallel over frames. Each frame writes
      // only its own z plane of each box it touches, so two frames of one box
      // never touch the same voxel.
      std::atomic<std::size_t> frames_done{0};
      std::atomic<std::size_t> bad_here{0};
      std::atomic<std::size_t> unread{0};
      std::vector<double> fetch_by_thread(std::max<std::size_t>(workers, 1),
                                          0.0);
      std::vector<double> decompress_by_thread(fetch_by_thread.size(), 0.0);
      std::vector<double> fill_by_thread(fetch_by_thread.size(), 0.0);
      const double t_region_start = now_wall();
      in_parallel_by_worker(
          frame_numbers.size(), [&](std::size_t which, std::size_t worker) {
            // A reader per thread, kept for the thread's life: it belongs to
            // one thread only, so outliving a call is harmless. Its TIMING lane
            // is the worker number, which is unique within the call; a
            // thread_local lane was shared between the calling thread and a new
            // one.
            thread_local std::unique_ptr<series::Reader> mine;
            if (!mine)
              mine = images->reader();
            const std::size_t slot = worker % fetch_by_thread.size();
            const std::int32_t z = frame_numbers[which];
            if (z < 0 || static_cast<std::size_t>(z) >= keys.size())
              return;
            series::Frame raw;
            const double t0 = now_wall();
            if (!mine->read(keys[static_cast<std::size_t>(z)], &raw)) {
              // A frame the writer never received: an unallocated chunk.
              // Counted, because a shoebox that spans it is missing a slice and
              // will integrate low, and silently dropping it leaves "frames
              // read" less than the number of images with no explanation.
              unread.fetch_add(1);
              return;
            }
            const double t1 = now_wall();
            const std::size_t height = static_cast<std::size_t>(raw.height);
            const std::size_t width = static_cast<std::size_t>(raw.width);
            const std::size_t bytes =
                decompress::frame_bytes(height, width, raw.bit_depth);
            std::vector<std::uint8_t> pixels(bytes);
            decompress::image(raw.data, raw.algorithm, raw.bit_depth, height,
                              width, {pixels.data(), bytes});
            const double t2 = now_wall();
            frames_done.fetch_add(1);

            const auto fill = [&](auto typed, std::uint32_t bad) {
              using Pixel = decltype(typed);
              const Pixel *const raw_pixels =
                  reinterpret_cast<const Pixel *>(pixels.data());
              std::size_t bad_count = 0;
              for (std::size_t i : touching[z]) {
                Shoebox &box = active[i];
                if (box.data.empty())
                  continue;
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
                    // The largest representable value is the bad-pixel marker,
                    // not a count: excluded from both sums rather than counted
                    // as zero, which would drag the background down wherever a
                    // module gap crosses a shoebox.
                    if (v == bad) {
                      // Clear VALIDITY and keep the region. A voxel in a module
                      // gap is still a foreground voxel, it just has no
                      // measurement in it -- and the profile fit needs to know
                      // that a part of the reflection is missing rather than
                      // simply not seeing it. Zeroing the whole mask made a
                      // reflection with a third of its foreground in a gap
                      // report two thirds of its intensity, which is exactly
                      // what profile fitting exists to avoid.
                      box.mask[into] &=
                          static_cast<std::uint8_t>(~shoebox_mask::kValid);
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
      // Only frames the series has. A box can reach past the last image, and
      // counting frames that cannot be read made a slice of 300 frames read
      // twice report "356 wanted, each read 1.69 times" instead of 2.00.
      for (const auto &entry : touching) {
        if (entry.first >= 0 &&
            static_cast<std::size_t>(entry.first) < keys.size()) {
          wanted.insert(entry.first);
        }
      }
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
        for (double x : v)
          all += x;
        return all;
      };
      t_fetch += summed(fetch_by_thread);
      t_decompress += summed(decompress_by_thread);
      t_fill += summed(fill_by_thread);

      // Close what is finished: every frame a box spans is now filled.
      std::vector<Shoebox> boxes;
      std::vector<std::size_t> rows;
      {
        std::vector<Shoebox> staying;
        std::vector<std::size_t> staying_rows;
        for (std::size_t i = 0; i < active.size(); ++i) {
          if (active[i].bbox[5] <= chunk_end) {
            boxes.push_back(std::move(active[i]));
            rows.push_back(active_rows[i]);
          } else {
            staying.push_back(std::move(active[i]));
            staying_rows.push_back(active_rows[i]);
          }
        }
        active.swap(staying);
        active_rows.swap(staying_rows);
      }
      const std::size_t count = boxes.size();
      chunk_start = chunk_end;

      if (pass == 0) {
        // Integrate, in parallel over boxes.
        const double t_close_start = now_wall();
        in_parallel(count, [&](std::size_t i) {
          if (boxes[i].data.empty())
            return;
          close(rows[i], &boxes[i]);
        });
        t_integrate += now_wall() - t_close_start;

        // Learn from the ones worth learning from, in parallel.
        //
        // The transform is the cost and it is per reflection with nothing
        // shared; only the accumulation into a cell's profile is shared, and
        // that is a few hundred adds against a few hundred thousand
        // multiplies. So each thread keeps its own set of profiles and they
        // are added together at the end -- which also makes the result
        // independent of the thread count, since addition of the same numbers
        // in a different order is the only thing that changes.
        if (fitting) {
          const double t0 = now_wall();
          // FIXED BLOCKS, NOT THREADS. The closing reflections are cut into a
          // fixed number of blocks by index, each summed in index order by one
          // task into its own partial, and the partials added in block order.
          // Which thread runs a block is left to the scheduler and changes
          // nothing: every addition happens in the same order on every run,
          // whatever the thread count. Refinement's reduction already worked
          // this way.
          //
          // It replaced partials per thread, which were a data race (see
          // in_parallel_by_worker) and, even without the race, added in an
          // order the scheduler chose: two four-thread runs differed by 3.6e-11
          // in intensity.prf.value.
          //
          // Sixteen blocks bounds the parallelism of this phase at sixteen,
          // which is a few per cent of a run, and costs sixteen partial sets of
          // profiles: about two megabytes each at 324 cells.
          constexpr std::size_t kBlocks = 16;
          const std::size_t blocks =
              std::min(kBlocks, std::max<std::size_t>(count, 1));
          std::vector<ReferenceProfiles> partial(blocks, reference);
          for (ReferenceProfiles &r : partial) {
            for (std::vector<double> &p : r.profile)
              std::fill(p.begin(), p.end(), 0.0);
            std::fill(r.spots.begin(), r.spots.end(), 0);
          }
          std::atomic<std::size_t> learned{0};
          in_parallel(blocks, [&](std::size_t b) {
            const std::size_t first = b * count / blocks;
            const std::size_t last = (b + 1) * count / blocks;
            for (std::size_t i = first; i < last; ++i) {
              if (boxes[i].data.empty())
                continue;
              const std::size_t row = rows[i];
              const double signal = isum.reals[row];
              const double sigma = std::sqrt(std::max(ivar.reals[row], 1e-12));
              if (!(signal > least_signal * sigma))
                continue;
              if (part_column.reals[row] < 0.99)
                continue;
              const Prediction &q = *planned[row].prediction;
              const Transformed t =
                  transform_shoebox(e, boxes[i], q.s1, q.phi, grid_spec);
              if (!t.valid || t.outside > 0.05)
                continue;
              const std::size_t region =
                  reference.region_of(panel, static_cast<std::size_t>(q.panel),
                                      q.px_fast, q.px_slow, q.z);
              if (add_reference(&partial[b], region, t))
                learned.fetch_add(1);
            }
          });
          for (const ReferenceProfiles &r : partial) {
            for (std::size_t g = 0; g < reference.profile.size(); ++g) {
              for (std::size_t k = 0; k < reference.profile[g].size(); ++k) {
                reference.profile[g][k] += r.profile[g][k];
              }
              reference.spots[g] += r.spots[g];
            }
          }
          references_used += learned.load();
          t_transform += now_wall() - t0;
        }
      } else {
        // Fit, in parallel over boxes: each writes only its own row.
        const double t0 = now_wall();
        in_parallel(count, [&](std::size_t i) {
          if (boxes[i].data.empty())
            return;
          const std::size_t row = rows[i];
          const Prediction &q = *planned[row].prediction;
          // The background the GLM found in the first pass, put back so the
          // fit subtracts the same thing the sum did.
          boxes[i].background.assign(boxes[i].size(),
                                     static_cast<float>(bmean.reals[row]));
          // No transform here. The second pass fits against the pixels, so
          // carrying the counts onto the grid is work whose answer is thrown
          // away -- and it is the same cost as carrying the profile back, so
          // doing both doubled this phase.
          //
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
          const std::vector<double> on_pixels =
              profile_on_pixels(e, boxes[i], q.s1, q.phi, grid_spec, local);
          if (on_pixels.empty())
            return;
          const ProfileFit fit =
              fit_on_pixels(boxes[i], on_pixels, integrate_options.gain);
          if (!fit.valid)
            return;
          // A fit is an extrapolation when part of the reflection is missing,
          // and past some point it is guesswork dressed as a measurement. The
          // intensity is still written, so it can be looked at; the flag that
          // says it was profile fitted is not, so nothing downstream merges it
          // by accident.
          measured.reals[row] = fit.measured;
          if (fit.measured < least_measured)
            return;
          iprf.reals[row] = fit.intensity;
          iprf_var.reals[row] = fit.variance;
          prf_cc.reals[row] = fit.correlation;
          flags.ints[row] |= flag::kIntegratedPrf;
        });
        t_fit += now_wall() - t0;
      }

      if (save) {
        for (std::size_t i = 0; i < count; ++i) {
          // Into DIALS' convention before it goes to the file.
          //
          // Here a bad pixel keeps its region flag and loses only Valid, so
          // the profile fit can tell a foreground voxel with nothing in it
          // from one that was never foreground -- which is what lets a
          // reflection crossing a module gap keep its intensity. DIALS' own
          // mask calculator sets Foreground and Background only on voxels that
          // are already Valid, so a region bit without Valid never occurs
          // there, and a table carrying them was rejected as an invalid
          // structure. The file is DIALS' format and follows DIALS'
          // convention: a voxel with no measurement is zero, as the spot
          // finder has always written it.
          to_dials_convention(&boxes[i]);
          saved[rows[i]] = std::move(boxes[i]);
        }
      }
      if (at >= planned.size() && active.empty() && pass == 0 && fitting) {
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
              for (double v : reference.profile[r])
                std::fprintf(f, "%.9g\n", v);
            }
            std::fclose(f);
            std::printf("wrote %s\n", profile_path.c_str());
          }
        }
        std::size_t empty = 0;
        for (std::size_t r = 0; r < reference.spots.size(); ++r) {
          if (reference.spots[r] < 10)
            ++empty;
        }
        std::printf(
            "Reference profiles: %zu learned from %zu reflections (%zu of %zu "
            "regions borrowed the detector average)\n",
            reference.region_count() - empty, references_used, empty,
            reference.region_count());
        if (references_used == 0) {
          std::fprintf(stderr,
                       "mxi_integrate: nothing to learn profiles from, so no "
                       "profile fitting: lower --reference-signal or check the "
                       "summation\n");
        } else {
          pass = 1;
          at = 0;
          chunk_start = planned.empty() ? 0 : planned[0].bbox[4];
          continue;
        }
      }
    }
    if (save)
      shoebox_bytes = encode_shoeboxes(saved);
    if (frames_missing > 0) {
      std::fprintf(
          stderr,
          "mxi_integrate: %zu reads found no frame: an unallocated chunk is a "
          "frame the writer never received, and a shoebox spanning one is "
          "missing a slice\n",
          frames_missing);
    }

    // What a user reads to judge the run: the summaries dials.integrate
    // prints, against resolution and overall, from the columns just written.
    {
      SummaryInput summary;
      const std::size_t n = planned.size();
      summary.d.assign(d_column.reals.begin(), d_column.reals.begin() + n);
      summary.flags.assign(flags.ints.begin(), flags.ints.begin() + n);
      summary.intensity_sum.assign(isum.reals.begin(), isum.reals.begin() + n);
      summary.variance_sum.assign(ivar.reals.begin(), ivar.reals.begin() + n);
      summary.intensity_prf.assign(iprf.reals.begin(), iprf.reals.begin() + n);
      summary.variance_prf.assign(iprf_var.reals.begin(),
                                  iprf_var.reals.begin() + n);
      summary.profile_correlation.assign(prf_cc.reals.begin(),
                                         prf_cc.reals.begin() + n);
      summary.background.assign(bmean.reals.begin(), bmean.reals.begin() + n);
      summary.partiality.assign(part_column.reals.begin(),
                                part_column.reals.begin() + n);
      summary.res_fast.resize(n);
      summary.res_slow.resize(n);
      // The OBSERVED centre against the prediction, over reflections that were
      // detected -- which is what dials.integrate's RMSD XY is, so the two can
      // be read side by side. Not xyzres.px: that is the unclipped centre,
      // honest about its noise and so twice as scattered for weak spots
      // (0.544 px against 0.267 on a 300 image insulin sweep, where DIALS gave
      // about 0.25). It is the right thing for judging predicted positions,
      // and the wrong one for a column people will set beside DIALS'.
      // xyzres.px being NaN is what says a reflection was not detected, since
      // xyzobs.px falls back to the prediction and would read as perfect.
      const double nan = std::numeric_limits<double>::quiet_NaN();
      for (std::size_t i = 0; i < n; ++i) {
        // Over reflections integrated by summation, as DIALS reports it: one
        // crossing a module gap has a centre pulled off by the gap, which
        // took the overall figure from 0.267 px to 0.343.
        const bool detected = std::isfinite(res_px.reals[i * 3 + 0]) &&
                              (flags.ints[i] & flag::kIntegratedSum) != 0;
        summary.res_fast[i] =
            detected ? obs_px.reals[i * 3 + 0] - cal_px.reals[i * 3 + 0] : nan;
        summary.res_slow[i] =
            detected ? obs_px.reals[i * 3 + 1] - cal_px.reals[i * 3 + 1] : nan;
      }
      print_summary(stdout, summarise(summary, 10));
    }
    std::printf("\n");

    if (save) {
      Table::Opaque column;
      column.type = "Shoebox<>";
      column.rows = planned.size();
      std::printf("Shoeboxes kept: %.1f MB\n", shoebox_bytes.size() / 1e6);
      column.bytes = std::move(shoebox_bytes);
      out.set_opaque("shoebox", std::move(column));
    }
    if (!e.identifier.empty())
      out.identifiers[0] = e.identifier;

    const std::string path = args.value("-o", "integrated.refl");
    const double t_write_start = now_wall();
    write_reflections(path, out);
    t_write = now_wall() - t_write_start;
    std::printf("Wrote %zu reflections to %s\n", planned.size(), path.c_str());

    if (args.has("--timing")) {
      const double total = now_wall() - t_start;
      const auto line = [&](const char *name, double seconds) {
        std::printf("  %-26s %8.3f s  %5.1f%%\n", name, seconds,
                    total > 0.0 ? 100.0 * seconds / total : 0.0);
      };
      std::printf("\nTiming: %zu threads, chunks of %zu frames, at most %zu "
                  "boxes opened a chunk, %zu open at once\n",
                  workers, window, max_boxes, most_open);
      std::printf("  %zu frames read, %zu wanted by a shoebox, each read %.2f "
                  "times; %zu voxels on masked pixels\n",
                  frames_read, wanted.size(),
                  wanted.empty() ? 0.0
                                 : static_cast<double>(frames_read) /
                                       static_cast<double>(wanted.size()),
                  bad_pixels);
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

} // namespace mxi

int main(int argc, char **argv) {
  // Mirrored to mxi_integrate.log in the working directory, as DIALS writes
  // dials.<program>.log; not for a run that only asks for help.
  if (!mxi::only_asks_for_help(argc, argv))
    mxi::mirror_to_log("mxi_integrate.log");
  return mxi::run_program(argc, argv);
}
