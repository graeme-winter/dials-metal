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
      "  --d-min D         resolution limit\n"
      "  --first-image N --last-image N   restrict to part of the scan\n"
      "  --gain G          detector gain, counts per photon (1)\n"

      "  --save-shoeboxes  keep the pixels and the mask in the output\n"
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
      "--gain",      "--save-shoeboxes", "--images", "--timing"};
  std::set<std::string> takes_value = known;
  takes_value.erase("--save-shoeboxes");
  takes_value.erase("--timing");
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

    PredictOptions predict_options;
    predict_options.d_min = args.number("--d-min", 0.0);
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
    std::size_t bad_pixels = 0;
    std::size_t integrated = 0;
    std::size_t most_open = 0;
    const Vec3 s0 = e.beam.s0();
    const Vec3 axis = e.goniometer.lab_axis();

    // The shoeboxes currently open, by their row in the table.
    std::map<std::size_t, Shoebox> open;
    std::size_t next = 0;

    const auto close = [&](std::size_t row, Shoebox *box) {
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
      if (r.valid) ++integrated;
      if (save) saved[row] = std::move(*box);
    };

    std::int64_t previous = -1;
    for (const std::string &key : keys) {
      const double t_frame_start = now_wall();
      if (!reader->read(key, &frame)) continue;
      t_fetch += now_wall() - t_frame_start;
      const std::int32_t z = static_cast<std::int32_t>(frame.number);
      if (frame.number <= previous) {
        throw std::runtime_error(
            "the image series returned frame " + std::to_string(frame.number) +
            " after " + std::to_string(previous) +
            "; this reads them in order and cannot take them out of it");
      }
      previous = frame.number;

      // Open everything that starts at or before this frame.
      const double t_open_start = now_wall();
      while (next < planned.size() && planned[next].bbox[4] <= z) {
        Shoebox box;
        BoxRejection ignored = BoxRejection::kNone;
        if (build_shoebox(e, *planned[next].prediction, mask_options, &box,
                          &ignored)) {
          box.data.assign(box.size(), 0.0f);
          open.emplace(next, std::move(box));
        }
        ++next;
      }
      t_open += now_wall() - t_open_start;
      most_open = std::max(most_open, open.size());

      {
        const std::size_t height = static_cast<std::size_t>(frame.height);
        const std::size_t width = static_cast<std::size_t>(frame.width);
        const std::size_t bytes =
            decompress::frame_bytes(height, width, frame.bit_depth);
        pixels.resize(bytes);
        const double t_decompress_start = now_wall();
        decompress::image(frame.data, frame.algorithm, frame.bit_depth, height,
                          width, {pixels.data(), bytes});
        ++frames_read;
        t_decompress += now_wall() - t_decompress_start;
        const double t_fill_start = now_wall();

        const auto fill = [&](auto typed, std::uint32_t bad) {
          using Pixel = decltype(typed);
          const Pixel *const raw = reinterpret_cast<const Pixel *>(pixels.data());
          for (auto &entry : open) {
            Shoebox &box = entry.second;
            if (z < box.bbox[4] || z >= box.bbox[5]) continue;
            const std::int32_t zi = z - box.bbox[4];
            for (std::int32_t y = 0; y < box.ny(); ++y) {
              const std::size_t row =
                  static_cast<std::size_t>(box.bbox[2] + y) * width;
              for (std::int32_t x = 0; x < box.nx(); ++x) {
                const std::size_t index =
                    row + static_cast<std::size_t>(box.bbox[0] + x);
                const std::size_t into = box.at(x, y, zi);
                const std::uint32_t v = static_cast<std::uint32_t>(raw[index]);
                // The largest representable value is the bad-pixel marker, not
                // a count: excluded from both sums rather than counted as zero,
                // which would drag the background down wherever a module gap
                // crosses a shoebox.
                if (v == bad) {
                  box.mask[into] = 0;
                  ++bad_pixels;
                } else {
                  box.data[into] = static_cast<float>(v);
                }
              }
            }
          }
        };
        if (frame.bit_depth == 16) {
          fill(std::uint16_t{}, 0xFFFFu);
        } else if (frame.bit_depth == 32) {
          fill(std::uint32_t{}, 0xFFFFFFFFu);
        } else {
          throw std::runtime_error("unsupported bit depth " +
                                   std::to_string(frame.bit_depth));
        }
        t_fill += now_wall() - t_fill_start;
      }

      // Close everything that ends here.
      const double t_close_start = now_wall();
      for (auto it = open.begin(); it != open.end();) {
        if (it->second.bbox[5] <= z + 1) {
          close(it->first, &it->second);
          it = open.erase(it);
        } else {
          ++it;
        }
      }
      t_integrate += now_wall() - t_close_start;
    }
    // Anything still open ran past the end of the images.
    for (auto &entry : open) close(entry.first, &entry.second);
    open.clear();
    if (save) shoebox_bytes = encode_shoeboxes(saved);
    std::printf("at most %zu shoeboxes open at once\n", most_open);
    std::printf("%zu frames read, %zu bad pixels masked\n", frames_read,
                bad_pixels);
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
      line("fetching frames", t_fetch);
      line("decompressing", t_decompress);
      line("filling shoeboxes", t_fill);
      line("background and summation", t_integrate);
      line("writing", t_write);
      std::printf("  %-26s %8.3f s\n", "total", total);
    }
    return 0;
  } catch (const std::exception &error) {
    std::fprintf(stderr, "mxi_integrate: %s\n", error.what());
    return 1;
  }
}
