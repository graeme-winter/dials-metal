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
      "  --block-size N    images held in memory at once (50). Every shoebox\n"
      "                    at once is 233 GB on ten rotations of insulin\n"
      "  --save-shoeboxes  keep the pixels and the mask in the output\n",
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
      "--gain",      "--save-shoeboxes", "--images", "--block-size"};
  std::set<std::string> takes_value = known;
  takes_value.erase("--save-shoeboxes");
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
    const ExperimentList experiments = read_experiments(args.positional[0]);
    if (experiments.size() == 0) {
      std::fprintf(stderr, "mxi_integrate: no experiments\n");
      return 1;
    }
    const Experiment &e = experiments[0];
    const Panel &panel = e.detector[0];

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
    const std::vector<Prediction> predictions = predict(e, predict_options);
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

    const std::size_t block_size =
        static_cast<std::size_t>(std::max(1.0, args.number("--block-size", 50.0)));

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

    const bool save = args.has("--save-shoeboxes");
    std::string shoebox_bytes;

    std::unique_ptr<series::Reader> reader = images->reader();
    const std::vector<std::string> keys = images->ready();
    // Frame number to key, so a block can ask for the frames it needs rather
    // than walking every key it does not.
    std::map<std::int64_t, std::string> key_of;
    for (const std::string &key : keys) {
      series::Frame probe;
      if (reader->read(key, &probe)) key_of[probe.number] = key;
    }

    series::Frame frame;
    std::vector<std::uint8_t> pixels;
    std::size_t frames_read = 0;
    std::size_t bad_pixels = 0;
    std::size_t integrated = 0;
    std::size_t at = 0;
    const Vec3 s0 = e.beam.s0();
    const Vec3 axis = e.goniometer.lab_axis();

    while (at < planned.size()) {
      const std::int32_t block_start = planned[at].bbox[4];
      const std::int32_t block_end =
          block_start + static_cast<std::int32_t>(block_size);
      std::size_t stop = at;
      std::int32_t highest = block_start;
      while (stop < planned.size() && planned[stop].bbox[4] < block_end) {
        highest = std::max(highest, planned[stop].bbox[5]);
        ++stop;
      }

      std::vector<Shoebox> boxes;
      boxes.reserve(stop - at);
      for (std::size_t i = at; i < stop; ++i) {
        Shoebox box;
        if (!build_shoebox(e, *planned[i].prediction, mask_options, &box)) {
          // It had a box a moment ago, so this cannot happen; if it ever does,
          // an empty box integrates to nothing rather than shifting every row
          // after it.
          box.panel = 0;
          for (int k = 0; k < 6; ++k) box.bbox[k] = planned[i].bbox[k];
        }
        box.data.assign(box.size(), 0.0f);
        boxes.push_back(std::move(box));
      }

      for (std::int32_t z = block_start; z < highest; ++z) {
        const auto found = key_of.find(z);
        if (found == key_of.end()) continue;
        if (!reader->read(found->second, &frame)) continue;
        const std::size_t height = static_cast<std::size_t>(frame.height);
        const std::size_t width = static_cast<std::size_t>(frame.width);
        const std::size_t bytes =
            decompress::frame_bytes(height, width, frame.bit_depth);
        pixels.resize(bytes);
        decompress::image(frame.data, frame.algorithm, frame.bit_depth, height,
                          width, {pixels.data(), bytes});
        ++frames_read;

        const auto fill = [&](auto typed, std::uint32_t bad) {
          using Pixel = decltype(typed);
          const Pixel *const raw =
              reinterpret_cast<const Pixel *>(pixels.data());
          for (std::size_t i = 0; i < boxes.size(); ++i) {
            Shoebox &box = boxes[i];
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
                // The largest representable value is the detector's bad-pixel
                // marker, not a count: 5.8 per cent of an Eiger frame, module
                // gaps and dead pixels. Excluded from both sums rather than
                // counted as zero, which would drag the background down
                // wherever a gap crosses a shoebox.
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
      }

      for (std::size_t i = at; i < stop; ++i) {
        Shoebox &box = boxes[i - at];
        const Prediction &p = *planned[i].prediction;
        const IntegratedReflection r =
            integrate_shoebox(&box, integrate_options);
        miller.ints[i * 3 + 0] = p.h;
        miller.ints[i * 3 + 1] = p.k;
        miller.ints[i * 3 + 2] = p.l;
        panel_column.ints[i] = static_cast<std::int64_t>(p.panel);
        id.ints[i] = 0;
        imageset.ints[i] = 0;
        flags.ints[i] =
            flag::kPredicted | (r.valid ? flag::kIntegratedSum : 0);
        entering.ints[i] = p.s1.dot(axis.cross(s0)) > 0.0 ? 1 : 0;
        for (int k = 0; k < 6; ++k) bbox.ints[i * 6 + k] = box.bbox[k];
        n_fg.ints[i] = static_cast<std::int64_t>(r.n_foreground);
        n_bg.ints[i] = static_cast<std::int64_t>(r.n_background);
        n_val.ints[i] = static_cast<std::int64_t>(r.n_valid);
        cal_px.reals[i * 3 + 0] = p.px_fast;
        cal_px.reals[i * 3 + 1] = p.px_slow;
        cal_px.reals[i * 3 + 2] = p.z;
        const auto mm = panel.px_to_mm(p.px_fast, p.px_slow);
        cal_mm.reals[i * 3 + 0] = mm.first;
        cal_mm.reals[i * 3 + 1] = mm.second;
        cal_mm.reals[i * 3 + 2] = p.phi;
        for (int k = 0; k < 3; ++k) s1_column.reals[i * 3 + k] = p.s1[k];
        isum.reals[i] = r.intensity;
        ivar.reals[i] = r.variance;
        bmean.reals[i] = r.background_mean;
        bsum.reals[i] = r.background_sum;
        bsumvar.reals[i] = r.background_sum_variance;
        qe_column.reals[i] = quantum_efficiency(panel, p.s1);
        lp_column.reals[i] = lorentz_polarization(e.beam, e.goniometer, p.s1);
        // Where the spot was, as against where the model put it. Falls back to
        // the prediction when there was no signal to find a centre of mass in,
        // because a NaN here stops dials.scale rather than being ignored by it.
        const double of = r.centroid_valid ? r.centroid_fast : p.px_fast;
        const double os = r.centroid_valid ? r.centroid_slow : p.px_slow;
        const double oz = r.centroid_valid ? r.centroid_z : p.z;
        obs_px.reals[i * 3 + 0] = of;
        obs_px.reals[i * 3 + 1] = os;
        obs_px.reals[i * 3 + 2] = oz;
        obs_px_var.reals[i * 3 + 0] = r.centroid_variance_fast;
        obs_px_var.reals[i * 3 + 1] = r.centroid_variance_slow;
        obs_px_var.reals[i * 3 + 2] = r.centroid_variance_z;
        const auto obs_mm_pair = panel.px_to_mm(of, os);
        obs_mm.reals[i * 3 + 0] = obs_mm_pair.first;
        obs_mm.reals[i * 3 + 1] = obs_mm_pair.second;
        obs_mm.reals[i * 3 + 2] = e.scan.phi_from_z(oz);
        if (r.valid) ++integrated;
      }
      if (save) shoebox_bytes += encode_shoeboxes(boxes);
      at = stop;
    }
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
    write_reflections(path, out);
    std::printf("wrote %s\n", path.c_str());
    return 0;
  } catch (const std::exception &error) {
    std::fprintf(stderr, "mxi_integrate: %s\n", error.what());
    return 1;
  }
}
