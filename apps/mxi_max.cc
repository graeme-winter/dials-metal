// mxi_max: the largest count on any valid pixel of a series, and whether it
// would fit in 16 bits.
//
// Metal's threshold takes 16-bit frames only, Apple's GPUs having no double
// precision, and many detectors write 32-bit frames whose counts, the dose
// kept down for radiation damage, never come near 65535. If no valid pixel of
// a series reaches 0xFFFD -- the values above it the 16-bit markers -- its
// frames can be narrowed to 16 bits and thresholded on the GPU,
// mxi_find --gpu --gpu-force, losing nothing. This says whether that is so: it
// reads and decompresses every frame on every core, as mxi_readtest does, and
// takes the largest value that is not a marker: max() - 1 a bad pixel, max() a
// tile join, in the pixel's width.

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <mutex>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "args.hh"
#include "decompress.hh"
#include "expt.hh"
#include "series.hh"

namespace {

double now() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

//: The first value 16-bit data reserves: 0xFFFD and above are markers.
constexpr std::uint64_t kSixteenBitLimit = 0xFFFD;

//: The largest value in a frame that is not a marker, and how many pixels are
//: each marker: max() - 1 a bad pixel, max() a tile join -- 0xfffe and 0xffff
//: in 16 bits, 0xfffffffe and 0xffffffff in 32 -- neither a count. One pass, a
//: loop with no branch the compiler can make SIMD.
template <typename T>
void scan(const std::uint8_t *bytes, std::size_t n, std::uint64_t *largest,
          std::uint64_t *bad, std::uint64_t *joins, std::uint64_t *over) {
  const T *pixels = reinterpret_cast<const T *>(bytes);
  const T join = static_cast<T>(~T{0}), dead = static_cast<T>(join - 1);
  T best = 0;
  std::uint64_t dead_n = 0, join_n = 0, above = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const T v = pixels[i];
    dead_n += v == dead;
    join_n += v == join;
    const T real = v >= dead ? T{0} : v;
    above += static_cast<std::uint64_t>(real) >= kSixteenBitLimit;
    best = std::max(best, real);
  }
  *largest = best;
  *bad = dead_n;
  *joins = join_n;
  *over = above;
}

bool looks_like_json(const std::string &path) {
  std::FILE *file = std::fopen(path.c_str(), "rb");
  if (file == nullptr)
    return false;
  int c = 0;
  while ((c = std::fgetc(file)) != EOF && std::isspace(c))
    ;
  std::fclose(file);
  return c == '{';
}

const char *kUsage =
    "usage: mxi_max EXPT|MASTER [-j N]\n"
    "  The largest count on any valid pixel of the images -- those an "
    "experiment\n"
    "  list names, or a master given itself -- and whether it fits in 16 "
    "bits:\n"
    "  below 0xFFFD, the values above being the 16-bit markers. If it does,\n"
    "  mxi_find --gpu --gpu-force can threshold 32-bit frames on a GPU that "
    "takes\n"
    "  16 bits only. Exits 0 if it fits, 1 if not.\n"
    "  -j, --threads N   N threads (every core)\n";

} // namespace

int main(int argc, char **argv) {
  using namespace mxi;
  const std::set<std::string> known = {"-j", "--threads", "--help"};
  const std::set<std::string> takes_value = {"-j", "--threads"};
  const Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (args.has("--help")) {
    std::printf("%s", kUsage);
    return 0;
  }
  if (!args.ok || args.positional.size() != 1) {
    if (!args.ok)
      std::fprintf(stderr, "mxi_max: %s\n", args.error.c_str());
    std::fprintf(stderr, "%s", kUsage);
    return 2;
  }
  try {
    std::string master = args.positional[0];
    if (looks_like_json(master)) {
      master = read_experiments(master).image_template;
      if (master.empty() || master.find('#') != std::string::npos)
        throw std::runtime_error(args.positional[0] +
                                 " names no single NXmx master");
    }
    std::size_t threads = static_cast<std::size_t>(
        args.number("--threads", args.number("-j", 0.0)));
    if (threads == 0)
      threads = std::max(1u, std::thread::hardware_concurrency());

    std::unique_ptr<series::Series> s = series::nxmx(master);
    series::Info info;
    if (!s->try_open(&info))
      throw std::runtime_error("cannot open " + master + " as an NXmx series");
    const std::uint64_t frames = info.images;

    std::atomic<std::uint64_t> next{0}, read_frames{0}, bad_total{0},
        join_total{0}, over_total{0};
    std::mutex best_mutex;
    std::uint64_t best = 0, best_frame = 0;
    unsigned bit_depth = 0;
    std::string failure;
    const double t0 = now();
    std::vector<std::thread> workers;
    for (std::size_t w = 0; w < threads; ++w)
      workers.emplace_back([&] {
        try {
          std::unique_ptr<series::Reader> reader = s->reader();
          series::Frame frame;
          std::vector<std::uint8_t> out;
          for (std::uint64_t f = next.fetch_add(1); f < frames;
               f = next.fetch_add(1)) {
            if (!reader->read(std::to_string(f), &frame))
              continue; // never written
            const std::size_t h = frame.height, wd = frame.width;
            const std::size_t bytes =
                decompress::frame_bytes(h, wd, frame.bit_depth);
            out.resize(bytes);
            decompress::image(frame.data, frame.algorithm, frame.bit_depth, h,
                              wd, {out.data(), out.size()});
            std::uint64_t largest = 0, bad = 0, joins = 0, over = 0;
            if (frame.bit_depth == 16)
              scan<std::uint16_t>(out.data(), h * wd, &largest, &bad, &joins,
                                  &over);
            else if (frame.bit_depth == 32)
              scan<std::uint32_t>(out.data(), h * wd, &largest, &bad, &joins,
                                  &over);
            else
              throw std::runtime_error(std::to_string(frame.bit_depth) +
                                       "-bit pixels");
            bad_total.fetch_add(bad);
            join_total.fetch_add(joins);
            over_total.fetch_add(over);
            read_frames.fetch_add(1);
            const std::lock_guard<std::mutex> lock(best_mutex);
            bit_depth = frame.bit_depth;
            if (largest > best || (largest == best && f < best_frame)) {
              best = largest;
              best_frame = f;
            }
          }
        } catch (const std::exception &e) {
          const std::lock_guard<std::mutex> lock(best_mutex);
          if (failure.empty())
            failure = e.what();
          next.store(frames);
        }
      });
    for (std::thread &t : workers)
      t.join();
    if (!failure.empty())
      throw std::runtime_error(failure);
    const double wall = now() - t0;

    const bool fits = best < kSixteenBitLimit;
    std::printf("mxi_max: %s, %llu frames of %u-bit pixels read in %.2f s on "
                "%zu thread%s\n",
                master.c_str(),
                static_cast<unsigned long long>(read_frames.load()), bit_depth,
                wall, threads, threads == 1 ? "" : "s");
    std::printf(
        "  the largest count on a valid pixel: %llu (0x%llx), on image %llu\n",
        static_cast<unsigned long long>(best),
        static_cast<unsigned long long>(best),
        static_cast<unsigned long long>(best_frame + 1));
    std::printf("  pixels marked bad: %llu; tile joins: %llu; valid pixels at "
                "0xFFFD or "
                "above: %llu\n",
                static_cast<unsigned long long>(bad_total.load()),
                static_cast<unsigned long long>(join_total.load()),
                static_cast<unsigned long long>(over_total.load()));
    if (bit_depth == 16)
      std::printf("  16-bit already: mxi_find --gpu takes it as it is\n");
    else if (fits)
      std::printf("  fits in 16 bits: mxi_find --gpu --gpu-force can take it "
                  "on a GPU of 16 "
                  "bits only\n");
    else
      std::printf("  does not fit in 16 bits: for a GPU of 16 bits only, "
                  "threshold it on the "
                  "CPU\n");
    return fits ? 0 : 1;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_max: %s\n", e.what());
    return 2;
  }
}
