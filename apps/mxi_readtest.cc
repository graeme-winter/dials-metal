// mxi_readtest: how fast an NXmx series can be read and decompressed, apart
// from everything a program then does with the pixels.
//
// Integration's frame reading on 32 threads fetched at 1.00 times the wall --
// one frame off the file at a time -- while `cat` read the same files from the
// page cache in 0.7 s. This reads every frame of the images an .expt names, on
// N threads, one of two ways, and decompresses each as the programs do:
//
//   --direct-chunk  through the series reader the programs use: HDF5's
//                   H5Dread_chunk, under the reader's lock
//   --pread         every chunk's place in its file found first, from HDF5,
//                   then each frame read with pread, outside any lock
//
// and reports the wall, the threads' time fetching and decompressing, and a
// checksum of every decompressed pixel, the same both ways or one of them is
// wrong. Nothing here is used by the pipeline; it is for measuring before the
// pipeline's reading is changed.

#include <hdf5.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <mutex>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

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

//: Where one frame's chunk is: the file descriptor to pread it from, its
//: address and size, and whether its filter was skipped.
struct Place {
  int fd = -1;
  std::uint64_t address = 0, size = 0;
  unsigned mask = 0;
};

//: Every frame's chunk, found through HDF5 before any frame is read: for each
//: block the dataset opened in its own file -- the external file, if the master
//: links to it -- its user block added to HDF5's addresses, and the file opened
//: for pread.
std::vector<Place> find_chunks(const std::vector<series::SourceBlock> &blocks,
                               std::uint64_t frames,
                               std::vector<int> *descriptors) {
  std::vector<Place> places(frames);
  for (const series::SourceBlock &b : blocks) {
    const hid_t file = H5Fopen(b.filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0)
      throw std::runtime_error("cannot open " + b.filename);
    const hid_t dataset = H5Dopen2(file, b.dataset.c_str(), H5P_DEFAULT);
    if (dataset < 0)
      throw std::runtime_error("cannot open " + b.dataset + " in " +
                               b.filename);
    const hid_t own = H5Iget_file_id(dataset);
    const ssize_t length = H5Fget_name(own, nullptr, 0);
    std::string name(static_cast<std::size_t>(length > 0 ? length : 0) + 1,
                     '\0');
    H5Fget_name(own, name.data(), name.size());
    name.resize(static_cast<std::size_t>(length > 0 ? length : 0));
    hsize_t user_block = 0;
    const hid_t fcpl = H5Fget_create_plist(own);
    H5Pget_userblock(fcpl, &user_block);
    H5Pclose(fcpl);
    const int fd = ::open(name.c_str(), O_RDONLY);
    if (fd < 0)
      throw std::runtime_error("cannot open " + name + " for pread");
    descriptors->push_back(fd);
    for (std::uint64_t f = b.first; f <= b.last && f < frames; ++f) {
      hsize_t coord[3] = {static_cast<hsize_t>(b.offset + (f - b.first)), 0, 0};
      unsigned mask = 0;
      haddr_t address = 0;
      hsize_t size = 0;
      if (H5Dget_chunk_info_by_coord(dataset, coord, &mask, &address, &size) <
          0)
        throw std::runtime_error(
            "H5Dget_chunk_info_by_coord failed for frame " + std::to_string(f));
      if (address != HADDR_UNDEF && size != 0)
        places[f] = {fd, static_cast<std::uint64_t>(address) + user_block,
                     static_cast<std::uint64_t>(size), mask};
    }
    H5Fclose(own);
    H5Dclose(dataset);
    H5Fclose(file);
  }
  return places;
}

bool pread_all(int fd, std::uint8_t *to, std::uint64_t size,
               std::uint64_t address) {
  std::uint64_t done = 0;
  while (done < size) {
    const ssize_t got =
        ::pread(fd, to + done, static_cast<std::size_t>(size - done),
                static_cast<off_t>(address + done));
    if (got < 0 && errno == EINTR)
      continue;
    if (got <= 0)
      return false;
    done += static_cast<std::uint64_t>(got);
  }
  return true;
}

//: FNV-1a over a frame's bytes; frames' hashes are summed, so the checksum does
//: not depend on which thread read which frame.
std::uint64_t fnv1a(const std::uint8_t *p, std::size_t n) {
  std::uint64_t h = 1469598103934665603ull;
  for (std::size_t i = 0; i < n; ++i) {
    h ^= p[i];
    h *= 1099511628211ull;
  }
  return h;
}

const char *kUsage =
    "usage: mxi_readtest EXPT (--direct-chunk | --pread) [options]\n"
    "  Reads and decompresses every frame of the images EXPT names, on N\n"
    "  threads, and reports how fast -- to measure reading apart from the\n"
    "  programs that use it.\n"
    "  --direct-chunk    through the series reader the programs use:\n"
    "                    H5Dread_chunk under the reader's lock\n"
    "  --pread           every chunk's place found first, then each frame\n"
    "                    read with pread, outside any lock\n"
    "  -j, --threads N   N threads (every core)\n"
    "  --frames N        the first N frames only\n"
    "  --no-decompress   fetch only; the checksum then over the compressed "
    "bytes\n"
    "  --images PATH     the NXmx master, not the .expt's imageset template\n";

} // namespace

int main(int argc, char **argv) {
  using namespace mxi;
  const std::set<std::string> known = {
      "--direct-chunk", "--pread",         "-j",       "--threads",
      "--frames",       "--no-decompress", "--images", "--help"};
  const std::set<std::string> takes_value = {"-j", "--threads", "--frames",
                                             "--images"};
  const Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (args.has("--help")) {
    std::printf("%s", kUsage);
    return 0;
  }
  if (!args.ok || args.positional.size() != 1 ||
      args.has("--direct-chunk") == args.has("--pread")) {
    if (!args.ok)
      std::fprintf(stderr, "mxi_readtest: %s\n", args.error.c_str());
    std::fprintf(stderr, "%s", kUsage);
    return 2;
  }
  try {
    const bool by_pread = args.has("--pread");
    const bool decompressing = !args.has("--no-decompress");
    std::string master = args.value("--images", "");
    if (master.empty())
      master = read_experiments(args.positional[0]).image_template;
    if (master.empty() || master.find('#') != std::string::npos)
      throw std::runtime_error("no single NXmx master to read; give --images");
    std::size_t threads = static_cast<std::size_t>(
        args.number("--threads", args.number("-j", 0.0)));
    if (threads == 0)
      threads = std::max(1u, std::thread::hardware_concurrency());

    // The geometry and compression, from the first frame through the reader
    // the programs use -- so that --pread decompresses as they would.
    std::unique_ptr<series::Series> s = series::nxmx(master);
    series::Info info;
    if (!s->try_open(&info))
      throw std::runtime_error("cannot open " + master + " as an NXmx series");
    std::uint64_t frames = info.images;
    if (args.has("--frames"))
      frames = std::min<std::uint64_t>(
          frames, static_cast<std::uint64_t>(args.number("--frames", 0.0)));
    series::Frame first;
    if (!s->reader()->read("0", &first))
      throw std::runtime_error("the first frame could not be read");
    const unsigned bit_depth = first.bit_depth;
    const std::size_t height = first.height, width = first.width;
    const std::string algorithm(first.algorithm);
    const std::size_t bytes = decompress::frame_bytes(height, width, bit_depth);

    std::vector<Place> places;
    std::vector<int> descriptors;
    double t_find = 0.0;
    if (by_pread) {
      const double t0 = now();
      series::Info again;
      places = find_chunks(series::nxmx_blocks(master, &again), frames,
                           &descriptors);
      t_find = now() - t0;
    }

    std::atomic<std::uint64_t> next{0}, checksum{0}, compressed{0},
        read_frames{0};
    std::vector<double> fetching(threads, 0.0), decompressing_s(threads, 0.0);
    const double t_start = now();
    std::mutex failure_mutex;
    std::string failure;
    std::vector<std::thread> workers;
    for (std::size_t w = 0; w < threads; ++w)
      workers.emplace_back([&, w] {
        try {
          std::unique_ptr<series::Reader> reader =
              by_pread ? nullptr : s->reader();
          std::vector<std::uint8_t> buffer, out(bytes);
          series::Frame frame;
          for (std::uint64_t f = next.fetch_add(1); f < frames;
               f = next.fetch_add(1)) {
            const double a = now();
            std::span<const std::uint8_t> data;
            std::string_view how = algorithm;
            if (by_pread) {
              const Place &p = places[f];
              if (p.size == 0)
                continue; // never written
              buffer.resize(p.size);
              if (!pread_all(p.fd, buffer.data(), p.size, p.address))
                throw std::runtime_error("pread failed for frame " +
                                         std::to_string(f));
              data = {buffer.data(), buffer.size()};
              if (p.mask & 1u)
                how = std::string_view();
            } else {
              if (!reader->read(std::to_string(f), &frame))
                continue;
              data = frame.data;
              how = frame.algorithm;
            }
            const double b = now();
            fetching[w] += b - a;
            compressed.fetch_add(data.size());
            std::uint64_t h = 0;
            if (decompressing) {
              decompress::image(data, how, bit_depth, height, width,
                                {out.data(), out.size()});
              h = fnv1a(out.data(), out.size());
              decompressing_s[w] += now() - b;
            } else {
              h = fnv1a(data.data(), data.size());
            }
            checksum.fetch_add(h);
            read_frames.fetch_add(1);
          }
        } catch (const std::exception &e) {
          const std::lock_guard<std::mutex> lock(failure_mutex);
          if (failure.empty())
            failure = e.what();
          next.store(frames); // the others stop
        }
      });
    for (std::thread &t : workers)
      t.join();
    if (!failure.empty())
      throw std::runtime_error(failure);
    const double wall = now() - t_start;
    for (int fd : descriptors)
      ::close(fd);

    double fetch = 0.0, decomp = 0.0;
    for (std::size_t w = 0; w < threads; ++w) {
      fetch += fetching[w];
      decomp += decompressing_s[w];
    }
    const double gb = static_cast<double>(compressed.load()) / 1e9;
    std::printf("mxi_readtest: %s, %llu frames of %zu x %zu, %u bit, %s\n",
                master.c_str(),
                static_cast<unsigned long long>(read_frames.load()), width,
                height, bit_depth,
                algorithm.empty() ? "uncompressed" : algorithm.c_str());
    std::printf("  read %s on %zu threads%s\n",
                by_pread ? "by pread, outside any lock"
                         : "by H5Dread_chunk, the reader's way",
                threads,
                decompressing ? ", and decompressed" : ", not decompressed");
    if (by_pread)
      std::printf("  finding the chunks     %8.3f s\n", t_find);
    std::printf("  wall                   %8.3f s   %.0f frames/s, %.2f GB/s "
                "compressed\n",
                wall, static_cast<double>(read_frames.load()) / wall,
                gb / wall);
    std::printf("  fetching               %8.3f thread-s   (%.2f x wall)\n",
                fetch, fetch / wall);
    if (decompressing)
      std::printf("  decompressing          %8.3f thread-s   (%.2f x wall)\n",
                  decomp, decomp / wall);
    std::printf("  compressed             %8.3f GB\n", gb);
    std::printf("  checksum               %016llx\n",
                static_cast<unsigned long long>(checksum.load()));
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_readtest: %s\n", e.what());
    return 1;
  }
}
