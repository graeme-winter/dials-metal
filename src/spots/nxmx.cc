// series::nxmx -- read a series out of an NXmx HDF5 file.
//
// The master file's /entry/data/data is a virtual dataset spanning a set of
// data files. A virtual dataset cannot be chunk-read, so the VDS mapping is
// unpacked instead: it names its own source files, datasets and offsets, which
// means nothing here has to assume how the writer laid the series out. Frames
// are then read with H5Dread_chunk straight out of the data file that holds
// them, so the compressed bytes reach the threshold without HDF5's filter
// pipeline ever running -- the same path the data took in, in reverse.

#include <hdf5.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

#include "series.hh"

namespace series {
namespace {

constexpr const char *kImageData = "/entry/data/data";

// HDF5 filter identifiers, as registered with The HDF Group.
constexpr H5Z_filter_t kBitshuffle = 32008;
constexpr H5Z_filter_t kLz4 = 32004;

// The HDF5 C library keeps global state -- identifier tables, free lists, the
// metadata cache -- and is not thread-safe unless it was built that way, which
// vcpkg's default port is not and several distributions are. Giving each thread
// its own file handle is not enough: two threads inside the library at once
// will corrupt something and crash, which is exactly what happens on a build
// without thread safety.
//
// So every call below is serialised, including the closes in Handle. Nothing is
// lost by doing it unconditionally, because a thread-safe HDF5 takes one global
// lock across its whole public API anyway -- the reads were already serialised
// on those builds. What still runs in parallel is the part that costs:
// decompression and analysis, downstream of here.
//
// Recursive because the scopes nest -- read() holds it across open(), and
// open() holds it while its temporary handles close.
std::recursive_mutex &hdf5_mutex() {
  static std::recursive_mutex mutex;
  return mutex;
}

using Guard = std::lock_guard<std::recursive_mutex>;

void ok(herr_t status, const char *what) {
  if (status < 0)
    throw std::runtime_error(std::string("HDF5: ") + what + " failed");
}

// HDF5 2.0 added a buffer-size argument to H5Dread_chunk -- in, the buffer's
// capacity; out, what the chunk needed -- and its versioned-API macros point
// the old name at the new function, so the call cannot be written once for
// both. Which shape this HDF5 has is settled by a compile probe in
// CMakeLists.txt, since vcpkg and the distributions are currently on different
// sides of it.
herr_t read_chunk(hid_t dataset, const hsize_t *offset, std::uint32_t *filters,
                  void *buffer, std::size_t capacity) {
#ifdef SPOTFINDER_H5DREAD_CHUNK_TAKES_SIZE
  std::size_t size = capacity;
  return H5Dread_chunk(dataset, H5P_DEFAULT, offset, filters, buffer, &size);
#else
  (void)capacity;
  return H5Dread_chunk(dataset, H5P_DEFAULT, offset, filters, buffer);
#endif
}

hid_t opened(hid_t id, const char *what) {
  if (id < 0)
    throw std::runtime_error(std::string("HDF5: ") + what + " failed");
  return id;
}

// HDF5 hands out integer handles with a different close function for each kind
// of object, so the closer travels with the handle.
class Handle {
public:
  Handle() = default;
  Handle(hid_t id, herr_t (*close)(hid_t)) : id_(id), close_(close) {}
  ~Handle() { reset(); }

  Handle(Handle &&other) noexcept : id_(other.id_), close_(other.close_) {
    other.id_ = -1;
  }
  Handle &operator=(Handle &&other) noexcept {
    if (this != &other) {
      reset();
      id_ = other.id_;
      close_ = other.close_;
      other.id_ = -1;
    }
    return *this;
  }
  Handle(const Handle &) = delete;
  Handle &operator=(const Handle &) = delete;

  hid_t get() const { return id_; }

private:
  void reset() {
    if (id_ >= 0 && close_ != nullptr) {
      const Guard guard(hdf5_mutex());
      close_(id_);
    }
    id_ = -1;
  }

  hid_t id_ = -1;
  herr_t (*close_)(hid_t) = nullptr;
};

// Where one run of frames lives. One of these per data file, in the usual case
// of a VDS built block by block.
struct Block {
  std::string filename; // resolved against the master file's directory
  std::string dataset;
  std::uint64_t first = 0;  // first frame index in the series
  std::uint64_t last = 0;   // last frame index, inclusive
  std::uint64_t offset = 0; // index of `first` within the data file
};

// The two strings that identify a source dataset, so it can be opened once per
// thread rather than once per frame.
std::string cache_key(const Block &block) {
  return block.filename + "\n" + block.dataset;
}

std::string plist_string(hid_t plist, int index,
                         ssize_t (*get)(hid_t, size_t, char *, size_t)) {
  const ssize_t size = get(plist, static_cast<size_t>(index), nullptr, 0);
  if (size < 0)
    throw std::runtime_error("HDF5: cannot size a virtual mapping name");
  // The size HDF5 reports excludes the terminator but the buffer size it wants
  // includes it, so the string is allocated one longer and trimmed after.
  std::string text(static_cast<std::size_t>(size) + 1, '\0');
  if (get(plist, static_cast<size_t>(index), text.data(), text.size()) < 0) {
    throw std::runtime_error("HDF5: cannot read a virtual mapping name");
  }
  text.resize(static_cast<std::size_t>(size));
  return text;
}

// Which compression the dataset declares. Only the identifier matters: the
// bytes are handed on untouched, so nothing here needs the filter to be
// installed.
std::string_view algorithm_of(hid_t dcpl) {
  const int filters = H5Pget_nfilters(dcpl);
  if (filters < 0)
    throw std::runtime_error("HDF5: cannot count filters");
  for (int i = 0; i < filters; i++) {
    unsigned flags = 0;
    size_t values = 0;
    const H5Z_filter_t filter =
        H5Pget_filter2(dcpl, static_cast<unsigned>(i), &flags, &values, nullptr,
                       0, nullptr, nullptr);
    if (filter == kBitshuffle)
      return "bslz4";
    if (filter == kLz4)
      return "lz4";
  }
  return {};
}

//: Where one frame's compressed chunk is in its file, read once from HDF5.
struct ChunkPlace {
  std::uint64_t address =
      0; //: from the start of the file, the user block included
  std::uint64_t size = 0; //: 0 for a chunk never written
  unsigned mask = 0;      //: filters skipped, as H5Dread_chunk reports them
};

//: A file descriptor, closed with its owner.
struct Descriptor {
  int fd = -1;
  Descriptor() = default;
  Descriptor(const Descriptor &) = delete;
  Descriptor &operator=(const Descriptor &) = delete;
  Descriptor(Descriptor &&other) noexcept : fd(other.fd) { other.fd = -1; }
  Descriptor &operator=(Descriptor &&other) noexcept {
    std::swap(fd, other.fd);
    return *this;
  }
  ~Descriptor() {
    if (fd >= 0)
      ::close(fd);
  }
};

struct SourceDataset {
  Handle file;
  Handle dataset;
  std::uint64_t height = 0;
  std::uint64_t width = 0;
  unsigned bit_depth = 0;
  std::string_view algorithm;
  // Reading without HDF5: every chunk's place, found once when the dataset is
  // opened, and the file it is in, opened for pread. Empty when that cannot
  // be done -- another file driver, a file that will not open -- and the
  // frames are read through HDF5 instead.
  std::vector<ChunkPlace> chunks;
  Descriptor direct;
};

std::atomic<std::uint64_t> g_direct_reads{0}, g_hdf5_reads{0};

//: Whether frames may be read without HDF5: MXI_HDF5_DIRECT=0 says not, to
//: compare against.
bool direct_reads_allowed() {
  const char *asked = std::getenv("MXI_HDF5_DIRECT");
  return !(asked && std::string(asked) == "0");
}

//: The whole of `size` bytes from `fd` at `address`, however many preads that
//: takes.
bool pread_all(int fd, std::uint8_t *to, std::uint64_t size,
               std::uint64_t address) {
  std::uint64_t done = 0;
  while (done < size) {
    const ssize_t got =
        ::pread(fd, to + done, static_cast<std::size_t>(size - done),
                static_cast<off_t>(address + done));
    if (got < 0) {
      if (errno == EINTR)
        continue;
      return false;
    }
    if (got == 0)
      return false;
    done += static_cast<std::uint64_t>(got);
  }
  return true;
}

class NxmxReader : public Reader {
public:
  explicit NxmxReader(std::vector<Block> blocks) : blocks_(std::move(blocks)) {}

  bool read(const std::string &key, Frame *frame) override {
    const std::uint64_t index = std::strtoull(key.c_str(), nullptr, 10);
    const Block &block = locate(index);
    SourceDataset *found_source = nullptr;
    {
      const Guard guard(hdf5_mutex());
      found_source = &open(block);
    }
    SourceDataset &source = *found_source;
    const std::uint64_t within = index - block.first + block.offset;

    // Without HDF5, and so without its lock: the chunk's place was found when
    // the dataset was opened, and the bytes are the ones H5Dread_chunk would
    // give. HDF5 serialises every call through one lock, so on 32 threads the
    // frames came off the file one at a time -- 43.6 thread-seconds fetching
    // in an integration of 1800 frames, 13 times the wall.
    if (source.direct.fd >= 0 && within < source.chunks.size()) {
      const ChunkPlace &place = source.chunks[static_cast<std::size_t>(within)];
      if (place.size == 0)
        return false; // never written, as below
      frame->storage.resize(static_cast<std::size_t>(place.size));
      if (!pread_all(source.direct.fd,
                     reinterpret_cast<std::uint8_t *>(frame->storage.data()),
                     place.size, place.address))
        throw std::runtime_error("reading frame " + std::to_string(index) +
                                 " from " + block.filename +
                                 " without HDF5 failed");
      frame->height = source.height;
      frame->width = source.width;
      frame->bit_depth = source.bit_depth;
      frame->algorithm =
          (place.mask & 1u) ? std::string_view() : source.algorithm;
      frame->data = {
          reinterpret_cast<const std::uint8_t *>(frame->storage.data()),
          frame->storage.size()};
      frame->number = static_cast<std::int64_t>(index);
      g_direct_reads.fetch_add(1, std::memory_order_relaxed);
      return true;
    }

    const Guard guard(hdf5_mutex());
    g_hdf5_reads.fetch_add(1, std::memory_order_relaxed);

    // The chunk holds one frame, so its coordinate is the frame's own index
    // within the data file, and the other two are zero by construction.
    hsize_t offset[3] = {static_cast<hsize_t>(within), 0, 0};

    unsigned mask = 0;
    haddr_t address = 0;
    hsize_t size = 0;
    ok(H5Dget_chunk_info_by_coord(source.dataset.get(), offset, &mask, &address,
                                  &size),
       "H5Dget_chunk_info_by_coord");

    // A chunk that was never allocated is a frame the writer never received.
    // It reads back as the fill value, which is not something to analyse.
    if (address == HADDR_UNDEF || size == 0)
      return false;

    frame->storage.resize(static_cast<std::size_t>(size));
    std::uint32_t filters = 0;
    ok(read_chunk(source.dataset.get(), offset, &filters, frame->storage.data(),
                  frame->storage.size()),
       "H5Dread_chunk");

    frame->height = source.height;
    frame->width = source.width;
    frame->bit_depth = source.bit_depth;
    // Bit i of the mask means filter i was skipped for this chunk, which is how
    // a writer stores a payload it could not compress. The algorithm string
    // belongs to the source dataset, which this reader owns and which outlives
    // every frame it fills.
    frame->algorithm = (filters & 1u) ? std::string_view() : source.algorithm;
    frame->data = {
        reinterpret_cast<const std::uint8_t *>(frame->storage.data()),
        frame->storage.size()};
    frame->number = static_cast<std::int64_t>(index);
    return true;
  }

private:
  const Block &locate(std::uint64_t index) const {
    for (const Block &block : blocks_) {
      if (index >= block.first && index <= block.last)
        return block;
    }
    throw std::runtime_error("frame " + std::to_string(index) +
                             " is outside every virtual mapping");
  }

  // Called with the HDF5 lock already held, which the recursive mutex allows.
  SourceDataset &open(const Block &block) {
    const Guard guard(hdf5_mutex());
    const std::string key = cache_key(block);
    const auto found = open_.find(key);
    if (found != open_.end())
      return found->second;

    SourceDataset source;
    source.file = Handle(
        opened(H5Fopen(block.filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT),
               ("opening " + block.filename).c_str()),
        H5Fclose);
    source.dataset = Handle(
        opened(H5Dopen2(source.file.get(), block.dataset.c_str(), H5P_DEFAULT),
               ("opening " + block.dataset).c_str()),
        H5Dclose);

    const Handle dcpl(opened(H5Dget_create_plist(source.dataset.get()),
                             "H5Dget_create_plist"),
                      H5Pclose);
    if (H5Pget_layout(dcpl.get()) != H5D_CHUNKED) {
      throw std::runtime_error(block.filename + " " + block.dataset +
                               " is not chunked, so it cannot be chunk-read");
    }
    hsize_t chunk[3] = {0, 0, 0};
    if (H5Pget_chunk(dcpl.get(), 3, chunk) != 3) {
      throw std::runtime_error("expected a three-dimensional chunk in " +
                               block.filename);
    }
    if (chunk[0] != 1) {
      throw std::runtime_error(
          "chunks in " + block.filename + " span " + std::to_string(chunk[0]) +
          " frames; only one frame per chunk can be read this way");
    }
    source.height = chunk[1];
    source.width = chunk[2];
    source.algorithm = algorithm_of(dcpl.get());

    const Handle type(opened(H5Dget_type(source.dataset.get()), "H5Dget_type"),
                      H5Tclose);
    source.bit_depth = static_cast<unsigned>(8 * H5Tget_size(type.get()));

    if (direct_reads_allowed())
      place_chunks(&source);

    return open_.emplace(key, std::move(source)).first->second;
  }

  // Every chunk's place in the dataset's own file, and that file opened for
  // pread -- or nothing, and the dataset read through HDF5. Only for the
  // default file driver, whose addresses are offsets into one file; the user
  // block, if any, is added, since HDF5's addresses start after it. The file is
  // the dataset's own, as HDF5 names it: the external file, when the master
  // links to its data.
  static void place_chunks(SourceDataset *source) {
    const hid_t dataset = source->dataset.get();
    const Handle fapl(H5Fget_access_plist(source->file.get()), H5Pclose);
    if (fapl.get() < 0 || H5Pget_driver(fapl.get()) != H5FD_SEC2)
      return;
    const hid_t file_of_dataset = H5Iget_file_id(dataset);
    if (file_of_dataset < 0)
      return;
    const Handle own_file(file_of_dataset, H5Fclose);
    const Handle own_fapl(H5Fget_access_plist(own_file.get()), H5Pclose);
    if (own_fapl.get() < 0 || H5Pget_driver(own_fapl.get()) != H5FD_SEC2)
      return;
    const ssize_t length = H5Fget_name(own_file.get(), nullptr, 0);
    if (length <= 0)
      return;
    std::string name(static_cast<std::size_t>(length) + 1, '\0');
    H5Fget_name(own_file.get(), name.data(), name.size());
    name.resize(static_cast<std::size_t>(length));
    hsize_t user_block = 0;
    {
      const Handle fcpl(H5Fget_create_plist(own_file.get()), H5Pclose);
      if (fcpl.get() < 0 || H5Pget_userblock(fcpl.get(), &user_block) < 0)
        return;
    }
    const Handle space(H5Dget_space(dataset), H5Sclose);
    hsize_t dims[3] = {0, 0, 0};
    if (space.get() < 0 || H5Sget_simple_extent_ndims(space.get()) != 3 ||
        H5Sget_simple_extent_dims(space.get(), dims, nullptr) != 3)
      return;
    std::vector<ChunkPlace> chunks(static_cast<std::size_t>(dims[0]));
    for (hsize_t i = 0; i < dims[0]; ++i) {
      hsize_t offset[3] = {i, 0, 0};
      unsigned mask = 0;
      haddr_t address = 0;
      hsize_t size = 0;
      if (H5Dget_chunk_info_by_coord(dataset, offset, &mask, &address, &size) <
          0)
        return;
      if (address == HADDR_UNDEF || size == 0)
        continue; // never written: size 0
      chunks[static_cast<std::size_t>(i)] = {
          static_cast<std::uint64_t>(address) + user_block,
          static_cast<std::uint64_t>(size), mask};
    }
    Descriptor direct;
    direct.fd = ::open(name.c_str(), O_RDONLY);
    if (direct.fd < 0)
      return;
    source->chunks = std::move(chunks);
    source->direct = std::move(direct);
  }

  std::vector<Block> blocks_;
  std::map<std::string, SourceDataset> open_;
};

class Nxmx : public Series {
public:
  explicit Nxmx(std::string master) : master_(std::move(master)) {}

  bool try_open(Info *info) override {
    std::error_code ignored;
    if (!std::filesystem::is_regular_file(master_, ignored))
      return false;

    const Guard guard(hdf5_mutex());

    const Handle file(
        opened(H5Fopen(master_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT),
               ("opening " + master_).c_str()),
        H5Fclose);
    const Handle data(opened(H5Dopen2(file.get(), kImageData, H5P_DEFAULT),
                             ("opening " + std::string(kImageData)).c_str()),
                      H5Dclose);

    const Handle space(opened(H5Dget_space(data.get()), "H5Dget_space"),
                       H5Sclose);
    if (H5Sget_simple_extent_ndims(space.get()) != 3) {
      throw std::runtime_error(std::string(kImageData) +
                               " is not a stack of images");
    }
    hsize_t dims[3] = {0, 0, 0};
    ok(H5Sget_simple_extent_dims(space.get(), dims, nullptr),
       "H5Sget_simple_extent_dims");

    const Handle dcpl(
        opened(H5Dget_create_plist(data.get()), "H5Dget_create_plist"),
        H5Pclose);
    blocks_ = H5Pget_layout(dcpl.get()) == H5D_VIRTUAL
                  ? unpack_virtual(dcpl.get(), dims[0])
                  : whole_dataset(dims[0]);

    info->name = std::filesystem::path(master_).stem().string();
    info->images = dims[0];
    info->height = dims[1];
    info->width = dims[2];
    return true;
  }

  // Every frame is there to be read as soon as the file opens: there is nothing
  // to wait for once the master file exists.
  bool restrict_frames(std::uint64_t first, std::uint64_t last) override {
    first_frame_ = first;
    last_frame_ = last;
    return true;
  }

  std::vector<std::string> ready() override {
    std::vector<std::string> keys;
    if (dispatched_)
      return keys;
    dispatched_ = true;
    keys.reserve(static_cast<std::size_t>(images_of_blocks()));
    for (std::uint64_t i = 0; i < images_of_blocks(); i++) {
      // Only the frames the scan covers, when one was given. The key is the
      // array index, so the range is a filter on it and costs nothing.
      if (i < first_frame_ || i > last_frame_)
        continue;
      char key[24];
      std::snprintf(key, sizeof(key), "%06llu",
                    static_cast<unsigned long long>(i));
      keys.emplace_back(key);
    }
    return keys;
  }

  bool finished() override { return true; }

  std::string describe() const override { return master_; }

  std::unique_ptr<Reader> reader() override {
    return std::unique_ptr<Reader>(new NxmxReader(blocks_));
  }

private:
  std::uint64_t images_of_blocks() const {
    std::uint64_t images = 0;
    for (const Block &block : blocks_)
      images = std::max(images, block.last + 1);
    return images;
  }

  // Read the mapping out of the dataset creation property list. Each mapping's
  // selection in the virtual dataset gives the range of frames it covers, and
  // the matching selection in the source gives where they start in the data
  // file; that is all that is needed to turn a frame number into a chunk
  // coordinate.
  // Called from try_open with the HDF5 lock held.
  std::vector<Block> unpack_virtual(hid_t dcpl, std::uint64_t images) {
    size_t count = 0;
    ok(H5Pget_virtual_count(dcpl, &count), "H5Pget_virtual_count");
    if (count == 0)
      throw std::runtime_error(std::string(kImageData) + " maps to nothing");

    const std::filesystem::path directory =
        std::filesystem::path(master_).parent_path();

    std::vector<Block> blocks;
    for (size_t i = 0; i < count; i++) {
      Block block;
      const std::string filename =
          plist_string(dcpl, static_cast<int>(i), H5Pget_virtual_filename);
      block.dataset =
          plist_string(dcpl, static_cast<int>(i), H5Pget_virtual_dsetname);

      // "." means the mapping points back into the master file itself, and
      // the master's path is already whatever the caller gave -- it must NOT
      // be resolved against its own directory a second time. Doing so turned
      // "../ins10_1.nxs" into "../../ins10_1.nxs" and failed to open, for any
      // master named by a relative path with a directory in it.
      //
      // Only a genuinely external filename, which HDF5 stores relative to the
      // master, needs the directory prepended.
      std::filesystem::path path;
      if (filename == ".") {
        path = std::filesystem::path(master_);
      } else {
        path = std::filesystem::path(filename);
        if (path.is_relative() && !directory.empty())
          path = directory / path;
      }
      block.filename = path.string();

      const Handle vspace(
          opened(H5Pget_virtual_vspace(dcpl, static_cast<unsigned>(i)),
                 "H5Pget_virtual_vspace"),
          H5Sclose);
      const Handle source(
          opened(H5Pget_virtual_srcspace(dcpl, static_cast<unsigned>(i)),
                 "H5Pget_virtual_srcspace"),
          H5Sclose);

      hsize_t start[3] = {0, 0, 0};
      hsize_t end[3] = {0, 0, 0};
      ok(H5Sget_select_bounds(vspace.get(), start, end),
         "H5Sget_select_bounds on a virtual selection");
      block.first = start[0];
      // An unlimited selection reports an enormous upper bound, so the
      // dataset's own extent is the authority on where the series stops.
      block.last = std::min<std::uint64_t>(end[0], images - 1);

      hsize_t source_start[3] = {0, 0, 0};
      hsize_t source_end[3] = {0, 0, 0};
      ok(H5Sget_select_bounds(source.get(), source_start, source_end),
         "H5Sget_select_bounds on a source selection");
      block.offset = source_start[0];

      blocks.push_back(block);
    }
    std::sort(blocks.begin(), blocks.end(),
              [](const Block &a, const Block &b) { return a.first < b.first; });
    return blocks;
  }

  std::vector<Block> whole_dataset(std::uint64_t images) {
    Block block;
    block.filename = master_;
    block.dataset = kImageData;
    block.first = 0;
    block.last = images == 0 ? 0 : images - 1;
    block.offset = 0;
    return {block};
  }

  std::string master_;
  std::vector<Block> blocks_;
  bool dispatched_ = false;
  std::uint64_t first_frame_ = 0;
  std::uint64_t last_frame_ = std::numeric_limits<std::uint64_t>::max();
};

} // namespace

void nxmx_read_counts(std::uint64_t *direct, std::uint64_t *through_hdf5) {
  *direct = g_direct_reads.load();
  *through_hdf5 = g_hdf5_reads.load();
}

std::unique_ptr<Series> nxmx(std::string master) {
  return std::unique_ptr<Series>(new Nxmx(std::move(master)));
}

} // namespace series
