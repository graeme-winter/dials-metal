// The reflection-table writer's output, held to what it produces today.
//
// This exists so that `spotfinder/src/refl.cc` can be replaced by the general
// writer in `src/refl.cc` without anyone having to take it on trust. A table
// that reads back is not evidence: the writer and the reader here share every
// assumption, which is exactly how five format bugs survived in this
// repository until DIALS read the output.
//
// If a rewrite changes these bytes, the difference has to be explained. If it
// cannot be byte identical -- the two writers order their fields differently
// -- then this test should be changed to decode both and compare columns,
// dtypes and shoebox contents. It must not be weakened to pass.

#include <cstdio>
#include <string>
#include <vector>

#include "dials_spots.hh"
#include "refl.hh"

namespace {

int failures = 0;

void check(bool condition, const std::string &what) {
  if (!condition) {
    std::printf("  FAIL %s\n", what.c_str());
    ++failures;
  }
}

// A fixed, arbitrary set of spots: repeatability is the point, not realism.
void build(std::vector<dials_spots::Spot> *spots,
           std::vector<dials_spots::Pixel> *pixels, std::size_t width) {
  std::size_t at = 0;
  for (int s = 0; s < 5; ++s) {
    dials_spots::Spot spot;
    spot.first = at;
    spot.n_signal = 0;
    spot.bbox[0] = 12 + s * 5;
    spot.bbox[1] = spot.bbox[0] + 3;
    spot.bbox[2] = 10 + s * 7;
    spot.bbox[3] = spot.bbox[2] + 3;
    spot.bbox[4] = s;
    spot.bbox[5] = s + 2 + s % 2;
    spot.position[0] = spot.bbox[0] + 1.5;
    spot.position[1] = spot.bbox[2] + 1.5;
    spot.position[2] = s + 0.5;
    spot.variance[0] = 0.4;
    spot.variance[1] = 0.5;
    spot.variance[2] = 0.6;
    spot.intensity = 0.0;
    for (int z = 0; z < 2 + s % 2; ++z) {
      for (int y = 0; y < 3; ++y) {
        for (int x = 0; x < 3; ++x) {
          dials_spots::Pixel p;
          p.index = static_cast<std::uint32_t>((10 + s * 7 + y) * width + 12 +
                                               s * 5 + x);
          p.frame = static_cast<std::int32_t>(z + s);
          p.value = static_cast<std::uint32_t>(3 + (x + y + z + s) % 7);
          pixels->push_back(p);
          ++at;
          ++spot.n_signal;
          spot.intensity += p.value;
        }
      }
    }
    spot.intensity_variance = spot.intensity;
    spots->push_back(spot);
  }
}

} // namespace

int main() {
  std::vector<dials_spots::Spot> spots;
  std::vector<dials_spots::Pixel> pixels;
  const std::size_t width = 64;
  build(&spots, &pixels, width);

  refl::Options options;
  options.panel = 0;
  options.id = 0;
  options.identifier = "0123abcd-0000-1111-2222-333344445555";
  options.shoeboxes = true;

  const std::string path = "test_refl_golden.refl";
  refl::write(path, spots, pixels, width, options);

  std::FILE *file = std::fopen(path.c_str(), "rb");
  check(file != nullptr, "the file was written");
  if (file == nullptr) {
    std::printf("FAIL: the golden reflection table, %d failures\n", ++failures);
    return 1;
  }
  std::string bytes;
  char buffer[4096];
  std::size_t got = 0;
  while ((got = std::fread(buffer, 1, sizeof(buffer), file)) > 0) {
    bytes.append(buffer, got);
  }
  std::fclose(file);

  // The size and a cheap checksum, which between them catch a changed field,
  // a changed dtype and a changed order without carrying a megabyte of hex in
  // the source.
  std::size_t checksum = 1469598103934665603ULL;
  for (unsigned char c : bytes) {
    checksum ^= c;
    checksum *= 1099511628211ULL;
  }
  const std::size_t expected_size = 2026UL;
  const std::size_t expected_checksum = 11545107909170840019UL;
  check(bytes.size() == expected_size,
        "size is " + std::to_string(bytes.size()) + ", was " +
            std::to_string(expected_size));
  check(checksum == expected_checksum, "checksum is " +
                                           std::to_string(checksum) + ", was " +
                                           std::to_string(expected_checksum));

  std::printf("%s: the golden reflection table, %d failures\n",
              failures == 0 ? "PASS" : "FAIL", failures);
  return failures == 0 ? 0 : 1;
}
