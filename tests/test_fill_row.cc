#include <cstdint>
#include <vector>

#include "../src/fill_row.hh"
#include "check.hh"

namespace mxi {

namespace {

template <typename Pixel> void check_row(const char *width) {
  const Pixel top = std::numeric_limits<Pixel>::max();
  // A count, the largest real count (max - 2), the bad pixel (max - 1), the
  // tile join (max), and a zero.
  const std::vector<Pixel> row = {7, static_cast<Pixel>(top - 2),
                                  static_cast<Pixel>(top - 1), top, 0};
  std::vector<float> to(row.size(), -1.0f);
  std::vector<std::uint8_t> mask(row.size(), shoebox_mask::kValid |
                                                 shoebox_mask::kForeground);
  const std::size_t markers =
      fill_row(row.data(), row.size(), to.data(), mask.data());
  const std::string w = width;
  check::equal(static_cast<long long>(markers), 2,
               w + ": the bad pixel and the tile join");
  check::is_true(to[0] == 7.0f && to[1] == static_cast<float>(top - 2) &&
                     to[4] == 0.0f,
                 w + ": counts are counts, the largest real one included");
  check::is_true(to[2] == 0.0f && to[3] == 0.0f,
                 w + ": a marker's voxel is 0, not a count");
  for (std::size_t x = 0; x < row.size(); ++x) {
    const bool valid = (mask[x] & shoebox_mask::kValid) != 0;
    const bool foreground = (mask[x] & shoebox_mask::kForeground) != 0;
    check::is_true(valid == (x != 2 && x != 3),
                   w + ": only the markers lose validity");
    check::is_true(foreground, w + ": the region kept");
  }
  const std::vector<Pixel> clean = {1, 2, 3};
  std::vector<float> to2(3);
  std::vector<std::uint8_t> mask2(3, shoebox_mask::kValid);
  check::equal(static_cast<long long>(
                   fill_row(clean.data(), 3, to2.data(), mask2.data())),
               0, w + ": a row with no marker");
  check::is_true(mask2 == std::vector<std::uint8_t>(3, shoebox_mask::kValid),
                 w + ": untouched");
}

} // namespace

TEST(both_markers_are_invalid_in_a_filled_row_and_counts_are_counts) {
  // A bad pixel is max() - 1 and a tile join max(), in 16 bits and 32. The
  // integrator once took only max() for a marker, and filled a bad pixel as a
  // count of 65534.
  check_row<std::uint16_t>("16 bits");
  check_row<std::uint32_t>("32 bits");
}

} // namespace mxi
