#pragma once

// One row of a shoebox filled from a frame: the counts as floats, and the
// detector's markers recognised. A detector marks two kinds of pixel with the
// largest values of its width -- max() - 1 a bad pixel, max() a tile join, the
// gap between modules: 0xfffe and 0xffff in 16 bits, 0xfffffffe and 0xffffffff
// in 32 -- and neither is a count. The spot finder has always masked both; the
// integrator took only max() for a marker, so a bad pixel was filled as a count
// of 65534, and in a reflection's foreground added that to its intensity.

#include <cstddef>
#include <cstdint>
#include <limits>

#include "shoebox.hh"

namespace mxi {

//: The first value of a pixel's width that is a marker, not a count.
template <typename Pixel> constexpr Pixel lowest_marker() {
  return static_cast<Pixel>(std::numeric_limits<Pixel>::max() - 1);
}

//: Converts nx pixels to floats in `to`; a marker's voxel is left at 0 --
//: excluded from both sums, not counted as zero, which would drag the
//: background down wherever a gap crosses a shoebox -- and its VALIDITY
//: cleared, the region kept, since a voxel in a gap is still foreground with no
//: measurement in it and the profile fit must know part of the reflection is
//: missing. One pass with no branch, which the compiler makes SIMD, and a
//: second only for a row with a marker in it. Returns how many markers.
template <typename Pixel>
std::size_t fill_row(const Pixel *from, std::size_t nx, float *to,
                     std::uint8_t *mask) {
  const Pixel marker = lowest_marker<Pixel>();
  bool any = false;
  for (std::size_t x = 0; x < nx; ++x) {
    to[x] = static_cast<float>(from[x]);
    any |= from[x] >= marker;
  }
  if (!any)
    return 0;
  std::size_t count = 0;
  for (std::size_t x = 0; x < nx; ++x)
    if (from[x] >= marker) {
      to[x] = 0.0f;
      mask[x] &= static_cast<std::uint8_t>(~shoebox_mask::kValid);
      ++count;
    }
  return count;
}

} // namespace mxi
