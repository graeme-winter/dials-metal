// Shoeboxes: the pixels DIALS carries in the reflection table itself.
//
// This package said for a long time that it had no pixels, and it was wrong.
// `strong.refl` carries a `Shoebox<>` column -- 13.6 MB of it for insulin, one
// box per strong spot -- which is exactly the data the profile model needs.
// Everything up to here read the table and skipped that column.
//
// The layout was derived from a real dials.find_spots file and checked against
// every one of its 13766 records: the bounding boxes match the table's own
// `bbox` column and the records consume the blob exactly, to the byte.
//
//     int32          panel
//     int32 x 6      bbox as x0, x1, y0, y1, z0, z1
//     uint8          a flag, 2 in every record seen
//     float32 x N    data          N = (x1-x0)(y1-y0)(z1-z0)
//     uint8  x N     mask
//     float32 x N    background    all zero out of dials.find_spots
//
// The bounds are half open, as everywhere else here: a box from 2010 to 2014
// is four pixels wide.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "refl.hh"

namespace mxi {

// DIALS' mask codes, as they appear in a shoebox. A strong spot out of
// dials.find_spots carries 0 and 5, which is Valid | Foreground.
namespace shoebox_mask {
constexpr std::uint8_t kValid = 1 << 0;
constexpr std::uint8_t kBackground = 1 << 1;
constexpr std::uint8_t kForeground = 1 << 2;
constexpr std::uint8_t kOverlapped = 1 << 3;
} // namespace shoebox_mask

struct Shoebox {
  std::int32_t panel = 0;
  //: x0, x1, y0, y1, z0, z1, half open.
  std::int32_t bbox[6] = {0, 0, 0, 0, 0, 0};
  std::uint8_t flag = 0;
  std::vector<float> data;
  std::vector<std::uint8_t> mask;
  std::vector<float> background;

  std::int32_t nx() const { return bbox[1] - bbox[0]; }
  std::int32_t ny() const { return bbox[3] - bbox[2]; }
  std::int32_t nz() const { return bbox[5] - bbox[4]; }
  std::size_t size() const {
    return static_cast<std::size_t>(nx()) * static_cast<std::size_t>(ny()) *
           static_cast<std::size_t>(nz());
  }
  //: Index into `data` for a voxel at (x, y, z) relative to the box corner.
  //: Slowest to fastest is z, y, x, which is the order the bytes are in.
  std::size_t at(std::int32_t x, std::int32_t y, std::int32_t z) const {
    return (static_cast<std::size_t>(z) * static_cast<std::size_t>(ny()) +
            static_cast<std::size_t>(y)) *
               static_cast<std::size_t>(nx()) +
           static_cast<std::size_t>(x);
  }
};

//: Decode the `shoebox` column of a table read from a file.
//:
//: Throws if the blob does not parse exactly: a layout that consumes all but a
//: few bytes is a layout that is wrong somewhere earlier, and reporting the
//: boxes it managed to read would be worse than refusing.
std::vector<Shoebox> decode_shoeboxes(const Table &table);

//: Encode, for tests and for writing a table with shoeboxes built here.
//: Put a shoebox's mask into DIALS' convention: a voxel with no measurement
//: in it is zero, with no region bit.
//:
//: Internally a bad pixel keeps its Foreground or Background bit and loses
//: only Valid, so a profile fit can tell a foreground voxel with nothing in it
//: from one that was never foreground. DIALS' mask calculator sets those bits
//: only on voxels that are already Valid, so the combination never occurs
//: there, and a table carrying it was rejected as an invalid structure. Call
//: this on anything going to a file DIALS will read.
//:
//: What is lost is recoverable: the foreground region is geometry, and can be
//: rebuilt from the bounding box and the profile model.
void to_dials_convention(Shoebox *box);

std::string encode_shoeboxes(const std::vector<Shoebox> &boxes);

} // namespace mxi
