#pragma once

// Spots per image as a text histogram, drawn as dials.find_spots draws it:
// the one part of a spot finder's report a user reads at a glance, since a dip
// or a dead stretch in the scan shows at once. Header-only and free of the
// spot finder's types, so it can be tested without linking the spot finder.

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

namespace spots {

//: The histogram's lines, top first. per_image[k] is the number of spots on
//: image first_image + k. Each of up to `width` columns is a run of WHOLE
//: images and holds the spots found in them, so the axis labels are exact.
//:
//: dials.find_spots bins differently -- over the data's own z range, with every
//: z moved by a quarter of an image, up or down alternately by its position in
//: the table, to break ties at the bin edges -- so on the same spots the
//: shapes agree and the maximum a bin differs: 269 here against DIALS' 286 on
//: a 300 image insulin sweep. The column heights are rounded as DIALS rounds
//: them.
inline std::vector<std::string>
spot_histogram(const std::vector<std::size_t> &per_image, long long first_image,
               std::size_t width = 60, std::size_t height = 10) {
  std::vector<std::string> out;
  const std::size_t images = per_image.size();
  if (images == 0 || width == 0 || height == 0)
    return out;
  const std::size_t columns = std::min(width, images);
  std::vector<std::size_t> bins(columns, 0);
  std::size_t total = 0;
  for (std::size_t k = 0; k < images; ++k) {
    // Images split among columns as evenly as integers allow.
    bins[k * columns / images] += per_image[k];
    total += per_image[k];
  }
  const std::size_t most = *std::max_element(bins.begin(), bins.end());
  out.push_back(std::to_string(total) + " spots found on " +
                std::to_string(images) + " images (max " +
                std::to_string(most) + " / bin)");
  // A column's height is its count scaled to the tallest and rounded, as
  // dials.find_spots draws it, so a column holding a twentieth of the maximum
  // or less shows no star at all.
  std::vector<std::size_t> heights(columns, 0);
  for (std::size_t c = 0; c < columns; ++c) {
    if (most > 0)
      heights[c] = (2 * bins[c] * height + most) / (2 * most);
  }
  for (std::size_t row = height; row >= 1; --row) {
    std::string line(columns, ' ');
    for (std::size_t c = 0; c < columns; ++c) {
      if (heights[c] >= row)
        line[c] = '*';
    }
    out.push_back(line);
  }
  // The axis: first image at the left, last at the right, "image" between.
  const std::string left = std::to_string(first_image);
  const std::string right =
      std::to_string(first_image + static_cast<long long>(images) - 1);
  std::string axis(columns, ' ');
  if (left.size() + right.size() + 1 <= columns) {
    axis.replace(0, left.size(), left);
    axis.replace(columns - right.size(), right.size(), right);
    const std::string label = "image";
    if (left.size() + right.size() + label.size() + 2 <= columns)
      axis.replace((columns - label.size()) / 2, label.size(), label);
  } else {
    axis = left + " to " + right;
  }
  out.push_back(axis);
  return out;
}

} // namespace spots
