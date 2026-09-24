#include "shoebox.h"

#include <cstring>

namespace mxi {

namespace {

template <typename T>
T take(const std::string &blob, std::size_t *at) {
  if (*at + sizeof(T) > blob.size()) {
    throw ReflError("the shoebox column ends in the middle of a record");
  }
  T value;
  std::memcpy(&value, blob.data() + *at, sizeof(T));
  *at += sizeof(T);
  return value;
}

template <typename T>
void put(std::string *out, T value) {
  const auto *bytes = reinterpret_cast<const char *>(&value);
  out->append(bytes, sizeof(T));
}

}  // namespace

std::vector<Shoebox> decode_shoeboxes(const Table &table) {
  std::vector<Shoebox> out;
  auto it = table.opaque().find("shoebox");
  if (it == table.opaque().end()) return out;
  const std::string &blob = it->second.bytes;

  std::size_t at = 0;
  out.reserve(table.nrows);
  for (std::size_t i = 0; i < table.nrows; ++i) {
    Shoebox box;
    box.panel = take<std::int32_t>(blob, &at);
    for (std::int32_t &v : box.bbox) v = take<std::int32_t>(blob, &at);
    box.flag = take<std::uint8_t>(blob, &at);
    if (box.nx() < 0 || box.ny() < 0 || box.nz() < 0) {
      throw ReflError("a shoebox has a bounding box with a negative extent");
    }
    const std::size_t n = box.size();
    if (at + 9 * n > blob.size()) {
      throw ReflError("the shoebox column ends in the middle of a record");
    }
    box.data.resize(n);
    std::memcpy(box.data.data(), blob.data() + at, 4 * n);
    at += 4 * n;
    box.mask.resize(n);
    std::memcpy(box.mask.data(), blob.data() + at, n);
    at += n;
    box.background.resize(n);
    std::memcpy(box.background.data(), blob.data() + at, 4 * n);
    at += 4 * n;
    out.push_back(std::move(box));
  }
  // Exactly, not nearly. A layout that leaves bytes over is wrong somewhere
  // earlier, and the boxes already read would be wrong with it.
  if (at != blob.size()) {
    throw ReflError("the shoebox column has " + std::to_string(blob.size() - at) +
                    " bytes left over after " + std::to_string(table.nrows) +
                    " records; the layout does not match this file");
  }
  return out;
}

void to_dials_convention(Shoebox *box) {
  if (box == nullptr) return;
  for (std::uint8_t &m : box->mask) {
    if ((m & shoebox_mask::kValid) == 0) m = 0;
  }
}

std::string encode_shoeboxes(const std::vector<Shoebox> &boxes) {
  std::string out;
  for (const Shoebox &box : boxes) {
    put(&out, box.panel);
    for (std::int32_t v : box.bbox) put(&out, v);
    put(&out, box.flag);
    const std::size_t n = box.size();
    if (box.data.size() != n || box.mask.size() != n || box.background.size() != n) {
      throw ReflError("a shoebox's arrays do not match its bounding box");
    }
    out.append(reinterpret_cast<const char *>(box.data.data()), 4 * n);
    out.append(reinterpret_cast<const char *>(box.mask.data()), n);
    out.append(reinterpret_cast<const char *>(box.background.data()), 4 * n);
  }
  return out;
}

}  // namespace mxi
