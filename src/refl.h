// Reading and writing DIALS reflection tables, hand-rolled msgpack.
//
// The layout, confirmed against files written by DIALS 3.x and by
// dials-metal-find-spots:
//
//   [ "dials::af::reflection_table", 2, { identifiers, nrows, data } ]
//   data[name] = [ type_name, [ nrows, blob ] ]
//
// Text is msgpack `str`, payloads are `bin`, blobs are little-endian with the
// components of a compound type adjacent. The count inside a payload is the
// number of ROWS, not scalars, and is checked rather than trusted, because the
// two disagreeing is what a truncated file looks like.
//
// Element widths: int 4, std::size_t 8, double 8, bool 1, vec2<double> 16,
// vec3<double> 24, int6 24, cctbx::miller::index<> 12.
//
// Columns of a type this does not understand -- a shoebox, most notably -- are
// read as opaque bytes and are droppable. They are not carried through on
// write, because nothing here can subset one correctly and a shoebox column
// that silently stops matching its table would be worse than its absence.

#pragma once

#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace mxi {

// The `flags` column is a bitmask, and dials.* filters on it. A table whose
// flags are never set processes perfectly and then cannot be selected,
// filtered or plotted by anything downstream, because every DIALS tool asks
// the flags which reflections it is looking at rather than inspecting the
// Miller indices.
//
// These values are dxtbx's, and the two that are used here were confirmed
// against a real DIALS indexed.refl: its indexed reflections carry 36, which
// is strong | indexed.
namespace flag {
constexpr std::int64_t kPredicted = 1 << 0;
constexpr std::int64_t kObserved = 1 << 1;
constexpr std::int64_t kIndexed = 1 << 2;
constexpr std::int64_t kUsedInRefinement = 1 << 3;
constexpr std::int64_t kStrong = 1 << 5;
//: Set on reflections the refinement drew in and then rejected. Identified
//: from a real DIALS indexed.refl rather than guessed: 8230 rows carry bit 17,
//: every one of them indexed, none of them also marked used_in_refinement, and
//: their median |xyzcal - xyzobs| is 0.843 px against 0.310 for the rest.
constexpr std::int64_t kCentroidOutlier = 1 << 17;
}  // namespace flag


class ReflError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

struct Column {
  std::string type;       // the C++ type name as it appears in the file
  std::size_t width = 1;  // components per row
  bool integral = false;  // which of the two stores below is in use
  std::vector<double> reals;
  std::vector<std::int64_t> ints;

  std::size_t rows() const {
    return (integral ? ints.size() : reals.size()) / (width ? width : 1);
  }
  double real(std::size_t row, std::size_t k = 0) const {
    return reals[row * width + k];
  }
  std::int64_t integer(std::size_t row, std::size_t k = 0) const {
    return ints[row * width + k];
  }
};

class Table {
 public:
  std::size_t nrows = 0;
  std::map<std::size_t, std::string> identifiers;
  int version = 2;

  bool has(const std::string &name) const { return columns_.count(name) > 0; }
  const Column &at(const std::string &name) const;
  std::vector<std::string> names() const;
  const std::vector<std::string> &dropped() const { return dropped_; }

  Column &real_column(const std::string &name, const std::string &type,
                      std::size_t width);
  // NOTE: these REPLACE an existing column of the same name with a zeroed one.
  // That is what a derived column wants -- xyzcal.px and the rest are
  // recomputed in full, and leaving stale values in the rows that are skipped
  // would be worse than clearing them. It is wrong for any column that has to
  // be read before it is written, and `flags` is one: setting the indexed bit
  // this way silently threw away the strong bit that dials.find_spots had set,
  // and the loss is invisible until something downstream filters on it.
  Column &int_column(const std::string &name, const std::string &type,
                     std::size_t width);

  // Returns the existing column if there is one, so it can be modified rather
  // than replaced. Creates a zeroed column if not.
  Column &modify_int_column(const std::string &name, const std::string &type,
                            std::size_t width);
  void set(const std::string &name, Column column) {
    columns_[name] = std::move(column);
  }

  void validate() const;

 private:
  std::map<std::string, Column> columns_;
  std::vector<std::string> dropped_;
  friend Table read_reflections(const std::string &);
};

Table read_reflections(const std::string &path);
void write_reflections(const std::string &path, const Table &table);

}  // namespace mxi
