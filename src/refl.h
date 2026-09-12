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
  Column &int_column(const std::string &name, const std::string &type,
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
