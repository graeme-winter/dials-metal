// A JSON reader and writer, in about three hundred lines.
//
// Hand-rolled for the same reason as everything else here: no dependencies.
// An .expt file is JSON of a very tame shape -- objects, arrays, numbers,
// strings -- and a parser for that is small enough that pulling in a library
// would cost more in build friction than it saves.
//
// Numbers are always double. The one thing that has to be right is round
// tripping: a detector origin read and written back must be the same bits, or
// a pipeline that passes an experiment through unchanged will appear to have
// moved the detector. Writing uses seventeen significant digits, which is
// enough to recover any double exactly.

#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace mxi {
namespace json {

class Value;
using Object = std::map<std::string, Value>;
using Array = std::vector<Value>;

class ParseError : public std::runtime_error {
public:
  ParseError(const std::string &what, std::size_t offset)
      : std::runtime_error(what + " at offset " + std::to_string(offset)),
        offset(offset) {}
  std::size_t offset;
};

class Value {
public:
  enum class Type { Null, Bool, Number, String, Array, Object };

  Value() = default;
  Value(bool b) : type_(Type::Bool), bool_(b) {}
  // Whether a number was written as an integer or as a float is not a
  // formatting detail here: dxtbx inspects the type of the first element of
  // some arrays and refuses the file if it is not what it expects. A scan's
  // oscillation array whose first value happens to be zero, written as `0`
  // rather than `0.0`, fails dials.refine with
  //
  //     DXTBX_ASSERT(obj_type == "float") failure
  //
  // and nothing in this package could have caught it, because it round-trips
  // through this reader perfectly. So the distinction is carried on the value
  // rather than guessed from its magnitude at the point of writing.
  Value(double d) : type_(Type::Number), number_(d) {}
  Value(int i) : type_(Type::Number), number_(i), integral_(true) {}
  Value(long i)
      : type_(Type::Number), number_(static_cast<double>(i)), integral_(true) {}
  Value(long long i)
      : type_(Type::Number), number_(static_cast<double>(i)), integral_(true) {}
  Value(std::size_t i)
      : type_(Type::Number), number_(static_cast<double>(i)), integral_(true) {}

  //: True when this number was built from an integer type and should be
  //: written without a decimal point.
  bool is_integral() const { return type_ == Type::Number && integral_; }
  Value(const char *s) : type_(Type::String), string_(s) {}
  Value(std::string s) : type_(Type::String), string_(std::move(s)) {}
  Value(Array a) : type_(Type::Array), array_(std::move(a)) {}
  Value(Object o) : type_(Type::Object), object_(std::move(o)) {}

  Type type() const { return type_; }
  bool is_null() const { return type_ == Type::Null; }
  bool is_number() const { return type_ == Type::Number; }
  bool is_string() const { return type_ == Type::String; }
  bool is_array() const { return type_ == Type::Array; }
  bool is_object() const { return type_ == Type::Object; }

  bool as_bool(bool fallback = false) const {
    return type_ == Type::Bool ? bool_ : fallback;
  }
  double as_number(double fallback = 0.0) const {
    return type_ == Type::Number ? number_ : fallback;
  }
  const std::string &as_string() const { return string_; }
  const Array &as_array() const { return array_; }
  const Object &as_object() const { return object_; }
  Array &as_array() { return array_; }
  Object &as_object() { return object_; }

  // Member lookup that returns a null Value rather than throwing. A missing
  // optional model in an .expt is normal, and forcing every caller to check
  // first turns the common path into noise.
  const Value &operator[](const std::string &key) const;
  const Value &operator[](std::size_t i) const;
  bool contains(const std::string &key) const {
    return type_ == Type::Object && object_.count(key) > 0;
  }
  std::size_t size() const {
    return type_ == Type::Array ? array_.size()
                                : (type_ == Type::Object ? object_.size() : 0);
  }

  // Numeric array helpers, which is what nearly every .expt field is.
  std::vector<double> numbers() const;
  bool numbers(std::size_t expected, double *out) const;

private:
  Type type_ = Type::Null;
  bool bool_ = false;
  double number_ = 0.0;
  bool integral_ = false;
  std::string string_;
  Array array_;
  Object object_;
};

Value parse(const std::string &text);
Value parse_file(const std::string &path);

// `indent` of zero writes it on one line. DIALS writes .expt files indented,
// and matching that makes a diff between one of ours and one of theirs
// readable, which is worth the bytes.
std::string dump(const Value &value, int indent = 1);
void dump_file(const std::string &path, const Value &value, int indent = 1);

} // namespace json
} // namespace mxi
