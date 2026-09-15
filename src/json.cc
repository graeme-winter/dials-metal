#include "json.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace mxi {
namespace json {

namespace {
const Value &null_value() {
  static const Value null;
  return null;
}
}  // namespace

const Value &Value::operator[](const std::string &key) const {
  if (type_ != Type::Object) return null_value();
  auto it = object_.find(key);
  return it == object_.end() ? null_value() : it->second;
}

const Value &Value::operator[](std::size_t i) const {
  if (type_ != Type::Array || i >= array_.size()) return null_value();
  return array_[i];
}

std::vector<double> Value::numbers() const {
  std::vector<double> out;
  if (type_ != Type::Array) return out;
  out.reserve(array_.size());
  for (const Value &v : array_) out.push_back(v.as_number());
  return out;
}

bool Value::numbers(std::size_t expected, double *out) const {
  if (type_ != Type::Array || array_.size() != expected) return false;
  for (std::size_t i = 0; i < expected; ++i) {
    if (!array_[i].is_number()) return false;
    out[i] = array_[i].as_number();
  }
  return true;
}

// --------------------------------------------------------------------------
// parsing
// --------------------------------------------------------------------------

namespace {

class Parser {
 public:
  explicit Parser(const std::string &text) : text_(text) {}

  Value parse() {
    skip_space();
    Value v = value();
    skip_space();
    if (at_ != text_.size()) fail("trailing content after the document");
    return v;
  }

 private:
  const std::string &text_;
  std::size_t at_ = 0;

  [[noreturn]] void fail(const std::string &what) { throw ParseError(what, at_); }

  char peek() const { return at_ < text_.size() ? text_[at_] : '\0'; }

  void skip_space() {
    while (at_ < text_.size()) {
      const char c = text_[at_];
      if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
        ++at_;
      } else {
        break;
      }
    }
  }

  void expect(char c) {
    if (peek() != c) fail(std::string("expected '") + c + "'");
    ++at_;
  }

  bool literal(const char *word) {
    const std::size_t n = std::char_traits<char>::length(word);
    if (text_.compare(at_, n, word) != 0) return false;
    at_ += n;
    return true;
  }

  Value value() {
    switch (peek()) {
      case '{':
        return object();
      case '[':
        return array();
      case '"':
        return Value(string());
      case 't':
        if (literal("true")) return Value(true);
        fail("expected true");
      case 'f':
        if (literal("false")) return Value(false);
        fail("expected false");
      case 'n':
        if (literal("null")) return Value();
        fail("expected null");
      // A bare Infinity or NaN is not JSON, but it does appear in files
      // written by software that let a division through, and refusing to read
      // such a file is less useful than reading it and letting the value be
      // obviously wrong further on.
      case 'N':
        if (literal("NaN")) return Value(std::nan(""));
        fail("expected NaN");
      case 'I':
        if (literal("Infinity")) return Value(HUGE_VAL);
        fail("expected Infinity");
      default:
        return number();
    }
  }

  Value object() {
    expect('{');
    Object out;
    skip_space();
    if (peek() == '}') {
      ++at_;
      return Value(std::move(out));
    }
    for (;;) {
      skip_space();
      std::string key = string();
      skip_space();
      expect(':');
      skip_space();
      out.emplace(std::move(key), value());
      skip_space();
      if (peek() == ',') {
        ++at_;
        continue;
      }
      expect('}');
      return Value(std::move(out));
    }
  }

  Value array() {
    expect('[');
    Array out;
    skip_space();
    if (peek() == ']') {
      ++at_;
      return Value(std::move(out));
    }
    for (;;) {
      skip_space();
      out.push_back(value());
      skip_space();
      if (peek() == ',') {
        ++at_;
        continue;
      }
      expect(']');
      return Value(std::move(out));
    }
  }

  std::string string() {
    expect('"');
    std::string out;
    while (at_ < text_.size()) {
      const char c = text_[at_++];
      if (c == '"') return out;
      if (c != '\\') {
        out.push_back(c);
        continue;
      }
      if (at_ >= text_.size()) fail("unterminated escape");
      const char e = text_[at_++];
      switch (e) {
        case '"': out.push_back('"'); break;
        case '\\': out.push_back('\\'); break;
        case '/': out.push_back('/'); break;
        case 'b': out.push_back('\b'); break;
        case 'f': out.push_back('\f'); break;
        case 'n': out.push_back('\n'); break;
        case 'r': out.push_back('\r'); break;
        case 't': out.push_back('\t'); break;
        case 'u': {
          if (at_ + 4 > text_.size()) fail("truncated \\u escape");
          unsigned code = 0;
          for (int i = 0; i < 4; ++i) {
            const char h = text_[at_++];
            code *= 16;
            if (h >= '0' && h <= '9') code += static_cast<unsigned>(h - '0');
            else if (h >= 'a' && h <= 'f') code += static_cast<unsigned>(h - 'a' + 10);
            else if (h >= 'A' && h <= 'F') code += static_cast<unsigned>(h - 'A' + 10);
            else fail("bad hex digit in \\u escape");
          }
          // UTF-8 encode. Surrogate pairs are not joined: no .expt field has
          // ever needed one, and a lone surrogate round trips as itself.
          if (code < 0x80) {
            out.push_back(static_cast<char>(code));
          } else if (code < 0x800) {
            out.push_back(static_cast<char>(0xC0 | (code >> 6)));
            out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
          } else {
            out.push_back(static_cast<char>(0xE0 | (code >> 12)));
            out.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
            out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
          }
          break;
        }
        default:
          fail("unknown escape");
      }
    }
    fail("unterminated string");
  }

  Value number() {
    const std::size_t start = at_;
    if (peek() == '-' || peek() == '+') ++at_;
    if (literal("Infinity")) {
      return Value(text_[start] == '-' ? -HUGE_VAL : HUGE_VAL);
    }
    while (at_ < text_.size()) {
      const char c = text_[at_];
      if ((c >= '0' && c <= '9') || c == '.' || c == 'e' || c == 'E' ||
          c == '+' || c == '-') {
        ++at_;
      } else {
        break;
      }
    }
    if (at_ == start) fail("expected a value");
    try {
      return Value(std::stod(text_.substr(start, at_ - start)));
    } catch (const std::exception &) {
      at_ = start;
      fail("malformed number");
    }
  }
};

void write(std::ostringstream &os, const Value &v, int indent, int depth) {
  const auto newline = [&](int d) {
    if (indent > 0) {
      os << "\n";
      os << std::string(static_cast<std::size_t>(indent * d), ' ');
    }
  };

  switch (v.type()) {
    case Value::Type::Null:
      os << "null";
      break;
    case Value::Type::Bool:
      os << (v.as_bool() ? "true" : "false");
      break;
    case Value::Type::Number: {
      const double d = v.as_number();
      if (std::isnan(d)) {
        os << "NaN";
      } else if (std::isinf(d)) {
        os << (d > 0 ? "Infinity" : "-Infinity");
      } else if (v.is_integral() && std::abs(d) < 1e15) {
        // Written without a decimal point only when the caller built it from
        // an integer type. Deciding by whether the value happens to be whole
        // would turn an oscillation that starts at zero into `0`, and dxtbx
        // reads the type of an array from its first element.
        os << static_cast<long long>(d);
      } else {
        // Seventeen significant digits recovers any double exactly, which is
        // what keeps a geometry from drifting when it passes through. A whole
        // value must still carry a decimal point, or it reads back as an
        // integer.
        char buffer[40];
        std::snprintf(buffer, sizeof(buffer), "%.17g", d);
        std::string text(buffer);
        if (text.find_first_of(".eEnN") == std::string::npos) text += ".0";
        os << text;
      }
      break;
    }
    case Value::Type::String: {
      os << '"';
      for (char c : v.as_string()) {
        switch (c) {
          case '"': os << "\\\""; break;
          case '\\': os << "\\\\"; break;
          case '\n': os << "\\n"; break;
          case '\r': os << "\\r"; break;
          case '\t': os << "\\t"; break;
          default:
            if (static_cast<unsigned char>(c) < 0x20) {
              char buffer[8];
              std::snprintf(buffer, sizeof(buffer), "\\u%04x", c);
              os << buffer;
            } else {
              os << c;
            }
        }
      }
      os << '"';
      break;
    }
    case Value::Type::Array: {
      const Array &a = v.as_array();
      if (a.empty()) {
        os << "[]";
        break;
      }
      os << '[';
      for (std::size_t i = 0; i < a.size(); ++i) {
        if (i) os << ',';
        newline(depth + 1);
        write(os, a[i], indent, depth + 1);
      }
      newline(depth);
      os << ']';
      break;
    }
    case Value::Type::Object: {
      const Object &o = v.as_object();
      if (o.empty()) {
        os << "{}";
        break;
      }
      os << '{';
      bool first = true;
      for (const auto &entry : o) {
        if (!first) os << ',';
        first = false;
        newline(depth + 1);
        write(os, Value(entry.first), indent, depth + 1);
        os << ": ";
        write(os, entry.second, indent, depth + 1);
      }
      newline(depth);
      os << '}';
      break;
    }
  }
}

}  // namespace

Value parse(const std::string &text) { return Parser(text).parse(); }

Value parse_file(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return parse(buffer.str());
}

std::string dump(const Value &value, int indent) {
  std::ostringstream os;
  write(os, value, indent, 0);
  return os.str();
}

void dump_file(const std::string &path, const Value &value, int indent) {
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("cannot write " + path);
  out << dump(value, indent) << "\n";
}

}  // namespace json
}  // namespace mxi
