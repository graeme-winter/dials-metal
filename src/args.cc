#include "args.hh"

#include <cstdlib>

namespace mxi {

double Arguments::number(const std::string &flag, double fallback) const {
  auto it = options.find(flag);
  if (it == options.end() || it->second.empty())
    return fallback;
  char *end = nullptr;
  const double v = std::strtod(it->second.c_str(), &end);
  return end == it->second.c_str() ? fallback : v;
}

Arguments parse_arguments(int argc, char **argv,
                          const std::set<std::string> &known,
                          const std::set<std::string> &takes_value,
                          const std::set<std::string> &optional_value) {
  Arguments out;
  const auto whole_number = [](const std::string &s) {
    if (s.empty())
      return false;
    for (char ch : s)
      if (ch < '0' || ch > '9')
        return false;
    return true;
  };
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      out.help = true;
      return out;
    }
    // A bare "-" is a filename by convention, not an option.
    if (arg.size() > 1 && arg[0] == '-') {
      if (!known.count(arg)) {
        out.ok = false;
        out.error = "unknown option '" + arg + "'";
        return out;
      }
      if (takes_value.count(arg)) {
        if (i + 1 >= argc) {
          out.ok = false;
          out.error = "option '" + arg + "' needs a value";
          return out;
        }
        out.options[arg] = argv[++i];
      } else if (optional_value.count(arg)) {
        if (i + 1 < argc && whole_number(argv[i + 1]))
          out.options[arg] = argv[++i];
        else
          out.options[arg] = "";
      } else {
        out.options[arg] = "";
      }
      continue;
    }
    out.positional.push_back(arg);
  }
  return out;
}

} // namespace mxi
