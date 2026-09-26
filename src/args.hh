// Splitting a command line into positional arguments and options.
//
// Written because all three programs originally read their file names from
// argv[1] and argv[2] and looked for options from argv[3] onwards, so
//
//     mxi_refine --beam indexed.expt indexed.refl
//
// took "--beam" as the experiment list, "indexed.expt" as the reflections, and
// then complained that "indexed.refl" was an unknown option. Options before
// positionals is what every other command line tool accepts, and the failure
// was confusing rather than merely inconvenient: the error named the wrong
// argument entirely.
//
// It lives here rather than in apps/ so that it can be tested. Argument
// parsing is exactly the kind of code that never gets tested because it feels
// too simple to get wrong.

#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

namespace mxi {

struct Arguments {
  std::vector<std::string> positional;
  //: Flag to value. A flag that takes no value maps to the empty string.
  std::map<std::string, std::string> options;
  bool ok = true;
  std::string error;
  bool help = false;

  bool has(const std::string &flag) const { return options.count(flag) > 0; }
  std::string value(const std::string &flag, const std::string &fallback = "") const {
    auto it = options.find(flag);
    return it == options.end() ? fallback : it->second;
  }
  double number(const std::string &flag, double fallback) const;
};

// `takes_value` names the flags that consume the following argument. Anything
// beginning with a dash and not in `known` is an error, which is deliberate:
// silently ignoring a misspelled option is how a run comes to use settings
// nobody chose.
Arguments parse_arguments(int argc, char **argv, const std::set<std::string> &known,
                          const std::set<std::string> &takes_value);

}  // namespace mxi
