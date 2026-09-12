// Command line parsing, which broke in a way that named the wrong argument.
//
//     mxi_refine --beam indexed.expt indexed.refl
//
// read argv[1] and argv[2] as the two file names and looked for options from
// argv[3], so it took "--beam" as the experiment list and then reported
// "unknown option 'indexed.refl'". Options before positionals is what every
// other command line tool accepts.

#include <set>
#include <string>
#include <vector>

#include "../src/args.h"
#include "check.h"

using namespace mxi;

namespace {

Arguments run(const std::vector<std::string> &words) {
  std::vector<char *> argv;
  std::vector<std::string> owned = words;
  argv.push_back(const_cast<char *>("prog"));
  for (std::string &w : owned) argv.push_back(const_cast<char *>(w.c_str()));
  static const std::set<std::string> known = {"--beam", "--strong-only",
                                              "--outlier-sigma", "--output-expt"};
  static const std::set<std::string> takes_value = {"--outlier-sigma",
                                                    "--output-expt"};
  return parse_arguments(static_cast<int>(argv.size()), argv.data(), known,
                         takes_value);
}

}  // namespace

TEST(options_may_come_before_the_file_names) {
  const Arguments a = run({"--beam", "indexed.expt", "indexed.refl"});
  check::is_true(a.ok, a.error.empty() ? "parsed" : a.error);
  check::equal(static_cast<long long>(a.positional.size()), 2, "two files");
  check::is_true(a.positional[0] == "indexed.expt", "first file");
  check::is_true(a.positional[1] == "indexed.refl", "second file");
  check::is_true(a.has("--beam"), "the flag was seen");
}

TEST(options_may_come_after_the_file_names) {
  const Arguments a = run({"indexed.expt", "indexed.refl", "--strong-only"});
  check::is_true(a.ok, "parsed");
  check::equal(static_cast<long long>(a.positional.size()), 2, "two files");
  check::is_true(a.has("--strong-only"), "the flag was seen");
}

TEST(options_may_be_interleaved) {
  const Arguments a =
      run({"--outlier-sigma", "3", "a.expt", "--beam", "b.refl"});
  check::is_true(a.ok, "parsed");
  check::is_true(a.positional[0] == "a.expt" && a.positional[1] == "b.refl",
                 "files found either side of the options");
  check::close(a.number("--outlier-sigma", 4.0), 3.0, 1e-12, "value taken");
  check::is_true(a.has("--beam"), "flag after a file");
}

TEST(a_value_is_not_mistaken_for_a_file) {
  const Arguments a = run({"--outlier-sigma", "3", "a.expt", "b.refl"});
  check::equal(static_cast<long long>(a.positional.size()), 2,
               "the 3 belongs to the option");
}

TEST(an_unknown_option_is_refused_rather_than_ignored) {
  const Arguments a = run({"a.expt", "b.refl", "--stromg-only"});
  check::is_true(!a.ok, "must be refused");
  check::is_true(a.error.find("--stromg-only") != std::string::npos,
                 "and must name the option, not a file");
}

TEST(a_missing_value_is_reported) {
  const Arguments a = run({"a.expt", "b.refl", "--outlier-sigma"});
  check::is_true(!a.ok, "must be refused");
  check::is_true(a.error.find("needs a value") != std::string::npos, "says why");
}

TEST(help_is_recognised_anywhere) {
  check::is_true(run({"--help"}).help, "alone");
  check::is_true(run({"a.expt", "--help", "b.refl"}).help, "among files");
  check::is_true(run({"-h"}).help, "short form");
}

TEST(missing_files_are_counted_not_guessed) {
  const Arguments a = run({"--beam"});
  check::is_true(a.ok, "the option itself is fine");
  check::equal(static_cast<long long>(a.positional.size()), 0, "no files");
}

TEST(defaults_are_returned_for_absent_options) {
  const Arguments a = run({"a.expt", "b.refl"});
  check::close(a.number("--outlier-sigma", 4.0), 4.0, 1e-12, "default kept");
  check::is_true(a.value("--output-expt", "refined.expt") == "refined.expt",
                 "string default kept");
}
