// mxi_background: fit backgrounds to pixel values read from a file.
//
// One shoebox per line: the pixel counts, whitespace separated. Prints one
// fitted background per line. It exists so the C++ the pipeline will use can
// be driven against a DIALS integrated.refl without a pixel reader in C++
// yet -- the comparison is of this code, not of a transcription of it.

#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "args.h"
#include "background.h"

using namespace mxi;

int main(int argc, char **argv) {
  const std::set<std::string> known = {"--tuning"};
  const Arguments args = parse_arguments(argc, argv, known, known);
  if (!args.ok) {
    std::fprintf(stderr, "mxi_background: %s\n", args.error.c_str());
    return 2;
  }
  BackgroundOptions options;
  options.tuning = args.number("--tuning", 1.345);

  std::istream *in = &std::cin;
  std::ifstream file;
  if (!args.positional.empty()) {
    file.open(args.positional[0]);
    if (!file) {
      std::fprintf(stderr, "mxi_background: cannot open %s\n",
                   args.positional[0].c_str());
      return 1;
    }
    in = &file;
  }

  std::string line;
  while (std::getline(*in, line)) {
    std::istringstream row(line);
    std::vector<double> values;
    double v = 0.0;
    while (row >> v) values.push_back(v);
    const BackgroundResult r = glm_background(values, options);
    std::printf("%.10g %d %d\n", r.mean, r.iterations, r.converged ? 1 : 0);
  }
  return 0;
}
