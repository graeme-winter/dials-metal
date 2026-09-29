#include <cstdio>
#include <thread>
#include <vector>

#include "../src/timing.hh"
#include "check.hh"

namespace mxi {

TEST(thread_seconds_sum_what_every_thread_adds) {
  ThreadSeconds total;
  std::vector<std::thread> threads;
  for (int t = 0; t < 8; ++t)
    threads.emplace_back([&] {
      for (int k = 0; k < 1000; ++k)
        total.add(0.001);
    });
  for (std::thread &t : threads)
    t.join();
  check::close(total.seconds(), 8.0, 1e-6, "8 threads x 1000 x 1 ms");
}

TEST(timing_prints_nothing_unless_asked) {
  std::FILE *f = std::tmpfile();
  Timing off(false);
  off.add("a phase", 1.0);
  off.report(f);
  check::equal(static_cast<long long>(std::ftell(f)), 0,
               "nothing without --timing");
  Timing on(true);
  on.add("a phase", 1.0);
  on.add("its part", 0.5, 1);
  on.report(f);
  std::rewind(f);
  char buffer[512] = {0};
  const std::size_t n = std::fread(buffer, 1, sizeof buffer - 1, f);
  std::fclose(f);
  const std::string text(buffer, n);
  check::is_true(text.find("\nTiming\n") != std::string::npos,
                 "a table with --timing");
  check::is_true(text.find("    its part") != std::string::npos,
                 "the part indented under it");
  check::is_true(text.find("  total") != std::string::npos,
                 "and the run's total");
}

} // namespace mxi
