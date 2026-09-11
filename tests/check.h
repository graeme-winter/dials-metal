// A test harness in a hundred lines, rather than a dependency.
//
// Registration is by static constructor, so a test file needs no list to be
// added to. Forgetting to register is the commonest way a test silently stops
// running, and this makes it impossible.
#pragma once

#include <cmath>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace check {

struct Case {
  std::string name;
  std::function<void()> body;
};

inline std::vector<Case> &registry() {
  static std::vector<Case> cases;
  return cases;
}

struct Failure {
  std::string what;
};

inline void fail(const std::string &message) { throw Failure{message}; }

inline void is_true(bool condition, const std::string &message) {
  if (!condition) fail(message);
}

inline void close(double a, double b, double tolerance,
                  const std::string &message) {
  if (!(std::abs(a - b) <= tolerance)) {
    char buffer[512];
    std::snprintf(buffer, sizeof(buffer), "%s: %.17g vs %.17g (tol %.3g)",
                  message.c_str(), a, b, tolerance);
    fail(buffer);
  }
}

inline void equal(long long a, long long b, const std::string &message) {
  if (a != b) {
    char buffer[512];
    std::snprintf(buffer, sizeof(buffer), "%s: %lld vs %lld", message.c_str(),
                  a, b);
    fail(buffer);
  }
}

struct Register {
  Register(const char *name, std::function<void()> body) {
    registry().push_back({name, std::move(body)});
  }
};

inline int run_all() {
  int failed = 0;
  for (const Case &c : registry()) {
    try {
      c.body();
      std::printf("  ok    %s\n", c.name.c_str());
    } catch (const Failure &f) {
      std::printf("  FAIL  %s\n        %s\n", c.name.c_str(), f.what.c_str());
      ++failed;
    } catch (const std::exception &e) {
      std::printf("  ERROR %s\n        %s\n", c.name.c_str(), e.what());
      ++failed;
    }
  }
  std::printf("%zu tests, %d failed\n", registry().size(), failed);
  return failed == 0 ? 0 : 1;
}

}  // namespace check

#define TEST(name)                                       \
  static void name();                                    \
  static ::check::Register register_##name(#name, name); \
  static void name()
