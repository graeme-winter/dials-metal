// Comparing two signal pixel lists, and saying by how much they differ.
//
// The device test wants a verdict; the benchmark wants to carry on and report
// timings even when the two implementations have parted company. Both want the
// same measurement, so it lives here rather than twice.
//
// Both lists are ascending by index -- the CPU emits in scan order, the device
// sorts before returning -- so this is a merge rather than a search, and a
// pixel found by one and not the other is distinguishable from one found by
// both with a different background. Which of those it is matters: a different
// background is arithmetic, a different set of pixels is a different answer.
//
// Written for -DENABLE_FAST_MATH, where the question stops being "are they
// identical" and becomes "how far apart are they".

#ifndef SPOTFINDER_COMPARE_SIGNAL_HH
#define SPOTFINDER_COMPARE_SIGNAL_HH

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

#include "signal_pixel.hh"

namespace compare_signal {

struct Difference {
  std::size_t cpu_pixels = 0;
  std::size_t gpu_pixels = 0;

  std::size_t only_cpu = 0; // found by the CPU and not the device
  std::size_t only_gpu = 0;
  std::size_t common = 0;

  // Of the pixels both found: the integer fields must agree, since a pixel's
  // value is a copy and its population is a count. A difference there is a
  // fault, not rounding.
  std::size_t integer_differs = 0;
  std::size_t background_differs = 0;

  double worst_absolute = 0.0;
  double worst_relative = 0.0;

  // The first pixel that differs in any way, for the detail line.
  bool have_example = false;
  SignalPixel cpu_example{};
  SignalPixel gpu_example{};

  bool identical() const {
    return only_cpu == 0 && only_gpu == 0 && integer_differs == 0 &&
           background_differs == 0;
  }

  // The same pixels, differing only in the float field. Under fast math this is
  // the tolerable outcome; without it, it is still a failure.
  bool same_pixels() const {
    return only_cpu == 0 && only_gpu == 0 && integer_differs == 0;
  }
};

inline Difference difference(const std::vector<SignalPixel> &cpu,
                             const std::vector<SignalPixel> &gpu) {
  Difference result;
  result.cpu_pixels = cpu.size();
  result.gpu_pixels = gpu.size();

  std::size_t a = 0, b = 0;
  while (a < cpu.size() || b < gpu.size()) {
    const bool take_cpu =
        b >= gpu.size() || (a < cpu.size() && cpu[a].index < gpu[b].index);
    const bool take_gpu =
        a >= cpu.size() || (b < gpu.size() && gpu[b].index < cpu[a].index);

    if (take_cpu) {
      result.only_cpu++;
      if (!result.have_example) {
        result.have_example = true;
        result.cpu_example = cpu[a];
      }
      a++;
      continue;
    }
    if (take_gpu) {
      result.only_gpu++;
      if (!result.have_example) {
        result.have_example = true;
        result.gpu_example = gpu[b];
      }
      b++;
      continue;
    }

    result.common++;
    const SignalPixel &left = cpu[a];
    const SignalPixel &right = gpu[b];
    const bool integers =
        left.value != right.value || left.population != right.population;
    const bool background = left.background != right.background;
    if (integers)
      result.integer_differs++;
    if (background) {
      result.background_differs++;
      const double absolute = std::fabs(static_cast<double>(left.background) -
                                        static_cast<double>(right.background));
      const double scale = std::fabs(static_cast<double>(left.background));
      if (absolute > result.worst_absolute)
        result.worst_absolute = absolute;
      if (scale > 0.0 && absolute / scale > result.worst_relative)
        result.worst_relative = absolute / scale;
    }
    if ((integers || background) && !result.have_example) {
      result.have_example = true;
      result.cpu_example = left;
      result.gpu_example = right;
    }
    a++;
    b++;
  }
  return result;
}

// Everything known about how the two differ, in a few lines. `what` names the
// frame or the configuration; `width` turns an index into a row and a column.
inline void report(const char *what, const Difference &d, std::size_t height,
                   std::size_t width) {
  std::printf("  %s: cpu %zu signal pixels, gpu %zu\n", what, d.cpu_pixels,
              d.gpu_pixels);
  std::printf("    %zu found by both, %zu only by the cpu, %zu only by the "
              "gpu\n",
              d.common, d.only_cpu, d.only_gpu);
  if (d.integer_differs > 0) {
    std::printf("    %zu of the shared pixels differ in value or population, "
                "which is not rounding\n",
                d.integer_differs);
  }
  if (d.background_differs > 0) {
    std::printf("    %zu differ in background: worst %.3g absolute, %.3g "
                "relative\n",
                d.background_differs, d.worst_absolute, d.worst_relative);
  }
  if (d.have_example) {
    const SignalPixel &a = d.cpu_example;
    const SignalPixel &b = d.gpu_example;
    std::printf("    first: cpu index %u value %u background %.9g population "
                "%u\n",
                a.index, a.value, static_cast<double>(a.background),
                a.population);
    std::printf("           gpu index %u value %u background %.9g population "
                "%u\n",
                b.index, b.value, static_cast<double>(b.background),
                b.population);
    const unsigned index =
        d.only_gpu > 0 && d.only_cpu == 0 ? b.index : a.index;
    std::printf("           row %u, column %u of %zu x %zu\n",
                index / static_cast<unsigned>(width),
                index % static_cast<unsigned>(width), height, width);
  }
}

} // namespace compare_signal

#endif // SPOTFINDER_COMPARE_SIGNAL_HH
