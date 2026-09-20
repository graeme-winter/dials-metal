// The three-dimensional transform through FFTW, when it is available.
//
// The built-in radix-2 in fft.cc is correct, has no dependencies and is tested
// against a direct transform. It is also the largest single phase of indexing
// on a fast machine -- 57 per cent of one measured run -- and a decomposition
// written in a few dozen lines is not going to compete with a library that has
// been tuned for twenty-five years.
//
// So this is a drop-in alternative, chosen at build time. Same signature, same
// normalisation, same answer to within rounding; `MXI_FFTW` in CMake decides
// which is compiled, and there is a test that the two agree.
//
// FFTW is GPL, which the licence of anything linking it must accommodate. That
// is why it is optional and off unless asked for, rather than a dependency.

#include <fftw3.h>

#include <complex>
#include <cstddef>
#include <mutex>
#include <vector>

#include "fft.h"

namespace mxi {

namespace {

// FFTW's planner is not thread safe, though the execution is. One lock around
// planning only, as the documentation asks.
std::mutex g_plan_lock;

}  // namespace

void fft3d(std::vector<std::complex<double>> &grid, std::size_t n, int sign) {
  if (n == 0) return;
  if (grid.size() != n * n * n) return;

  // fftw_complex and std::complex<double> are required to have the same
  // layout, which is what makes this a cast rather than a copy.
  static_assert(sizeof(std::complex<double>) == sizeof(fftw_complex),
                "std::complex<double> must be layout compatible with "
                "fftw_complex for this cast to be legitimate");
  auto *data = reinterpret_cast<fftw_complex *>(grid.data());

  fftw_plan plan;
  {
    const std::lock_guard<std::mutex> held(g_plan_lock);
    // FFTW_ESTIMATE rather than FFTW_MEASURE: measuring writes to the array
    // and takes longer than the transform it is planning, for a grid used a
    // handful of times per run.
    plan = fftw_plan_dft_3d(static_cast<int>(n), static_cast<int>(n),
                            static_cast<int>(n), data, data,
                            // FFTW_FORWARD is exp(-i...), and the built-in
                            // uses exp(+i * sign * ...), so a positive sign
                            // here is FFTW's BACKWARD. Getting this the
                            // natural-looking way round conjugates every
                            // peak's position through the origin, which on a
                            // centrosymmetric lattice still looks like a
                            // lattice -- the agreement test caught it, nothing
                            // else would have.
                            sign >= 0 ? FFTW_BACKWARD : FFTW_FORWARD,
                            FFTW_ESTIMATE);
  }
  if (plan == nullptr) return;
  fftw_execute(plan);
  {
    const std::lock_guard<std::mutex> held(g_plan_lock);
    fftw_destroy_plan(plan);
  }
}

}  // namespace mxi
