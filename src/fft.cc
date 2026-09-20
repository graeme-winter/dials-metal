#include "fft.h"

#include <cmath>
#include <stdexcept>

namespace mxi {

// Defined here rather than in fft_fftw.cc so that it exists whichever
// transform is compiled, and setting it is harmless when nothing reads it.
std::size_t g_fft_threads = 0;

std::size_t next_power_of_two(std::size_t n) {
  std::size_t p = 1;
  while (p < n) p <<= 1;
  return p;
}

namespace {

// Iterative Cooley-Tukey on a strided line, so the same routine serves all
// three axes without transposing the grid between passes.
void transform_line(std::complex<double> *data, std::size_t n,
                    std::size_t stride, int sign) {
  // Bit reversal.
  for (std::size_t i = 1, j = 0; i < n; ++i) {
    std::size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(data[i * stride], data[j * stride]);
  }
  for (std::size_t len = 2; len <= n; len <<= 1) {
    const double angle = sign * 2.0 * M_PI / static_cast<double>(len);
    const std::complex<double> step(std::cos(angle), std::sin(angle));
    for (std::size_t i = 0; i < n; i += len) {
      std::complex<double> w(1.0, 0.0);
      for (std::size_t k = 0; k < len / 2; ++k) {
        std::complex<double> &u = data[(i + k) * stride];
        std::complex<double> &v = data[(i + k + len / 2) * stride];
        const std::complex<double> t = w * v;
        v = u - t;
        u = u + t;
        w *= step;
      }
    }
  }
}

}  // namespace

void fft3d_builtin(std::vector<std::complex<double>> &grid, std::size_t n,
                   int sign) {
  if ((n & (n - 1)) != 0 || n == 0) {
    throw std::invalid_argument("fft3d needs a power-of-two grid size");
  }
  if (grid.size() != n * n * n) {
    throw std::invalid_argument("fft3d grid is the wrong size");
  }
  const std::size_t nn = n * n;
  // Slowest axis: stride n*n.
  for (std::size_t j = 0; j < n; ++j) {
    for (std::size_t k = 0; k < n; ++k) {
      transform_line(&grid[j * n + k], n, nn, sign);
    }
  }
  // Middle axis: stride n.
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = 0; k < n; ++k) {
      transform_line(&grid[i * nn + k], n, n, sign);
    }
  }
  // Fastest axis: contiguous.
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      transform_line(&grid[i * nn + j * n], n, 1, sign);
    }
  }
}

}  // namespace mxi

#ifndef MXI_USE_FFTW
namespace mxi {
// No library asked for, so the built-in is the transform.
void fft3d(std::vector<std::complex<double>> &grid, std::size_t n, int sign) {
  fft3d_builtin(grid, n, sign);
}
}  // namespace mxi
#endif
