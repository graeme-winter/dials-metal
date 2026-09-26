#include "linalg.hh"

#include <cmath>

namespace mxi {

bool solve_spd(double *a, double *b, std::size_t n) {
  // Cholesky, in place: a becomes its own lower factor, b becomes the answer.
  // Used for refinement normal equations, where a failure to factorise means
  // the problem is rank deficient -- some parameter the data does not
  // constrain -- and the caller must stop rather than take a pseudo-inverse
  // and report a confident wrong answer.
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      double sum = a[i * n + j];
      for (std::size_t k = 0; k < j; ++k)
        sum -= a[i * n + k] * a[j * n + k];
      if (i == j) {
        if (!(sum > 0.0))
          return false;
        a[i * n + j] = std::sqrt(sum);
      } else {
        a[i * n + j] = sum / a[j * n + j];
      }
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    double sum = b[i];
    for (std::size_t k = 0; k < i; ++k)
      sum -= a[i * n + k] * b[k];
    b[i] = sum / a[i * n + i];
  }
  for (std::size_t ii = n; ii-- > 0;) {
    double sum = b[ii];
    for (std::size_t k = ii + 1; k < n; ++k)
      sum -= a[k * n + ii] * b[k];
    b[ii] = sum / a[ii * n + ii];
  }
  return true;
}

} // namespace mxi
