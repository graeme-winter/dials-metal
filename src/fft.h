// A three-dimensional complex FFT, radix two.
//
// Written rather than linked because it is sixty lines and the alternative is
// a dependency on a laptop and on a beamline machine. A 256^3 transform is
// 16.7 million points and takes a couple of seconds here, which is a small
// part of indexing; if it ever stops being small, this is a self-contained
// thing to replace.
//
// Note for a device port: the 1D passes along each axis are independent, so
// the natural decomposition is one thread per line. That is also the reason
// this is separated out rather than written inline in the indexer.

#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace mxi {

// Smallest power of two not less than n.
std::size_t next_power_of_two(std::size_t n);

// In-place transform of an n x n x n grid stored with the last index fastest.
// `sign` is -1 for the forward transform and +1 for the inverse; no scaling is
// applied in either direction, because every use here looks at the modulus of
// the result and a constant factor never matters.
void fft3d(std::vector<std::complex<double>> &grid, std::size_t n, int sign);

}  // namespace mxi
