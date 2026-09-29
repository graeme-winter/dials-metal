// The extended dispersion threshold's three stages fused into one pass over a
// tile, for the GPU: each 32 x 32 tile of output reads its pixels once, with
// the halo the three windows need, and keeps the two masks between the stages
// in shared memory instead of writing each to device memory and reading it
// back. Every window is summed separably, rows then columns, by all of the
// block's threads rather than by one thread a row.
//
// The phases are functions of one item each, compiled for the host and the
// device alike: the CUDA kernel spreads them over its threads, and
// dext_fused_emulated() below runs them in order on the CPU, which is how the
// tile arithmetic -- which halo pixel is which -- is tested against dext()
// without a GPU. The sums are integers, which add to the same value in any
// order, wrapping included; the tests on them are written as the three
// kernels write them, so that the answer is theirs, bit for bit.
//
// Regions, in pixels relative to the tile's origin:
//   R0, 52 x 52 from -10: the pixels stage 0 reads;
//   R1, 46 x 46 from -7:  stage 0's mask, which stage 1 reads;
//   R2, 42 x 42 from -5:  stage 1's mask and the pixels stage 2 reads.
// A pixel outside the frame is loaded as masked, which is exact: a masked
// pixel adds nothing to any sum, as a pixel outside the frame adds nothing;
// stage 0 marks it 1, so that stage 1's erode cannot grow through it, as stage
// 1 skips a neighbour outside the frame; and stage 2 leaves it out of the
// background, as it leaves out one outside the frame.

#ifndef SPOTFINDER_DEXT_FUSED_HH
#define SPOTFINDER_DEXT_FUSED_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "dext.hh"
#include "signal_pixel.hh"

#ifdef __CUDACC__
#define DEXT_FUSED_HD __host__ __device__
#else
#define DEXT_FUSED_HD
#endif

namespace dext_fused {

constexpr int kTile = 32;
constexpr int kR2 = kTile + 2 * 5; // 42: stage 2's 11 x 11 window
constexpr int kR1 = kR2 + 2 * 2;   // 46: stage 1's 5 x 5 window
constexpr int kR0 = kR1 + 2 * 3;   // 52: stage 0's 7 x 7 window
constexpr int kOrigin0 = -10, kOrigin1 = -7, kOrigin2 = -5;

template <typename T> DEXT_FUSED_HD constexpr T masked_value() {
  return static_cast<T>(~T(0) - 1);
}

DEXT_FUSED_HD inline float root(float value) {
#ifdef __CUDA_ARCH__
  return sqrtf(value);
#else
  return std::sqrt(value);
#endif
}
DEXT_FUSED_HD inline double root(double value) {
#ifdef __CUDA_ARCH__
  return sqrt(value);
#else
  return std::sqrt(value);
#endif
}

//: What one block keeps in shared memory. Stage 0's row sums and stage 2's are
//: never alive together, so they share their storage.
template <typename T> struct Tile {
  using W = accumulator_t<T>;
  T pixel[kR0][kR0];
  std::uint8_t mask0[kR1][kR1];
  std::uint8_t erode_row[kR1][kR2];
  std::uint8_t mask1[kR2][kR2];
  union {
    struct {
      std::uint8_t count[kR0][kR1];
      W sum[kR0][kR1], sum2[kR0][kR1];
    } s0;
    struct {
      std::uint8_t count[kR2][kTile];
      W sum[kR2][kTile];
    } s2;
  } rows;
};

//: Load R0's pixel (r, c), for a tile at (ti, tj) of a frame height x width.
template <typename T>
DEXT_FUSED_HD inline void load(Tile<T> &t, const T *frame, int height,
                               int width, int ti, int tj, int r, int c) {
  const int i = ti + kOrigin0 + r, j = tj + kOrigin0 + c;
  t.pixel[r][c] = (i >= 0 && i < height && j >= 0 && j < width)
                      ? frame[i * width + j]
                      : masked_value<T>();
}

//: Stage 0 along a row: the count, sum and sum of squares of the valid pixels
//: in R0's row r, columns c to c + 6; r < 52, c < 46.
template <typename T>
DEXT_FUSED_HD inline void stage0_row(Tile<T> &t, int r, int c) {
  using W = accumulator_t<T>;
  const T masked = masked_value<T>();
  unsigned count = 0;
  W sum = 0, sum2 = 0;
  for (int d = 0; d < 7; ++d) {
    const T pixel = t.pixel[r][c + d];
    const W valid = pixel >= masked ? W(0) : W(1);
    const W p = valid * static_cast<W>(pixel);
    count += static_cast<unsigned>(valid);
    sum += p;
    sum2 += p * p;
  }
  t.rows.s0.count[r][c] = static_cast<std::uint8_t>(count);
  t.rows.s0.sum[r][c] = sum;
  t.rows.s0.sum2[r][c] = sum2;
}

//: Stage 0's test at R1's (r, c), down the column of row sums: R0's rows r to
//: r + 6 at column c are the 7 x 7 window around R0's (r + 3, c + 3).
template <typename T>
DEXT_FUSED_HD inline void stage0_test(Tile<T> &t, int r, int c) {
  using W = accumulator_t<T>;
  using R = real_t<T>;
  const R sigma_b = R(6.0);
  W m_sum = 0, i_sum = 0, i2_sum = 0;
  for (int d = 0; d < 7; ++d) {
    m_sum += static_cast<W>(t.rows.s0.count[r + d][c]);
    i_sum += t.rows.s0.sum[r + d][c];
    i2_sum += t.rows.s0.sum2[r + d][c];
  }
  bool signal = false;
  if (m_sum >= 2) {
    // n * sum(i^2) - sum(i)^2 - (n - 1) * sum(i)
    //     > sigma_b * sum(i) * sqrt(2 * (n - 1)), as stage0 writes it
    const R m = static_cast<R>(m_sum);
    const R s = static_cast<R>(i_sum);
    const R s2 = static_cast<R>(i2_sum);
    signal = (m * s2 - s * s - s * (m - R(1))) >
             (s * sigma_b * root(R(2) * (m - R(1))));
  }
  t.mask0[r][c] =
      (signal || t.pixel[r + 3][c + 3] >= masked_value<T>()) ? 1 : 0;
}

//: Stage 1 along a row: whether R1's row r, columns c to c + 4, are all set;
//: r < 46, c < 42.
template <typename T>
DEXT_FUSED_HD inline void stage1_row(Tile<T> &t, int r, int c) {
  std::uint8_t all = 1;
  for (int d = 0; d < 5; ++d)
    all &= t.mask0[r][c + d];
  t.erode_row[r][c] = all;
}

//: Stage 1's erode at R2's (r, c): a masked pixel is 1, and otherwise the 5 x 5
//: neighbourhood of R1's (r + 2, c + 2) must all be set -- which includes the
//: centre, so a centre of 0 gives 0, as stage1 has it.
template <typename T>
DEXT_FUSED_HD inline void stage1_erode(Tile<T> &t, int r, int c) {
  if (t.pixel[r + 5][c + 5] >= masked_value<T>()) {
    t.mask1[r][c] = 1;
    return;
  }
  std::uint8_t all = 1;
  for (int d = 0; d < 5; ++d)
    all &= t.erode_row[r + d][c];
  t.mask1[r][c] = all;
}

//: Stage 2 along a row: the count and sum of the background pixels -- valid
//: and not flagged by stage 1 -- in R2's row r, columns c to c + 10; r < 42,
//: c < 32.
template <typename T>
DEXT_FUSED_HD inline void stage2_row(Tile<T> &t, int r, int c) {
  using W = accumulator_t<T>;
  const T masked = masked_value<T>();
  unsigned count = 0;
  W sum = 0;
  for (int d = 0; d < 11; ++d) {
    const T pixel = t.pixel[r + 5][c + d + 5];
    const W valid = (pixel >= masked) || t.mask1[r][c + d] ? W(0) : W(1);
    count += static_cast<unsigned>(valid);
    sum += valid * static_cast<W>(pixel);
  }
  t.rows.s2.count[r][c] = static_cast<std::uint8_t>(count);
  t.rows.s2.sum[r][c] = sum;
}

//: What stage 2 decides at the tile's output (r, c).
struct Decision {
  bool signal = false;
  std::uint32_t value = 0;
  float background = 0.0f;
  std::uint16_t population = 0;
};

//: Stage 2's Poisson test at the tile's (r, c), R2's (r + 5, c + 5), down the
//: column of row sums: R2's rows r to r + 10 at column c. `in_frame` is
//: whether that pixel is in the frame at all.
template <typename T>
DEXT_FUSED_HD inline Decision stage2_test(const Tile<T> &t, int r, int c,
                                          bool in_frame) {
  using W = accumulator_t<T>;
  using R = real_t<T>;
  const R sigma_s = R(3.0);
  W m_sum = 0, i_sum = 0;
  for (int d = 0; d < 11; ++d) {
    m_sum += static_cast<W>(t.rows.s2.count[r + d][c]);
    i_sum += t.rows.s2.sum[r + d][c];
  }
  Decision out;
  if (!in_frame)
    return out;
  const T raw = t.pixel[r + 10][c + 10];
  const T p = raw >= masked_value<T>() ? T(0) : raw;
  if (p > 0 && t.mask1[r + 5][c + 5]) {
    const R mean =
        m_sum >= 2 ? static_cast<R>(i_sum) / static_cast<R>(m_sum) : R(0);
    out.signal = static_cast<R>(p) >= (mean + sigma_s * root(mean));
    out.value = static_cast<std::uint32_t>(p);
    out.background = static_cast<float>(mean);
    out.population = static_cast<std::uint16_t>(m_sum);
  }
  return out;
}

//: The fused threshold on the CPU, tile by tile, each phase run in the order
//: the kernel runs it: for testing the tile arithmetic against dext(), not for
//: speed. Returns what dext() returns, the survivors ascending by index.
template <typename T>
int dext_fused_emulated(const T *frame, std::vector<SignalPixel> &out,
                        std::size_t height, std::size_t width) {
  out.clear();
  if (frame == nullptr || height == 0 || width == 0 ||
      height > 0x7fffffffu / width)
    return -1;
  const int h = static_cast<int>(height), w = static_cast<int>(width);
  const auto tile = std::make_unique<Tile<T>>();
  Tile<T> &t = *tile;
  for (int ti = 0; ti < h; ti += kTile)
    for (int tj = 0; tj < w; tj += kTile) {
      for (int r = 0; r < kR0; ++r)
        for (int c = 0; c < kR0; ++c)
          load(t, frame, h, w, ti, tj, r, c);
      for (int r = 0; r < kR0; ++r)
        for (int c = 0; c < kR1; ++c)
          stage0_row(t, r, c);
      for (int r = 0; r < kR1; ++r)
        for (int c = 0; c < kR1; ++c)
          stage0_test(t, r, c);
      for (int r = 0; r < kR1; ++r)
        for (int c = 0; c < kR2; ++c)
          stage1_row(t, r, c);
      for (int r = 0; r < kR2; ++r)
        for (int c = 0; c < kR2; ++c)
          stage1_erode(t, r, c);
      for (int r = 0; r < kR2; ++r)
        for (int c = 0; c < kTile; ++c)
          stage2_row(t, r, c);
      for (int r = 0; r < kTile; ++r)
        for (int c = 0; c < kTile; ++c) {
          const int i = ti + r, j = tj + c;
          const Decision d = stage2_test(t, r, c, i < h && j < w);
          if (d.signal) {
            SignalPixel pixel;
            pixel.index = static_cast<std::uint32_t>(i * w + j);
            pixel.value = d.value;
            pixel.background = d.background;
            pixel.population = d.population;
            pixel.reserved = 0;
            out.push_back(pixel);
          }
        }
    }
  std::sort(out.begin(), out.end(),
            [](const SignalPixel &a, const SignalPixel &b) {
              return a.index < b.index;
            });
  return 0;
}

} // namespace dext_fused

#endif // SPOTFINDER_DEXT_FUSED_HH
