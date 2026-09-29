// The fused threshold's tile arithmetic, against dext(), without a GPU.
//
// dext_fused_emulated() runs the phases the CUDA kernel runs, in its order,
// tile by tile on the CPU. Each frame here is thresholded both ways and every
// signal pixel must agree in every field, the background's float bits
// included: frames whose sizes are not multiples of the tile, spots on the
// edges and corners, masked and saturated pixels, a masked row and column as a
// module gap makes, and both pixel types.

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "dext.hh"
#include "dext_fused.hh"

namespace {

int failures = 0;
int frames = 0;

template <typename T>
std::vector<T> make(std::size_t h, std::size_t w, unsigned seed, double mean) {
  std::mt19937 rng(seed);
  std::poisson_distribution<int> background(mean);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  std::vector<T> f(h * w);
  for (auto &p : f)
    p = static_cast<T>(background(rng));
  // Spots, some on the edges and corners.
  const int spots = static_cast<int>(h * w / 400) + 3;
  for (int s = 0; s < spots; ++s) {
    const double ci = s == 0   ? 0.0
                      : s == 1 ? static_cast<double>(h - 1)
                               : u(rng) * static_cast<double>(h);
    const double cj = s == 0   ? 0.0
                      : s == 1 ? static_cast<double>(w - 1)
                               : u(rng) * static_cast<double>(w);
    const double peak = 20.0 + 500.0 * u(rng);
    for (int di = -4; di <= 4; ++di)
      for (int dj = -4; dj <= 4; ++dj) {
        const long i = static_cast<long>(ci) + di,
                   j = static_cast<long>(cj) + dj;
        if (i < 0 || j < 0 || i >= static_cast<long>(h) ||
            j >= static_cast<long>(w))
          continue;
        const double v = peak * std::exp(-(di * di + dj * dj) / 3.0);
        T &p = f[static_cast<std::size_t>(i) * w + static_cast<std::size_t>(j)];
        p = static_cast<T>(p + static_cast<T>(v));
      }
  }
  // Masked and saturated pixels, a module gap's row and column, and a patch
  // of zeroes.
  const T masked = static_cast<T>(~T(0) - 1), top = static_cast<T>(~T(0));
  for (auto &p : f) {
    const double x = u(rng);
    if (x < 0.005)
      p = masked;
    else if (x < 0.007)
      p = top;
  }
  if (h > 20)
    for (std::size_t j = 0; j < w; ++j)
      f[(h / 2) * w + j] = masked;
  if (w > 20)
    for (std::size_t i = 0; i < h; ++i)
      f[i * w + w / 3] = masked;
  for (std::size_t i = 0; i < std::min<std::size_t>(h, 6); ++i)
    for (std::size_t j = 0; j < std::min<std::size_t>(w, 9); ++j)
      f[i * w + j] = 0;
  return f;
}

template <typename T>
void compare(std::size_t h, std::size_t w, unsigned seed, double mean) {
  const std::vector<T> f = make<T>(h, w, seed, mean);
  std::vector<SignalPixel> a, b;
  const int ra = dext<T>(f.data(), a, h, w);
  const int rb = dext_fused::dext_fused_emulated<T>(f.data(), b, h, w);
  ++frames;
  const std::string what = std::to_string(h) + " x " + std::to_string(w) +
                           ", " + std::to_string(8 * sizeof(T)) +
                           "-bit, seed " + std::to_string(seed);
  if (ra != rb || a.size() != b.size()) {
    std::printf("FAIL %s: dext %d with %zu pixels, fused %d with %zu\n",
                what.c_str(), ra, a.size(), rb, b.size());
    ++failures;
    return;
  }
  for (std::size_t k = 0; k < a.size(); ++k)
    if (std::memcmp(&a[k], &b[k], sizeof(SignalPixel)) != 0) {
      std::printf("FAIL %s: pixel %zu differs: index %u/%u value %u/%u "
                  "background %.9g/%.9g "
                  "population %u/%u\n",
                  what.c_str(), k, a[k].index, b[k].index, a[k].value,
                  b[k].value, static_cast<double>(a[k].background),
                  static_cast<double>(b[k].background), a[k].population,
                  b[k].population);
      ++failures;
      return;
    }
}

} // namespace

int main() {
  const std::size_t sizes[][2] = {{1, 1},    {5, 7},    {32, 32},   {33, 65},
                                  {100, 37}, {67, 300}, {257, 511}, {64, 96}};
  unsigned seed = 1;
  std::size_t signal = 0;
  for (const auto &s : sizes)
    for (double mean : {0.3, 2.0, 12.0}) {
      compare<std::uint16_t>(s[0], s[1], seed, mean);
      compare<std::uint32_t>(s[0], s[1], seed, mean);
      ++seed;
    }
  // And one frame at a detector's size in 16 bits, the case that matters.
  {
    const std::vector<std::uint16_t> f =
        make<std::uint16_t>(3262, 3108, 99, 1.5);
    std::vector<SignalPixel> a, b;
    dext<std::uint16_t>(f.data(), a, 3262, 3108);
    dext_fused::dext_fused_emulated<std::uint16_t>(f.data(), b, 3262, 3108);
    ++frames;
    signal = a.size();
    if (a.size() != b.size() ||
        (!a.empty() && std::memcmp(a.data(), b.data(),
                                   a.size() * sizeof(SignalPixel)) != 0)) {
      std::printf("FAIL 3262 x 3108: dext %zu pixels, fused %zu\n", a.size(),
                  b.size());
      ++failures;
    }
  }
  std::printf("%s: the fused threshold against dext, %d frames (%zu signal "
              "pixels on the "
              "largest), %d failures\n",
              failures == 0 ? "PASS" : "FAILED", frames, signal, failures);
  return failures == 0 ? 0 : 1;
}
