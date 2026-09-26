// Narrow spots at KNOWN sub-image positions. Each image integrates the rocking
// curve over its own oscillation, so image k holds the part of the Gaussian in
// phi between k and k+1 -- the physics that matters when a spot is thinner
// than an image.
#include "integrate.hh"
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>
namespace mxi {

static double Phi(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

int run_program() {
  std::mt19937 rng(3);
  const double background = 0.3;
  const double signal = 2000.0;
  const double sxy = 1.3;
  const int nz = 9;
  std::printf("error = centre minus truth, in images, by where in its image "
              "the reflection lies\n\n");
  for (double offset : {0.0, 0.1, 0.2}) {
    std::printf(
        "=== foreground window centred %.1f images from the truth ===\n",
        offset);
    for (double sz : {0.15, 0.3, 0.5, 1.0}) {
      std::printf("sigma %.2f images ", sz);
      double tot[3] = {0, 0, 0};
      int ntot = 0;
      for (int bin = 0; bin < 10; ++bin) {
        double sum[3] = {0, 0, 0};
        int n = 0;
        for (int trial = 0; trial < 300; ++trial) {
          std::uniform_real_distribution<double> u(0.0, 0.1);
          const double phase = 0.1 * bin + u(rng); // sub-image position
          const double zt = 4.0 + phase;           // truth, in images
          const double cx = 10.0 + 0.3, cy = 10.0 - 0.2;
          Shoebox box;
          box.panel = 0;
          box.bbox[0] = 0;
          box.bbox[1] = 21;
          box.bbox[2] = 0;
          box.bbox[3] = 21;
          box.bbox[4] = 0;
          box.bbox[5] = nz;
          const std::size_t N = box.size();
          box.data.assign(N, 0.0f);
          box.background.assign(N, 0.0f);
          box.mask.assign(N, shoebox_mask::kValid | shoebox_mask::kBackground);
          // expected counts and a Poisson draw
          std::vector<double> mu(N, 0.0);
          double norm_xy = 0.0;
          for (int y = 0; y < 21; ++y)
            for (int x = 0; x < 21; ++x) {
              const double dx = (x + 0.5 - cx) / sxy, dy = (y + 0.5 - cy) / sxy;
              norm_xy += std::exp(-0.5 * (dx * dx + dy * dy));
            }
          for (int z = 0; z < nz; ++z) {
            const double fz =
                Phi((z + 1 - zt) / sz) -
                Phi((z - zt) / sz); // the image integrates the curve
            for (int y = 0; y < 21; ++y)
              for (int x = 0; x < 21; ++x) {
                const double dx = (x + 0.5 - cx) / sxy,
                             dy = (y + 0.5 - cy) / sxy;
                const std::size_t at = box.at(x, y, z);
                mu[at] = background + signal * fz *
                                          std::exp(-0.5 * (dx * dx + dy * dy)) /
                                          norm_xy;
                // The mask's rule: a voxel is foreground when its image centre
                // is within n_sigma of the PREDICTED position along the scan.
                const double half = std::max(3.0 * sz, 0.5);
                if (std::fabs(x + 0.5 - cx) < 5 &&
                    std::fabs(y + 0.5 - cy) < 5 &&
                    std::fabs((z + 0.5) - (zt + offset)) <= half)
                  box.mask[at] =
                      shoebox_mask::kValid | shoebox_mask::kForeground;
              }
          }
          for (std::size_t i = 0; i < N; ++i) {
            std::poisson_distribution<int> p(mu[i]);
            box.data[i] = (float)p(rng);
          }

          // The spot finder: raw counts over pixels above threshold, frame at
          // k + 0.5, and a single-frame spot put at the middle of its frame.
          const double thresh = background + 3.0 * std::sqrt(background);
          double w = 0, wz = 0;
          int zlo = nz, zhi = -1;
          for (int z = 0; z < nz; ++z)
            for (int y = 0; y < 21; ++y)
              for (int x = 0; x < 21; ++x) {
                const double c = box.data[box.at(x, y, z)];
                if (c > thresh && std::fabs(x + 0.5 - cx) < 5 &&
                    std::fabs(y + 0.5 - cy) < 5) {
                  w += c;
                  wz += c * (z + 0.5);
                  zlo = std::min(zlo, z);
                  zhi = std::max(zhi, z);
                }
              }
          if (!(w > 0))
            continue;
          double spot_z = wz / w;
          if (zhi == zlo)
            spot_z = zlo + 0.5;

          IntegrateOptions o;
          o.background.tuning = 1e6;
          const IntegratedReflection r = integrate_shoebox(&box, o);
          if (!r.centroid_valid || !r.unbiased_valid)
            continue;
          const double err[3] = {spot_z - zt, r.centroid_z - zt,
                                 r.unbiased_z - zt};
          for (int k = 0; k < 3; ++k) {
            sum[k] += err[k];
            tot[k] += err[k];
          }
          ++n;
          ++ntot;
        }
      }
      std::printf(" mean error: spot finder %+.4f  integrator clipped %+.4f  "
                  "unclipped %+.4f\n",
                  tot[0] / ntot, tot[1] / ntot, tot[2] / ntot);
    }
  }
  // main may fall off its end and return 0; run_program may not, and did
  // not have to while it was main.
  return 0;
}

} // namespace mxi

int main() { return mxi::run_program(); }
