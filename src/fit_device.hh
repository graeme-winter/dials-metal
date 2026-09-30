// Profile fitting in single precision, for a device: the reference profile
// carried onto a shoebox's pixels, and the weighted least squares against
// them -- what reference.cc's profile_on_pixels and fit_on_pixels do in double,
// written again in float for a GPU.
//
// Each step is a function of one item, compiled for the host and for a device
// alike: a kernel spreads them over its threads, and fit_box_emulated() below
// runs them in order on the CPU, which is how the arithmetic is tested against
// the double fit without a device. Two departures from the CPU, both because
// of single precision:
//
// * no exact recomputation near a cell boundary. The CPU interpolates a
//   subdivision's position from the pixel's corners and computes it exactly
//   when it falls within a thousandth of a cell of a boundary, so that its
//   cells are the direct computation's; in float the exact computation is not
//   exact either, so it is left out, and a subdivision near a boundary may fall
//   in the other cell -- measured, not assumed, by the emulation;
// * sums in float. Over a few thousand voxels a box that costs a few parts in
//   ten million, against counting statistics of parts in a hundred.
//
// Laid out as the CPU lays it out: a box of nx x ny x nz voxels, x fastest;
// the reference a side^3 grid, i1 fastest, a slice of the detector's face
// side^2 and the planes along the rotation; a pixel's subdivisions sub x sub.

#ifndef MXI_FIT_DEVICE_HH
#define MXI_FIT_DEVICE_HH

#include <cmath>
#include <cstdint>

#ifdef __CUDACC__
#define FIT_HD __host__ __device__
#else
#define FIT_HD
#endif

namespace fitdev {

// shoebox_mask's bits and reference.cc's least variance; the host checks they
// agree (fit_batch.cc).
constexpr std::uint8_t kValid = 1u << 0, kBackground = 1u << 1,
                       kForeground = 1u << 2;
constexpr float kLeastVariance = 0.05f;
constexpr int kMaxSide = 17; //: grids up to n = 8
constexpr int kMaxSub = 7;   //: subdivisions up to 7 x 7
constexpr float kDegrees = 57.29577951308232f;

struct F3 {
  float x = 0.0f, y = 0.0f, z = 0.0f;
};
FIT_HD inline F3 operator+(F3 a, F3 b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
FIT_HD inline F3 operator-(F3 a, F3 b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
FIT_HD inline F3 operator*(F3 a, float s) {
  return {a.x * s, a.y * s, a.z * s};
}
FIT_HD inline float dot(F3 a, F3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
FIT_HD inline F3 cross(F3 a, F3 b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
FIT_HD inline float f_sqrt(float v) {
#ifdef __CUDA_ARCH__
  return sqrtf(v);
#else
  return std::sqrt(v);
#endif
}
FIT_HD inline float f_exp(float v) {
#ifdef __CUDA_ARCH__
  return expf(v);
#else
  return std::exp(v);
#endif
}
FIT_HD inline float f_floor(float v) {
#ifdef __CUDA_ARCH__
  return floorf(v);
#else
  return std::floor(v);
#endif
}
FIT_HD inline float f_abs(float v) { return v < 0.0f ? -v : v; }
FIT_HD inline float norm(F3 a) { return f_sqrt(dot(a, a)); }

//: A detector panel, as Panel keeps it: millimetres, and the parallax
//: correction's sensor.
struct PanelF {
  F3 origin, fast, slow, normal;
  float pixel_fast = 0.0f, pixel_slow = 0.0f;
  float mu = 0.0f, thickness = 0.0f;
  int parallax = 0, conditional = 0;
};

//: What every box of a batch shares: the beam, the rotation, the scan, the
//: grid and the fit's settings.
struct Setup {
  F3 s0;       //: the beam, 1/A, towards the detector
  F3 axis_lab; //: the rotation axis in the laboratory, unit length
  float osc_start = 0.0f, osc = 0.0f,
        z_offset = 0.0f; //: radians, radians an image
  int side = 9, sub = 5;
  float span_d = 0.0f, step_d = 0.0f, span_m = 0.0f, step_m = 0.0f; //: degrees
  float gain = 1.0f;
  int iterations = 3;
};

//: One box: where it is, its reflection, its constant background, and where its
//: voxels and its local reference are in the batch's arrays.
struct BoxF {
  int x0 = 0, y0 = 0, z0 = 0, nx = 0, ny = 0, nz = 0, panel = 0;
  F3 s1;
  float phi = 0.0f, background = 0.0f;
  std::uint32_t voxels_at = 0;    //: into data, mask and the profile scratch
  std::uint32_t reference_at = 0; //: into the local references, side^3 each
  std::uint32_t corners_at = 0;   //: into the corner scratch, (nx + 1)(ny + 1)
};

//: What a fit gives back, as ProfileFit.
struct FitF {
  int valid = 0, iterations = 0;
  float intensity = 0.0f, variance = 0.0f, correlation = 0.0f, measured = 1.0f;
};

//: The laboratory position of a point on a panel given in pixels, parallax
//: included: Panel::px_to_mm and lab_coord_mm.
FIT_HD inline F3 lab_px(const PanelF &p, float px, float py) {
  const float raw_f = px * p.pixel_fast, raw_s = py * p.pixel_slow;
  float off_f = 0.0f, off_s = 0.0f;
  if (p.parallax && p.mu > 0.0f && p.thickness > 0.0f) {
    const F3 lab = p.origin + p.fast * raw_f + p.slow * raw_s;
    const F3 u = lab * (1.0f / norm(lab));
    const float cos_theta = f_abs(dot(u, p.normal));
    if (cos_theta > 0.0f) {
      const float attenuation = 1.0f / p.mu;
      const float path = p.thickness / cos_theta;
      const float transmitted = f_exp(-p.mu * path);
      float depth = attenuation - (path + attenuation) * transmitted;
      if (p.conditional) {
        const float absorbed = 1.0f - transmitted;
        if (absorbed > 0.0f)
          depth /= absorbed;
      }
      off_f = depth * dot(u, p.fast);
      off_s = depth * dot(u, p.slow);
    }
  }
  return p.origin + p.fast * (raw_f - off_f) + p.slow * (raw_s - off_s);
}

//: The Kabsch frame of a reflection: kabsch_frame.
struct FrameF {
  int valid = 0;
  F3 e1, e2, s1;
  float zeta = 0.0f, s1_length = 0.0f;
};
FIT_HD inline FrameF frame_of(const Setup &setup, F3 s1) {
  FrameF out;
  const F3 c = cross(s1, setup.s0);
  const float lc = norm(c), ls = norm(s1);
  if (!(lc > 0.0f) || !(ls > 0.0f))
    return out;
  out.e1 = c * (1.0f / lc);
  const F3 second = cross(s1, out.e1);
  const float l2 = norm(second);
  if (!(l2 > 0.0f))
    return out;
  out.e2 = second * (1.0f / l2);
  out.s1 = s1;
  out.s1_length = ls;
  out.zeta = dot(setup.axis_lab, out.e1);
  out.valid = 1;
  return out;
}

//: A point's e1 and e2, degrees: epsilon_of.
FIT_HD inline void epsilon12(const FrameF &frame, const PanelF &p, float px,
                             float py, float *e1, float *e2) {
  const F3 lab = lab_px(p, px, py);
  const float length = norm(lab);
  if (!(length > 0.0f)) {
    *e1 = *e2 = 1e30f;
    return;
  }
  const F3 d = lab * (frame.s1_length / length) - frame.s1;
  *e1 = kDegrees * dot(frame.e1, d) / frame.s1_length;
  *e2 = kDegrees * dot(frame.e2, d) / frame.s1_length;
}

//: Step 1: one corner of the box's pixels, (cx, cy) of (nx + 1)(ny + 1).
FIT_HD inline void corner(const FrameF &frame, const PanelF &p, const BoxF &box,
                          int cx, int cy, float *c1, float *c2) {
  const int at = cy * (box.nx + 1) + cx;
  epsilon12(frame, p, static_cast<float>(box.x0 + cx),
            static_cast<float>(box.y0 + cy), &c1[at], &c2[at]);
}

//: Step 2: the cells one pixel's subdivisions fall in, and how many fall in
//: each, in order of first appearance: pixel_cells, without its exact
//: recomputation near a boundary. Returns how many cells.
FIT_HD inline int pixel_cells(const Setup &setup, const BoxF &box,
                              const float *c1, const float *c2, int x, int y,
                              std::uint16_t *ids, std::uint16_t *counts) {
  const int cw = box.nx + 1;
  const int a00 = y * cw + x, a10 = a00 + 1, a01 = a00 + cw, a11 = a01 + 1;
  const int sub = setup.sub;
  int found = 0;
  for (int sy = 0; sy < sub; ++sy) {
    const float v = (static_cast<float>(sy) + 0.5f) / static_cast<float>(sub);
    for (int sx = 0; sx < sub; ++sx) {
      const float u = (static_cast<float>(sx) + 0.5f) / static_cast<float>(sub);
      const float w00 = (1.0f - u) * (1.0f - v), w10 = u * (1.0f - v),
                  w01 = (1.0f - u) * v, w11 = u * v;
      const float e1 =
          w00 * c1[a00] + w10 * c1[a10] + w01 * c1[a01] + w11 * c1[a11];
      const float e2 =
          w00 * c2[a00] + w10 * c2[a10] + w01 * c2[a01] + w11 * c2[a11];
      const int i1 =
          static_cast<int>(f_floor((e1 + setup.span_d) / setup.step_d));
      const int i2 =
          static_cast<int>(f_floor((e2 + setup.span_d) / setup.step_d));
      if (i1 < 0 || i1 >= setup.side || i2 < 0 || i2 >= setup.side)
        continue;
      const std::uint16_t id = static_cast<std::uint16_t>(i2 * setup.side + i1);
      int q = 0;
      while (q < found && ids[q] != id)
        ++q;
      if (q == found) {
        ids[found] = id;
        counts[found] = 0;
        ++found;
      }
      ++counts[q];
    }
  }
  return found;
}

//: Step 3: which of the grid's planes one z-slice of the box covers, and by
//: how much: plane_weights. Returns how many planes; -1 if the slice spans no
//: rotation at all.
FIT_HD inline int plane_weights(const Setup &setup, const FrameF &frame,
                                const BoxF &box, int z, int *js, float *ws) {
  const float image = static_cast<float>(box.z0 + z);
  const float phi_low = setup.osc_start + (image - setup.z_offset) * setup.osc;
  const float phi_high = phi_low + setup.osc;
  float lo = kDegrees * (frame.zeta * (phi_low - box.phi));
  float hi = kDegrees * (frame.zeta * (phi_high - box.phi));
  if (lo > hi) {
    const float t = lo;
    lo = hi;
    hi = t;
  }
  const float span = hi - lo;
  if (!(span > 0.0f))
    return -1;
  const int j_low =
      static_cast<int>(f_floor((lo + setup.span_m) / setup.step_m));
  const int j_high =
      static_cast<int>(f_floor((hi + setup.span_m) / setup.step_m));
  int n = 0;
  for (int j = j_low < 0 ? 0 : j_low;
       j <= (j_high < setup.side - 1 ? j_high : setup.side - 1); ++j) {
    const float plane_low =
        -setup.span_m + static_cast<float>(j) * setup.step_m;
    const float plane_high = plane_low + setup.step_m;
    const float overlap =
        (hi < plane_high ? hi : plane_high) - (lo > plane_low ? lo : plane_low);
    if (overlap > 0.0f) {
      js[n] = j;
      ws[n] = overlap / span;
      ++n;
    }
  }
  return n;
}

//: Step 4: element q of a slice's reference, its planes blended.
FIT_HD inline float blended(const Setup &setup, const float *reference,
                            const int *js, const float *ws, int planes, int q) {
  const int slice = setup.side * setup.side;
  float value = 0.0f;
  for (int k = 0; k < planes; ++k)
    value += ws[k] * reference[js[k] * slice + q];
  return value;
}

//: Step 5: one voxel of the profile on the pixels, from its pixel's cells and
//: its slice's blended reference.
FIT_HD inline float pixel_value(const std::uint16_t *ids,
                                const std::uint16_t *counts, int found,
                                const float *slice, float share) {
  float value = 0.0f;
  for (int q = 0; q < found; ++q)
    value += static_cast<float>(counts[q]) * slice[ids[q]];
  return share * value;
}

//: Step 6: the weighted least squares against the profile on the pixels:
//: fit_on_pixels, its sums serial here; a kernel's are the same terms, summed
//: across its threads.
FIT_HD inline FitF fit_profile(const Setup &setup, const BoxF &box,
                               const float *data, const std::uint8_t *mask,
                               const float *profile) {
  FitF out;
  const int voxels = box.nx * box.ny * box.nz;
  const float b = box.background;
  float profile_sum = 0.0f, seen = 0.0f;
  for (int i = 0; i < voxels; ++i) {
    if ((mask[i] & kForeground) == 0)
      continue;
    profile_sum += profile[i];
    if (mask[i] & kValid)
      seen += profile[i];
  }
  if (!(profile_sum > 0.0f) || !(seen > 0.0f))
    return out;
  out.measured = seen / profile_sum;
  float scale = 0.0f;
  for (int i = 0; i < voxels; ++i)
    if ((mask[i] & kForeground) && (mask[i] & kValid))
      scale += data[i] - b;
  scale /= profile_sum;
  float variance = 0.0f;
  const int rounds = setup.iterations > 1 ? setup.iterations : 1;
  for (int round = 0; round < rounds; ++round) {
    float numerator = 0.0f, denominator = 0.0f;
    for (int i = 0; i < voxels; ++i) {
      if (!(mask[i] & kForeground) || !(mask[i] & kValid))
        continue;
      const float expected = b + (scale > 0.0f ? scale : 0.0f) * profile[i];
      const float v =
          setup.gain * (expected > kLeastVariance ? expected : kLeastVariance);
      numerator += profile[i] * (data[i] - b) / v;
      denominator += profile[i] * profile[i] / v;
    }
    if (!(denominator > 0.0f))
      return out;
    scale = numerator / denominator;
    variance = 1.0f / denominator;
    out.iterations = round + 1;
  }
  out.intensity = scale * profile_sum;
  out.variance = variance * profile_sum * profile_sum;
  float foreground = 0.0f, background_pixels = 0.0f, background_sum = 0.0f;
  for (int i = 0; i < voxels; ++i) {
    if ((mask[i] & kValid) == 0)
      continue;
    if (mask[i] & kForeground) {
      foreground += 1.0f;
      background_sum += b;
    } else if (mask[i] & kBackground) {
      background_pixels += 1.0f;
    }
  }
  if (background_pixels > 0.0f)
    out.variance +=
        setup.gain * (foreground / background_pixels) * background_sum;
  float n = 0.0f, sp = 0.0f, sd = 0.0f, spp = 0.0f, sdd = 0.0f, spd = 0.0f;
  for (int i = 0; i < voxels; ++i) {
    if (!(mask[i] & kForeground) || !(mask[i] & kValid))
      continue;
    const float a = profile[i], d = data[i] - b;
    n += 1.0f;
    sp += a;
    sd += d;
    spp += a * a;
    sdd += d * d;
    spd += a * d;
  }
  if (n > 1.0f) {
    const float cov = spd / n - (sp / n) * (sd / n);
    const float va = spp / n - (sp / n) * (sp / n),
                vb = sdd / n - (sd / n) * (sd / n);
    if (va > 0.0f && vb > 0.0f)
      out.correlation = cov / f_sqrt(va * vb);
  }
  out.valid = 1;
  return out;
}

//: One box on the CPU, the steps in the order a kernel runs them: for testing
//: the arithmetic against the double fit, and for running --gpu's path where
//: there is no device. `corners` holds (nx + 1)(ny + 1) floats twice, `profile`
//: the box's voxels, `blend_scratch` side^2 floats a slice.
inline FitF fit_box_emulated(const Setup &setup, const PanelF *panels,
                             const BoxF &box, const float *data,
                             const std::uint8_t *mask, const float *reference,
                             float *corners, float *profile,
                             float *blend_scratch) {
  FitF none;
  const FrameF frame = frame_of(setup, box.s1);
  if (!frame.valid)
    return none;
  const PanelF &p = panels[box.panel];
  const int cn = (box.nx + 1) * (box.ny + 1);
  float *c1 = corners, *c2 = corners + cn;
  for (int cy = 0; cy <= box.ny; ++cy)
    for (int cx = 0; cx <= box.nx; ++cx)
      corner(frame, p, box, cx, cy, c1, c2);
  // As the kernels do it: every slice's blended reference, each element with
  // its own slice's plane weights; then each pixel's cells found once and
  // walked down all its slices. The same sums in the same order as a slice at
  // a time.
  const float share = 1.0f / static_cast<float>(setup.sub * setup.sub);
  const int slice_size = setup.side * setup.side;
  std::uint16_t ids[kMaxSub * kMaxSub], counts[kMaxSub * kMaxSub];
  int js[kMaxSide];
  float ws[kMaxSide];
  float *blend = blend_scratch;
  for (int e = 0; e < box.nz * slice_size; ++e) {
    const int z = e / slice_size, q = e % slice_size;
    const int planes = plane_weights(setup, frame, box, z, js, ws);
    blend[e] = planes > 0 ? blended(setup, reference, js, ws, planes, q) : 0.0f;
  }
  const int pixels = box.nx * box.ny;
  for (int k = 0; k < pixels; ++k) {
    const int found =
        pixel_cells(setup, box, c1, c2, k % box.nx, k / box.nx, ids, counts);
    for (int z = 0; z < box.nz; ++z)
      profile[z * pixels + k] =
          pixel_value(ids, counts, found, blend + z * slice_size, share);
  }
  return fit_profile(setup, box, data, mask, profile);
}

} // namespace fitdev

#endif // MXI_FIT_DEVICE_HH
