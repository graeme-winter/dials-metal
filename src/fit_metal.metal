// Profile fitting on a Metal device, in single precision: the steps of
// fit_device.hh, transcribed -- the shading language has no <cmath> and wants
// every pointer's address space, so it cannot include the header -- and
// spread over a threadgroup a box as fit_cuda.cu spreads them. fit_device.hh
// says what each step is and how it departs from the CPU's double fit;
// --gpu-emulate runs the same steps in order on the CPU. The structs match
// fit_device.hh's field for field, all four-byte fields in the same order, and
// fit_batch.cc checks their sizes.

#include <metal_stdlib>
using namespace metal;

constant uchar kValid = 1, kBackground = 2, kForeground = 4;
constant float kLeastVariance = 0.05f;
constant float kDegrees = 57.29577951308232f;
#define K_MAX_SIDE 17
#define K_MAX_SUB 7
#define K_THREADS 128
// Every slice's blended reference at once, in threadgroup memory, when it fits:
// 6144 floats, 24 KB, of the 32 a threadgroup has.
#define K_BLEND 6144

struct F3 { float x, y, z; };
static inline F3 add(F3 a, F3 b) { return F3{a.x + b.x, a.y + b.y, a.z + b.z}; }
static inline F3 sub(F3 a, F3 b) { return F3{a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline F3 scale3(F3 a, float s) { return F3{a.x * s, a.y * s, a.z * s}; }
static inline float dot3(F3 a, F3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline F3 cross3(F3 a, F3 b) {
  return F3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
static inline float norm3(F3 a) { return sqrt(dot3(a, a)); }

struct PanelF {
  F3 origin, fast, slow, normal;
  float pixel_fast, pixel_slow, mu, thickness;
  int parallax, conditional;
};
struct Setup {
  F3 s0, axis_lab;
  float osc_start, osc, z_offset;
  int side, sub;
  float span_d, step_d, span_m, step_m;
  float gain;
  int iterations;
};
struct BoxF {
  int x0, y0, z0, nx, ny, nz, panel;
  F3 s1;
  float phi, background;
  uint voxels_at, reference_at, corners_at;
};
struct FitF {
  int valid, iterations;
  float intensity, variance, correlation, measured;
};
struct FrameF {
  int valid;
  F3 e1, e2, s1;
  float zeta, s1_length;
};

static F3 lab_px(PanelF p, float px, float py) {
  const float raw_f = px * p.pixel_fast, raw_s = py * p.pixel_slow;
  float off_f = 0.0f, off_s = 0.0f;
  if (p.parallax != 0 && p.mu > 0.0f && p.thickness > 0.0f) {
    const F3 lab = add(add(p.origin, scale3(p.fast, raw_f)), scale3(p.slow, raw_s));
    const F3 u = scale3(lab, 1.0f / norm3(lab));
    const float cos_theta = fabs(dot3(u, p.normal));
    if (cos_theta > 0.0f) {
      const float attenuation = 1.0f / p.mu;
      const float path = p.thickness / cos_theta;
      const float transmitted = exp(-p.mu * path);
      float depth = attenuation - (path + attenuation) * transmitted;
      if (p.conditional != 0) {
        const float absorbed = 1.0f - transmitted;
        if (absorbed > 0.0f)
          depth /= absorbed;
      }
      off_f = depth * dot3(u, p.fast);
      off_s = depth * dot3(u, p.slow);
    }
  }
  return add(add(p.origin, scale3(p.fast, raw_f - off_f)), scale3(p.slow, raw_s - off_s));
}

static FrameF frame_of(constant Setup &setup, F3 s1) {
  FrameF out;
  out.valid = 0;
  out.zeta = 0.0f;
  out.s1_length = 0.0f;
  out.e1 = out.e2 = out.s1 = F3{0.0f, 0.0f, 0.0f};
  const F3 c = cross3(s1, setup.s0);
  const float lc = norm3(c), ls = norm3(s1);
  if (!(lc > 0.0f) || !(ls > 0.0f))
    return out;
  out.e1 = scale3(c, 1.0f / lc);
  const F3 second = cross3(s1, out.e1);
  const float l2 = norm3(second);
  if (!(l2 > 0.0f))
    return out;
  out.e2 = scale3(second, 1.0f / l2);
  out.s1 = s1;
  out.s1_length = ls;
  out.zeta = dot3(setup.axis_lab, out.e1);
  out.valid = 1;
  return out;
}

static void epsilon12(FrameF frame, PanelF p, float px, float py, thread float &e1,
                      thread float &e2) {
  const F3 lab = lab_px(p, px, py);
  const float length = norm3(lab);
  if (!(length > 0.0f)) {
    e1 = e2 = 1e30f;
    return;
  }
  const F3 d = sub(scale3(lab, frame.s1_length / length), frame.s1);
  e1 = kDegrees * dot3(frame.e1, d) / frame.s1_length;
  e2 = kDegrees * dot3(frame.e2, d) / frame.s1_length;
}

static int pixel_cells(constant Setup &setup, BoxF box, device const float *c1,
                       device const float *c2, int x, int y, thread ushort *ids,
                       thread ushort *counts) {
  const int cw = box.nx + 1;
  const int a00 = y * cw + x, a10 = a00 + 1, a01 = a00 + cw, a11 = a01 + 1;
  const int sub = setup.sub;
  int found = 0;
  for (int sy = 0; sy < sub; ++sy) {
    const float v = (float(sy) + 0.5f) / float(sub);
    for (int sx = 0; sx < sub; ++sx) {
      const float u = (float(sx) + 0.5f) / float(sub);
      const float w00 = (1.0f - u) * (1.0f - v), w10 = u * (1.0f - v), w01 = (1.0f - u) * v,
                  w11 = u * v;
      const float e1 = w00 * c1[a00] + w10 * c1[a10] + w01 * c1[a01] + w11 * c1[a11];
      const float e2 = w00 * c2[a00] + w10 * c2[a10] + w01 * c2[a01] + w11 * c2[a11];
      const int i1 = int(floor((e1 + setup.span_d) / setup.step_d));
      const int i2 = int(floor((e2 + setup.span_d) / setup.step_d));
      if (i1 < 0 || i1 >= setup.side || i2 < 0 || i2 >= setup.side)
        continue;
      const ushort id = ushort(i2 * setup.side + i1);
      int q = 0;
      while (q < found && ids[q] != id)
        ++q;
      if (q == found) {
        ids[found] = id;
        counts[found] = 0;
        ++found;
      }
      counts[q] += 1;
    }
  }
  return found;
}

static int plane_weights(constant Setup &setup, FrameF frame, BoxF box, int z,
                         threadgroup int *js, threadgroup float *ws) {
  const float image = float(box.z0 + z);
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
  const int j_low = int(floor((lo + setup.span_m) / setup.step_m));
  const int j_high = int(floor((hi + setup.span_m) / setup.step_m));
  int n = 0;
  const int first = j_low < 0 ? 0 : j_low;
  const int last = j_high < setup.side - 1 ? j_high : setup.side - 1;
  for (int j = first; j <= last; ++j) {
    const float plane_low = -setup.span_m + float(j) * setup.step_m;
    const float plane_high = plane_low + setup.step_m;
    const float overlap = min(hi, plane_high) - max(lo, plane_low);
    if (overlap > 0.0f) {
      js[n] = j;
      ws[n] = overlap / span;
      ++n;
    }
  }
  return n;
}

// The same, into a thread's own arrays.
static int plane_weights_here(constant Setup &setup, FrameF frame, BoxF box, int z,
                         thread int *js, thread float *ws) {
  const float image = float(box.z0 + z);
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
  const int j_low = int(floor((lo + setup.span_m) / setup.step_m));
  const int j_high = int(floor((hi + setup.span_m) / setup.step_m));
  int n = 0;
  const int first = j_low < 0 ? 0 : j_low;
  const int last = j_high < setup.side - 1 ? j_high : setup.side - 1;
  for (int j = first; j <= last; ++j) {
    const float plane_low = -setup.span_m + float(j) * setup.step_m;
    const float plane_high = plane_low + setup.step_m;
    const float overlap = min(hi, plane_high) - max(lo, plane_low);
    if (overlap > 0.0f) {
      js[n] = j;
      ws[n] = overlap / span;
      ++n;
    }
  }
  return n;
}

static float block_sum(float value, threadgroup float *scratch, uint t) {
  scratch[t] = value;
  threadgroup_barrier(mem_flags::mem_threadgroup);
  for (uint stride = K_THREADS / 2; stride > 0; stride /= 2) {
    if (t < stride)
      scratch[t] += scratch[t + stride];
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }
  const float total = scratch[0];
  threadgroup_barrier(mem_flags::mem_threadgroup);
  return total;
}

kernel void fit_kernel(constant Setup &setup [[buffer(0)]],
                       device const PanelF *panels [[buffer(1)]],
                       device const BoxF *boxes [[buffer(2)]],
                       constant int &n_boxes [[buffer(3)]],
                       device const float *data [[buffer(4)]],
                       device const uchar *mask [[buffer(5)]],
                       device const float *reference [[buffer(6)]],
                       device float *corners [[buffer(7)]],
                       device float *profile_all [[buffer(8)]],
                       device FitF *out [[buffer(9)]],
                       uint b [[threadgroup_position_in_grid]],
                       uint t [[thread_position_in_threadgroup]]) {
  threadgroup FrameF frame;
  threadgroup int planes;
  threadgroup int js[K_MAX_SIDE];
  threadgroup float ws[K_MAX_SIDE];
  threadgroup float slice[K_MAX_SIDE * K_MAX_SIDE];
  threadgroup float scratch[K_THREADS];
  threadgroup float blend[K_BLEND];
  if (int(b) >= n_boxes)
    return;
  const BoxF box = boxes[b];
  if (t == 0)
    frame = frame_of(setup, box.s1);
  threadgroup_barrier(mem_flags::mem_threadgroup);
  const FrameF f = frame;
  if (f.valid == 0) {
    if (t == 0)
      out[b] = FitF{0, 0, 0.0f, 0.0f, 0.0f, 1.0f};
    return;
  }
  const PanelF p = panels[box.panel];
  const int cw = box.nx + 1, cn = cw * (box.ny + 1);
  device float *c1 = corners + box.corners_at;
  device float *c2 = c1 + cn;
  for (int n = int(t); n < cn; n += K_THREADS) {
    float e1, e2;
    epsilon12(f, p, float(box.x0 + n % cw), float(box.y0 + n / cw), e1, e2);
    c1[n] = e1;
    c2[n] = e2;
  }
  threadgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup);

  const float share = 1.0f / float(setup.sub * setup.sub);
  const int pixels = box.nx * box.ny, slice_size = setup.side * setup.side;
  device const float *ref = reference + box.reference_at;
  device float *profile = profile_all + box.voxels_at;
  if (box.nz * slice_size <= K_BLEND) {
    // Every slice's blended reference at once, each element with its own
    // slice's plane weights -- no one thread working while the rest wait --
    // and then each pixel's cells found once and walked down all its slices.
    // The same sums in the same order as a slice at a time: the same values.
    for (int e = int(t); e < box.nz * slice_size; e += K_THREADS) {
      const int z = e / slice_size, q = e % slice_size;
      int my_js[K_MAX_SIDE];
      float my_ws[K_MAX_SIDE];
      const int np = plane_weights_here(setup, f, box, z, my_js, my_ws);
      float value = 0.0f;
      for (int k = 0; k < np; ++k)
        value += my_ws[k] * ref[my_js[k] * slice_size + q];
      blend[e] = value;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (int k = int(t); k < pixels; k += K_THREADS) {
      ushort ids[K_MAX_SUB * K_MAX_SUB], counts[K_MAX_SUB * K_MAX_SUB];
      const int found = pixel_cells(setup, box, c1, c2, k % box.nx, k / box.nx, ids, counts);
      for (int z = 0; z < box.nz; ++z) {
        threadgroup const float *one = blend + z * slice_size;
        float value = 0.0f;
        for (int q = 0; q < found; ++q)
          value += float(counts[q]) * one[ids[q]];
        profile[z * pixels + k] = value * share;
      }
    }
    threadgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup);
  } else {
    for (int z = 0; z < box.nz; ++z) {
      if (t == 0)
        planes = plane_weights(setup, f, box, z, js, ws);
      threadgroup_barrier(mem_flags::mem_threadgroup);
      const int np = planes;
      if (np > 0)
        for (int q = int(t); q < slice_size; q += K_THREADS) {
          float value = 0.0f;
          for (int k = 0; k < np; ++k)
            value += ws[k] * ref[js[k] * slice_size + q];
          slice[q] = value;
        }
      threadgroup_barrier(mem_flags::mem_threadgroup);
      for (int k = int(t); k < pixels; k += K_THREADS) {
        float value = 0.0f;
        if (np > 0) {
          ushort ids[K_MAX_SUB * K_MAX_SUB], counts[K_MAX_SUB * K_MAX_SUB];
          const int found = pixel_cells(setup, box, c1, c2, k % box.nx, k / box.nx, ids, counts);
          for (int q = 0; q < found; ++q)
            value += float(counts[q]) * slice[ids[q]];
          value *= share;
        }
        profile[z * pixels + k] = value;
      }
      threadgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup);
    }
  }

  const int voxels = pixels * box.nz;
  device const float *d = data + box.voxels_at;
  device const uchar *m = mask + box.voxels_at;
  const float bg = box.background;
  float mine_sum = 0.0f, mine_seen = 0.0f, mine_scale = 0.0f;
  for (int i = int(t); i < voxels; i += K_THREADS) {
    if ((m[i] & kForeground) == 0)
      continue;
    mine_sum += profile[i];
    if ((m[i] & kValid) != 0) {
      mine_seen += profile[i];
      mine_scale += d[i] - bg;
    }
  }
  const float profile_sum = block_sum(mine_sum, scratch, t);
  const float seen = block_sum(mine_seen, scratch, t);
  float fit_scale = block_sum(mine_scale, scratch, t);
  if (!(profile_sum > 0.0f) || !(seen > 0.0f)) {
    if (t == 0)
      out[b] = FitF{0, 0, 0.0f, 0.0f, 0.0f, 1.0f};
    return;
  }
  fit_scale /= profile_sum;
  float variance = 0.0f;
  int iterations = 0;
  const int rounds = setup.iterations > 1 ? setup.iterations : 1;
  for (int round = 0; round < rounds; ++round) {
    float num = 0.0f, den = 0.0f;
    for (int i = int(t); i < voxels; i += K_THREADS) {
      if ((m[i] & kForeground) == 0 || (m[i] & kValid) == 0)
        continue;
      const float expected = bg + max(fit_scale, 0.0f) * profile[i];
      const float v = setup.gain * max(expected, kLeastVariance);
      num += profile[i] * (d[i] - bg) / v;
      den += profile[i] * profile[i] / v;
    }
    const float numerator = block_sum(num, scratch, t);
    const float denominator = block_sum(den, scratch, t);
    if (!(denominator > 0.0f)) {
      if (t == 0)
        out[b] = FitF{0, 0, 0.0f, 0.0f, 0.0f, 1.0f};
      return;
    }
    fit_scale = numerator / denominator;
    variance = 1.0f / denominator;
    iterations = round + 1;
  }
  float fg = 0.0f, bgp = 0.0f, n = 0.0f, sp = 0.0f, sd = 0.0f, spp = 0.0f, sdd = 0.0f, spd = 0.0f;
  for (int i = int(t); i < voxels; i += K_THREADS) {
    if ((m[i] & kValid) == 0)
      continue;
    if ((m[i] & kForeground) != 0) {
      fg += 1.0f;
      const float a = profile[i], x = d[i] - bg;
      n += 1.0f;
      sp += a;
      sd += x;
      spp += a * a;
      sdd += x * x;
      spd += a * x;
    } else if ((m[i] & kBackground) != 0) {
      bgp += 1.0f;
    }
  }
  const float foreground = block_sum(fg, scratch, t);
  const float background_pixels = block_sum(bgp, scratch, t);
  const float N = block_sum(n, scratch, t), SP = block_sum(sp, scratch, t);
  const float SD = block_sum(sd, scratch, t), SPP = block_sum(spp, scratch, t);
  const float SDD = block_sum(sdd, scratch, t), SPD = block_sum(spd, scratch, t);
  if (t == 0) {
    FitF r;
    r.valid = 1;
    r.iterations = iterations;
    r.measured = seen / profile_sum;
    r.intensity = fit_scale * profile_sum;
    r.variance = variance * profile_sum * profile_sum;
    if (background_pixels > 0.0f)
      r.variance += setup.gain * (foreground / background_pixels) * (foreground * bg);
    r.correlation = 0.0f;
    if (N > 1.0f) {
      const float cov = SPD / N - (SP / N) * (SD / N);
      const float va = SPP / N - (SP / N) * (SP / N), vb = SDD / N - (SD / N) * (SD / N);
      if (va > 0.0f && vb > 0.0f)
        r.correlation = cov / sqrt(va * vb);
    }
    out[b] = r;
  }
}
