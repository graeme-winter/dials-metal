// Profile fitting on a CUDA device, in single precision: the steps of
// fit_device.hh spread over a block of threads a box. fit_device.hh says what
// each step is and how it departs from the CPU's double-precision fit, and
// fit_batch_emulated() runs the same steps in order on the CPU, which is how
// their arithmetic is tested without a device.
//
// A block of 128 threads a box: the corners of its pixels in parallel; then a
// z-slice at a time, one thread the slice's plane weights, all of them its
// blended reference into shared memory, each its pixels' values; then the
// least squares, each of its sums a block reduction. A pixel's cells are worked
// out again for every slice rather than kept: that is arithmetic, which a
// device has, where keeping them is memory, which it has less of.

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>

#include <cuda_runtime.h>

#include "fit_device.hh"

namespace fitdev {

namespace {

constexpr int kThreads = 128;
// Every slice's blended reference at once, in shared memory, when it fits: 6144
// floats, 24 KB.
constexpr int kBlend = 6144;

// The device's own clock around the upload, the kernel and the download,
// summed across batches.
cudaEvent_t g_marks[4];
bool g_marks_made = false;
double g_seconds[3] = {0.0, 0.0, 0.0};

__device__ float block_sum(float value, float *scratch) {
  const int t = static_cast<int>(threadIdx.x);
  scratch[t] = value;
  __syncthreads();
  for (int stride = kThreads / 2; stride > 0; stride /= 2) {
    if (t < stride)
      scratch[t] += scratch[t + stride];
    __syncthreads();
  }
  const float total = scratch[0];
  __syncthreads();
  return total;
}

__global__ void fit_kernel(Setup setup, const PanelF *panels, const BoxF *boxes,
                           int n_boxes, const float *data,
                           const std::uint8_t *mask, const float *reference,
                           float *corners, float *profile_all, FitF *out) {
  const int b = static_cast<int>(blockIdx.x);
  if (b >= n_boxes)
    return;
  const int t = static_cast<int>(threadIdx.x);
  __shared__ FrameF frame;
  __shared__ int planes;
  __shared__ int js[kMaxSide];
  __shared__ float ws[kMaxSide];
  __shared__ float slice[kMaxSide * kMaxSide];
  __shared__ float scratch[kThreads];
  __shared__ float blend[kBlend];
  const BoxF box = boxes[b];
  if (t == 0)
    frame = frame_of(setup, box.s1);
  __syncthreads();
  if (!frame.valid) {
    if (t == 0)
      out[b] = FitF{};
    return;
  }
  const PanelF p = panels[box.panel];
  const int cw = box.nx + 1, cn = cw * (box.ny + 1);
  float *c1 = corners + box.corners_at, *c2 = c1 + cn;
  for (int n = t; n < cn; n += kThreads)
    corner(frame, p, box, n % cw, n / cw, c1, c2);
  __syncthreads();
  const float share = 1.0f / static_cast<float>(setup.sub * setup.sub);
  const int pixels = box.nx * box.ny, slice_size = setup.side * setup.side;
  const float *ref = reference + box.reference_at;
  float *profile = profile_all + box.voxels_at;
  if (box.nz * slice_size <= kBlend) {
    // Every slice's blended reference at once, each element with its own
    // slice's plane weights, and then each pixel's cells found once and walked
    // down all its slices: the same sums in the same order, the same values.
    for (int e = t; e < box.nz * slice_size; e += kThreads) {
      const int z = e / slice_size, q = e % slice_size;
      int my_js[kMaxSide];
      float my_ws[kMaxSide];
      const int np = plane_weights(setup, frame, box, z, my_js, my_ws);
      blend[e] = np > 0 ? blended(setup, ref, my_js, my_ws, np, q) : 0.0f;
    }
    __syncthreads();
    for (int k = t; k < pixels; k += kThreads) {
      std::uint16_t ids[kMaxSub * kMaxSub], counts[kMaxSub * kMaxSub];
      const int found =
          pixel_cells(setup, box, c1, c2, k % box.nx, k / box.nx, ids, counts);
      for (int z = 0; z < box.nz; ++z)
        profile[z * pixels + k] =
            pixel_value(ids, counts, found, blend + z * slice_size, share);
    }
    __syncthreads();
  } else {
    for (int z = 0; z < box.nz; ++z) {
      if (t == 0)
        planes = plane_weights(setup, frame, box, z, js, ws);
      __syncthreads();
      if (planes > 0)
        for (int q = t; q < slice_size; q += kThreads)
          slice[q] = blended(setup, ref, js, ws, planes, q);
      __syncthreads();
      for (int k = t; k < pixels; k += kThreads) {
        float value = 0.0f;
        if (planes > 0) {
          std::uint16_t ids[kMaxSub * kMaxSub], counts[kMaxSub * kMaxSub];
          const int found = pixel_cells(setup, box, c1, c2, k % box.nx,
                                        k / box.nx, ids, counts);
          value = pixel_value(ids, counts, found, slice, share);
        }
        profile[z * pixels + k] = value;
      }
      __syncthreads();
    }
  }
  const int voxels = pixels * box.nz;
  const float *d = data + box.voxels_at;
  const std::uint8_t *m = mask + box.voxels_at;
  const float bg = box.background;
  float mine_sum = 0.0f, mine_seen = 0.0f, mine_scale = 0.0f;
  for (int i = t; i < voxels; i += kThreads) {
    if ((m[i] & kForeground) == 0)
      continue;
    mine_sum += profile[i];
    if (m[i] & kValid) {
      mine_seen += profile[i];
      mine_scale += d[i] - bg;
    }
  }
  const float profile_sum = block_sum(mine_sum, scratch);
  const float seen = block_sum(mine_seen, scratch);
  float scale = block_sum(mine_scale, scratch);
  if (!(profile_sum > 0.0f) || !(seen > 0.0f)) {
    if (t == 0)
      out[b] = FitF{};
    return;
  }
  scale /= profile_sum;
  float variance = 0.0f;
  int iterations = 0;
  const int rounds = setup.iterations > 1 ? setup.iterations : 1;
  for (int round = 0; round < rounds; ++round) {
    float num = 0.0f, den = 0.0f;
    for (int i = t; i < voxels; i += kThreads) {
      if (!(m[i] & kForeground) || !(m[i] & kValid))
        continue;
      const float expected = bg + (scale > 0.0f ? scale : 0.0f) * profile[i];
      const float v =
          setup.gain * (expected > kLeastVariance ? expected : kLeastVariance);
      num += profile[i] * (d[i] - bg) / v;
      den += profile[i] * profile[i] / v;
    }
    const float numerator = block_sum(num, scratch),
                denominator = block_sum(den, scratch);
    if (!(denominator > 0.0f)) {
      if (t == 0)
        out[b] = FitF{};
      return;
    }
    scale = numerator / denominator;
    variance = 1.0f / denominator;
    iterations = round + 1;
  }
  float fg = 0.0f, bgp = 0.0f, n = 0.0f, sp = 0.0f, sd = 0.0f, spp = 0.0f,
        sdd = 0.0f, spd = 0.0f;
  for (int i = t; i < voxels; i += kThreads) {
    if ((m[i] & kValid) == 0)
      continue;
    if (m[i] & kForeground) {
      fg += 1.0f;
      const float a = profile[i], x = d[i] - bg;
      n += 1.0f;
      sp += a;
      sd += x;
      spp += a * a;
      sdd += x * x;
      spd += a * x;
    } else if (m[i] & kBackground) {
      bgp += 1.0f;
    }
  }
  const float foreground = block_sum(fg, scratch),
              background_pixels = block_sum(bgp, scratch);
  const float N = block_sum(n, scratch), SP = block_sum(sp, scratch),
              SD = block_sum(sd, scratch);
  const float SPP = block_sum(spp, scratch), SDD = block_sum(sdd, scratch),
              SPD = block_sum(spd, scratch);
  if (t == 0) {
    FitF f;
    f.measured = seen / profile_sum;
    f.iterations = iterations;
    f.intensity = scale * profile_sum;
    f.variance = variance * profile_sum * profile_sum;
    if (background_pixels > 0.0f)
      f.variance +=
          setup.gain * (foreground / background_pixels) * (foreground * bg);
    if (N > 1.0f) {
      const float cov = SPD / N - (SP / N) * (SD / N);
      const float va = SPP / N - (SP / N) * (SP / N),
                  vb = SDD / N - (SD / N) * (SD / N);
      if (va > 0.0f && vb > 0.0f)
        f.correlation = cov / sqrtf(va * vb);
    }
    f.valid = 1;
    out[b] = f;
  }
}

//: Page-locked host memory, which the device reads at the link's full speed,
//: kept and grown as a Buffer is. A copy from ordinary memory is staged by the
//: driver through its own pinned buffers, in pieces, well below that.
struct Pinned {
  void *p = nullptr;
  std::size_t size = 0;
  bool reserve(std::size_t bytes) {
    if (bytes <= size)
      return true;
    if (p)
      cudaFreeHost(p);
    p = nullptr;
    size = 0;
    if (cudaMallocHost(&p, bytes > 0 ? bytes : 1) != cudaSuccess)
      return false;
    size = bytes;
    return true;
  }
};

struct Buffer {
  void *p = nullptr;
  std::size_t size = 0;
  bool reserve(std::size_t bytes) {
    if (bytes <= size)
      return true;
    if (p)
      cudaFree(p);
    p = nullptr;
    size = 0;
    if (cudaMalloc(&p, bytes > 0 ? bytes : 1) != cudaSuccess)
      return false;
    size = bytes;
    return true;
  }
};

bool ok(cudaError_t e, const char *what) {
  if (e == cudaSuccess)
    return true;
  std::fprintf(stderr, "CUDA profile fitting: %s: %s\n", what,
               cudaGetErrorString(e));
  return false;
}

} // namespace

bool fit_cuda_run(const Setup &setup, const PanelF *panels, int n_panels,
                  const BoxF *boxes, int n_boxes, const float *data,
                  const std::uint8_t *mask, std::size_t voxels,
                  const float *reference, std::size_t reference_floats,
                  std::size_t corner_floats, FitF *out) {
  static Buffer d_panels, d_boxes, d_data, d_mask, d_reference, d_corners,
      d_profile, d_out;
  if (!d_panels.reserve(sizeof(PanelF) * static_cast<std::size_t>(n_panels)) ||
      !d_boxes.reserve(sizeof(BoxF) * static_cast<std::size_t>(n_boxes)) ||
      !d_data.reserve(sizeof(float) * voxels) || !d_mask.reserve(voxels) ||
      !d_reference.reserve(sizeof(float) * reference_floats) ||
      !d_corners.reserve(sizeof(float) * corner_floats) ||
      !d_profile.reserve(sizeof(float) * voxels) ||
      !d_out.reserve(sizeof(FitF) * static_cast<std::size_t>(n_boxes)))
    return ok(cudaErrorMemoryAllocation, "allocating device memory");
  if (!g_marks_made) {
    for (cudaEvent_t &e : g_marks)
      cudaEventCreate(&e);
    g_marks_made = true;
  }
  cudaEventRecord(g_marks[0]);
  // The pixels, masks and references -- nearly all of what goes -- through
  // pinned memory: a host copy into it, then the device reads it at the link's
  // speed.
  static Pinned h_data, h_mask, h_reference;
  if (!h_data.reserve(sizeof(float) * voxels) || !h_mask.reserve(voxels) ||
      !h_reference.reserve(sizeof(float) * reference_floats))
    return ok(cudaErrorMemoryAllocation, "allocating pinned memory");
  std::memcpy(h_data.p, data, sizeof(float) * voxels);
  std::memcpy(h_mask.p, mask, voxels);
  std::memcpy(h_reference.p, reference, sizeof(float) * reference_floats);
  if (!ok(cudaMemcpy(d_panels.p, panels, sizeof(PanelF) * n_panels,
                     cudaMemcpyHostToDevice),
          "copying panels") ||
      !ok(cudaMemcpy(d_boxes.p, boxes, sizeof(BoxF) * n_boxes,
                     cudaMemcpyHostToDevice),
          "copying boxes") ||
      !ok(cudaMemcpyAsync(d_data.p, h_data.p, sizeof(float) * voxels,
                          cudaMemcpyHostToDevice),
          "copying pixels") ||
      !ok(cudaMemcpyAsync(d_mask.p, h_mask.p, voxels, cudaMemcpyHostToDevice),
          "copying masks") ||
      !ok(cudaMemcpyAsync(d_reference.p, h_reference.p,
                          sizeof(float) * reference_floats,
                          cudaMemcpyHostToDevice),
          "copying references"))
    return false;
  cudaEventRecord(g_marks[1]);
  fit_kernel<<<n_boxes, kThreads>>>(
      setup, static_cast<const PanelF *>(d_panels.p),
      static_cast<const BoxF *>(d_boxes.p), n_boxes,
      static_cast<const float *>(d_data.p),
      static_cast<const std::uint8_t *>(d_mask.p),
      static_cast<const float *>(d_reference.p),
      static_cast<float *>(d_corners.p), static_cast<float *>(d_profile.p),
      static_cast<FitF *>(d_out.p));
  if (!ok(cudaGetLastError(), "launching the fit"))
    return false;
  cudaEventRecord(g_marks[2]);
  const bool copied = ok(
      cudaMemcpy(out, d_out.p, sizeof(FitF) * n_boxes, cudaMemcpyDeviceToHost),
      "copying the fits back");
  cudaEventRecord(g_marks[3]);
  cudaEventSynchronize(g_marks[3]);
  for (int s = 0; s < 3; ++s) {
    float ms = 0.0f;
    if (cudaEventElapsedTime(&ms, g_marks[s], g_marks[s + 1]) == cudaSuccess)
      g_seconds[s] += 1e-3 * static_cast<double>(ms);
  }
  return copied;
}

void fit_cuda_times(double out[3]) {
  for (int s = 0; s < 3; ++s)
    out[s] = g_seconds[s];
}

const char *fit_cuda_name() {
  static char name[256] = {0};
  static int state = 0; // 0 unknown, 1 present, 2 absent
  if (state == 0) {
    int count = 0;
    cudaDeviceProp prop;
    if (cudaGetDeviceCount(&count) == cudaSuccess && count > 0 &&
        cudaGetDeviceProperties(&prop, 0) == cudaSuccess) {
      std::snprintf(name, sizeof name, "%s (CUDA)", prop.name);
      state = 1;
    } else {
      state = 2;
    }
  }
  return state == 1 ? name : nullptr;
}

} // namespace fitdev
