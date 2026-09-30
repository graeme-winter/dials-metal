// Profile fitting on a Metal device: the batch's arrays into buffers the GPU
// shares with the CPU -- on Apple silicon one memory, so nothing crosses a bus
// -- one threadgroup of 128 a box running fit_metal.metal, and the fits read
// back where the GPU wrote them.

#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>
#include <dispatch/dispatch.h>

#include <chrono>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <string>

#include "fit_batch.hh"

// The compiled shader, embedded by CMake so that the program carries its own.
extern "C" const unsigned char mxi_fit_metallib[];
extern "C" const unsigned long mxi_fit_metallib_size;

namespace mxi {

namespace {

constexpr int kThreads = 128; // fit_metal.metal's K_THREADS

class Pool {
public:
  Pool() : pool_(NS::AutoreleasePool::alloc()->init()) {}
  ~Pool() { pool_->release(); }
  Pool(const Pool &) = delete;
  Pool &operator=(const Pool &) = delete;

private:
  NS::AutoreleasePool *pool_;
};

//: A shared buffer, kept and grown between batches.
struct Shared {
  MTL::Buffer *buffer = nullptr;
  std::size_t size = 0;
  bool reserve(MTL::Device *device, std::size_t bytes) {
    if (bytes <= size && buffer != nullptr)
      return true;
    if (buffer)
      buffer->release();
    buffer = device->newBuffer(bytes > 0 ? bytes : 16,
                               MTL::ResourceStorageModeShared);
    size = buffer ? bytes : 0;
    return buffer != nullptr;
  }
  void *contents() const { return buffer->contents(); }
};

struct FitMetal {
  MTL::Device *device = nullptr;
  MTL::CommandQueue *queue = nullptr;
  MTL::ComputePipelineState *pipeline = nullptr;
  std::string name;
  Shared panels, boxes, data, mask, reference, corners, profile, out;
  double seconds[3] = {0.0, 0.0, 0.0}; // copying in, the kernel, reading back

  bool ready() const { return pipeline != nullptr; }

  void build() {
    Pool pool;
    device = MTL::CreateSystemDefaultDevice();
    if (device == nullptr)
      return;
    dispatch_data_t bytes =
        dispatch_data_create(mxi_fit_metallib, mxi_fit_metallib_size, nullptr,
                             DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    NS::Error *error = nullptr;
    MTL::Library *library = device->newLibrary(bytes, &error);
    dispatch_release(bytes);
    if (library == nullptr) {
      std::fprintf(
          stderr, "Metal profile fitting: could not load the shader library\n");
      return;
    }
    MTL::Function *function = library->newFunction(
        NS::String::string("fit_kernel", NS::UTF8StringEncoding));
    if (function == nullptr) {
      std::fprintf(stderr,
                   "Metal profile fitting: no fit_kernel in the library\n");
      library->release();
      return;
    }
    MTL::ComputePipelineState *state =
        device->newComputePipelineState(function, &error);
    function->release();
    library->release();
    if (state == nullptr) {
      std::fprintf(stderr,
                   "Metal profile fitting: could not make the pipeline\n");
      return;
    }
    if (state->maxTotalThreadsPerThreadgroup() <
        static_cast<NS::UInteger>(kThreads)) {
      std::fprintf(
          stderr,
          "Metal profile fitting: the pipeline takes fewer than %d threads\n",
          kThreads);
      state->release();
      return;
    }
    queue = device->newCommandQueue();
    if (queue == nullptr) {
      state->release();
      return;
    }
    pipeline = state;
    name = std::string(device->name()->utf8String()) + " (Metal)";
  }

  static FitMetal &instance() {
    static FitMetal it;
    static std::once_flag once;
    std::call_once(once, [] { it.build(); });
    return it;
  }
};

double now() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

} // namespace

bool fit_batch_device(const FitBatch &batch, std::vector<fitdev::FitF> *out) {
  FitMetal &m = FitMetal::instance();
  if (!m.ready())
    return false;
  out->assign(batch.boxes.size(), fitdev::FitF{});
  if (batch.boxes.empty())
    return true;
  Pool pool;
  const int n_boxes = static_cast<int>(batch.boxes.size());
  const double t0 = now();
  if (!m.panels.reserve(m.device,
                        sizeof(fitdev::PanelF) * batch.panels.size()) ||
      !m.boxes.reserve(m.device, sizeof(fitdev::BoxF) * batch.boxes.size()) ||
      !m.data.reserve(m.device, sizeof(float) * batch.voxels) ||
      !m.mask.reserve(m.device, batch.voxels) ||
      !m.reference.reserve(m.device, sizeof(float) * batch.reference.size()) ||
      !m.corners.reserve(m.device, sizeof(float) * batch.corner_floats) ||
      !m.profile.reserve(m.device, sizeof(float) * batch.voxels) ||
      !m.out.reserve(m.device, sizeof(fitdev::FitF) * batch.boxes.size()))
    return false;
  std::memcpy(m.panels.contents(), batch.panels.data(),
              sizeof(fitdev::PanelF) * batch.panels.size());
  std::memcpy(m.boxes.contents(), batch.boxes.data(),
              sizeof(fitdev::BoxF) * batch.boxes.size());
  std::memcpy(m.data.contents(), batch.data.data(),
              sizeof(float) * batch.voxels);
  std::memcpy(m.mask.contents(), batch.mask.data(), batch.voxels);
  std::memcpy(m.reference.contents(), batch.reference.data(),
              sizeof(float) * batch.reference.size());
  const double t1 = now();

  MTL::CommandBuffer *command = m.queue->commandBuffer();
  if (command == nullptr)
    return false;
  MTL::ComputeCommandEncoder *encoder = command->computeCommandEncoder();
  if (encoder == nullptr)
    return false;
  encoder->setComputePipelineState(m.pipeline);
  encoder->setBytes(&batch.setup, sizeof(fitdev::Setup), 0);
  encoder->setBuffer(m.panels.buffer, 0, 1);
  encoder->setBuffer(m.boxes.buffer, 0, 2);
  encoder->setBytes(&n_boxes, sizeof(int), 3);
  encoder->setBuffer(m.data.buffer, 0, 4);
  encoder->setBuffer(m.mask.buffer, 0, 5);
  encoder->setBuffer(m.reference.buffer, 0, 6);
  encoder->setBuffer(m.corners.buffer, 0, 7);
  encoder->setBuffer(m.profile.buffer, 0, 8);
  encoder->setBuffer(m.out.buffer, 0, 9);
  // dispatchThreadgroups, a threadgroup a box, as the kernel indexes them.
  encoder->dispatchThreadgroups(
      MTL::Size(static_cast<NS::UInteger>(n_boxes), 1, 1),
      MTL::Size(kThreads, 1, 1));
  encoder->endEncoding();
  command->commit();
  command->waitUntilCompleted();
  if (command->status() != MTL::CommandBufferStatusCompleted) {
    std::fprintf(stderr, "Metal profile fitting: the command buffer failed\n");
    return false;
  }
  const double kernel = command->GPUEndTime() - command->GPUStartTime();
  const double t2 = now();
  std::memcpy(out->data(), m.out.contents(),
              sizeof(fitdev::FitF) * batch.boxes.size());
  const double t3 = now();
  m.seconds[0] += t1 - t0;
  m.seconds[1] += kernel;
  m.seconds[2] += t3 - t2;
  return true;
}

const char *fit_device_name() {
  FitMetal &m = FitMetal::instance();
  return m.ready() ? m.name.c_str() : nullptr;
}

void fit_device_times(double out[3]) {
  FitMetal &m = FitMetal::instance();
  for (int s = 0; s < 3; ++s)
    out[s] = m.seconds[s];
}

} // namespace mxi
