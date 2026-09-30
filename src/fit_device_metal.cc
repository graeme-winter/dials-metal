// Profile fitting on a Metal device: the batch's arrays into buffers the GPU
// shares with the CPU -- on Apple silicon one memory, so nothing crosses a bus
// -- one threadgroup of 128 a box running fit_metal.metal, and the fits read
// back where the GPU wrote them. Two batches may be in flight, each in its own
// buffers, so that one fits while the CPU packs the next.

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
  //: One batch in flight each: its own buffers, and its command buffer until
  //: it is collected, so that one fits while the next is packed and copied.
  struct Slot {
    Shared panels, boxes, data, mask, reference, corners, profile, out;
    MTL::CommandBuffer *command = nullptr;
    std::size_t n_boxes = 0;
  };
  Slot slots[kFitSlots];
  int next_ticket = 0;
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

int fit_batch_submit(const FitBatch &batch) {
  FitMetal &m = FitMetal::instance();
  if (!m.ready() || batch.boxes.empty())
    return -1;
  const int ticket = m.next_ticket;
  FitMetal::Slot &s = m.slots[ticket % kFitSlots];
  if (s.command != nullptr)
    return -1; // not collected: the caller collects before reusing a slot
  Pool pool;
  const int n_boxes = static_cast<int>(batch.boxes.size());
  const double t0 = now();
  if (!s.panels.reserve(m.device,
                        sizeof(fitdev::PanelF) * batch.panels.size()) ||
      !s.boxes.reserve(m.device, sizeof(fitdev::BoxF) * batch.boxes.size()) ||
      !s.data.reserve(m.device, sizeof(float) * batch.voxels) ||
      !s.mask.reserve(m.device, batch.voxels) ||
      !s.reference.reserve(m.device, sizeof(float) * batch.reference.size()) ||
      !s.corners.reserve(m.device, sizeof(float) * batch.corner_floats) ||
      !s.profile.reserve(m.device, sizeof(float) * batch.voxels) ||
      !s.out.reserve(m.device, sizeof(fitdev::FitF) * batch.boxes.size()))
    return -1;
  std::memcpy(s.panels.contents(), batch.panels.data(),
              sizeof(fitdev::PanelF) * batch.panels.size());
  std::memcpy(s.boxes.contents(), batch.boxes.data(),
              sizeof(fitdev::BoxF) * batch.boxes.size());
  std::memcpy(s.data.contents(), batch.data.data(),
              sizeof(float) * batch.voxels);
  std::memcpy(s.mask.contents(), batch.mask.data(), batch.voxels);
  std::memcpy(s.reference.contents(), batch.reference.data(),
              sizeof(float) * batch.reference.size());
  m.seconds[0] += now() - t0;

  MTL::CommandBuffer *command = m.queue->commandBuffer();
  if (command == nullptr)
    return -1;
  MTL::ComputeCommandEncoder *encoder = command->computeCommandEncoder();
  if (encoder == nullptr)
    return -1;
  encoder->setComputePipelineState(m.pipeline);
  encoder->setBytes(&batch.setup, sizeof(fitdev::Setup), 0);
  encoder->setBuffer(s.panels.buffer, 0, 1);
  encoder->setBuffer(s.boxes.buffer, 0, 2);
  encoder->setBytes(&n_boxes, sizeof(int), 3);
  encoder->setBuffer(s.data.buffer, 0, 4);
  encoder->setBuffer(s.mask.buffer, 0, 5);
  encoder->setBuffer(s.reference.buffer, 0, 6);
  encoder->setBuffer(s.corners.buffer, 0, 7);
  encoder->setBuffer(s.profile.buffer, 0, 8);
  encoder->setBuffer(s.out.buffer, 0, 9);
  // dispatchThreadgroups, a threadgroup a box, as the kernel indexes them.
  encoder->dispatchThreadgroups(
      MTL::Size(static_cast<NS::UInteger>(n_boxes), 1, 1),
      MTL::Size(kThreads, 1, 1));
  encoder->endEncoding();
  command->commit();
  // Kept past this pool's end, until it is collected.
  s.command = command->retain();
  s.n_boxes = batch.boxes.size();
  ++m.next_ticket;
  return ticket;
}

bool fit_batch_collect(int ticket, std::vector<fitdev::FitF> *out) {
  FitMetal &m = FitMetal::instance();
  FitMetal::Slot &s = m.slots[ticket % kFitSlots];
  if (s.command == nullptr)
    return false;
  Pool pool;
  s.command->waitUntilCompleted();
  const bool done = s.command->status() == MTL::CommandBufferStatusCompleted;
  if (done)
    m.seconds[1] += s.command->GPUEndTime() - s.command->GPUStartTime();
  else
    std::fprintf(stderr, "Metal profile fitting: the command buffer failed\n");
  s.command->release();
  s.command = nullptr;
  if (!done)
    return false;
  const double t0 = now();
  out->resize(s.n_boxes);
  std::memcpy(out->data(), s.out.contents(), sizeof(fitdev::FitF) * s.n_boxes);
  m.seconds[2] += now() - t0;
  return true;
}

bool fit_batch_device(const FitBatch &batch, std::vector<fitdev::FitF> *out) {
  if (batch.boxes.empty()) {
    out->clear();
    return true;
  }
  const int ticket = fit_batch_submit(batch);
  return ticket >= 0 && fit_batch_collect(ticket, out);
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
