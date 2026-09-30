// Profile fitting on a CUDA device: the batch's arrays to fit_cuda_submit and
// fit_cuda_collect, in fit_cuda.cu, which is all nvcc compiles.

#include "fit_batch.hh"

namespace fitdev {
bool fit_cuda_submit(int slot, const Setup &setup, const PanelF *panels,
                     int n_panels, const BoxF *boxes, int n_boxes,
                     const float *data, const std::uint8_t *mask,
                     std::size_t voxels, const float *reference,
                     std::size_t reference_floats, std::size_t corner_floats);
bool fit_cuda_collect(int slot, FitF *out);
const char *fit_cuda_name();
void fit_cuda_times(double out[3]);
} // namespace fitdev

namespace mxi {

namespace {
int g_next_ticket = 0;
std::size_t g_boxes[kFitSlots] = {0, 0};
} // namespace

int fit_batch_submit(const FitBatch &batch) {
  if (batch.boxes.empty() || fit_device_name() == nullptr)
    return -1;
  const int ticket = g_next_ticket;
  if (!fitdev::fit_cuda_submit(
          ticket % kFitSlots, batch.setup, batch.panels.data(),
          static_cast<int>(batch.panels.size()), batch.boxes.data(),
          static_cast<int>(batch.boxes.size()), batch.data.data(),
          batch.mask.data(), batch.voxels, batch.reference.data(),
          batch.reference.size(), batch.corner_floats))
    return -1;
  g_boxes[ticket % kFitSlots] = batch.boxes.size();
  ++g_next_ticket;
  return ticket;
}

bool fit_batch_collect(int ticket, std::vector<fitdev::FitF> *out) {
  out->assign(g_boxes[ticket % kFitSlots], fitdev::FitF{});
  return fitdev::fit_cuda_collect(ticket % kFitSlots, out->data());
}

bool fit_batch_device(const FitBatch &batch, std::vector<fitdev::FitF> *out) {
  if (batch.boxes.empty()) {
    out->clear();
    return true;
  }
  const int ticket = fit_batch_submit(batch);
  return ticket >= 0 && fit_batch_collect(ticket, out);
}

const char *fit_device_name() { return fitdev::fit_cuda_name(); }
void fit_device_times(double out[3]) { fitdev::fit_cuda_times(out); }

} // namespace mxi
