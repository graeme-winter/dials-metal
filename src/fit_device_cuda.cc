// Profile fitting on a CUDA device: the batch's arrays to fit_cuda_run, in
// fit_cuda.cu, which is all nvcc compiles.

#include "fit_batch.hh"

namespace fitdev {
bool fit_cuda_run(const Setup &setup, const PanelF *panels, int n_panels,
                  const BoxF *boxes, int n_boxes, const float *data,
                  const std::uint8_t *mask, std::size_t voxels,
                  const float *reference, std::size_t reference_floats,
                  std::size_t corner_floats, FitF *out);
const char *fit_cuda_name();
void fit_cuda_times(double out[3]);
} // namespace fitdev

namespace mxi {

bool fit_batch_device(const FitBatch &batch, std::vector<fitdev::FitF> *out) {
  out->assign(batch.boxes.size(), fitdev::FitF{});
  if (batch.boxes.empty())
    return true;
  return fitdev::fit_cuda_run(
      batch.setup, batch.panels.data(), static_cast<int>(batch.panels.size()),
      batch.boxes.data(), static_cast<int>(batch.boxes.size()),
      batch.data.data(), batch.mask.data(), batch.voxels,
      batch.reference.data(), batch.reference.size(), batch.corner_floats,
      out->data());
}

const char *fit_device_name() { return fitdev::fit_cuda_name(); }
void fit_device_times(double out[3]) { fitdev::fit_cuda_times(out); }

} // namespace mxi
