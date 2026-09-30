// No device for profile fitting in this build: --gpu fits in single precision
// on the CPU instead, and says so.

#include "fit_batch.hh"

namespace mxi {

bool fit_batch_device(const FitBatch &, std::vector<fitdev::FitF> *) {
  return false;
}
const char *fit_device_name() { return nullptr; }
void fit_device_times(double out[3]) { out[0] = out[1] = out[2] = 0.0; }

} // namespace mxi
