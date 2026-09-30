// No device for profile fitting in this build: --gpu fits in single precision
// on the CPU instead, and says so.

#include "fit_batch.hh"

namespace mxi {

bool fit_batch_device(const FitBatch &, std::vector<fitdev::FitF> *) {
  return false;
}
const char *fit_device_name() { return nullptr; }

} // namespace mxi
