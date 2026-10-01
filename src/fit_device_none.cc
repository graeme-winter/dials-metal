// No device for profile fitting in this build: --gpu fits on the CPU, in
// double, and says so.
//
// Except with MXI_FIT_FAKE_DEVICE=1, for testing: then this poses as a device
// of kFitSlots slots, a submitted batch fitted by the CPU emulation into its
// slot and handed back when collected, so that --gpu runs mxi_integrate's own
// path for a device -- tickets, slots reused, every fit collected before the
// table is written -- where no device is.

#include "fit_batch.hh"

#include <cstdlib>
#include <string>

namespace mxi {

namespace {

bool fake() {
  const char *asked = std::getenv("MXI_FIT_FAKE_DEVICE");
  return asked && std::string(asked) == "1";
}

struct FakeSlot {
  bool busy = false;
  std::vector<fitdev::FitF> fits;
  // Its own memory, offered to pack into, as a real device's is.
  std::vector<float> data, reference;
  std::vector<std::uint8_t> mask;
};
FakeSlot g_fake[kFitSlots];
int g_next_ticket = 0;

} // namespace

bool fit_batch_destination(std::size_t voxels, std::size_t reference_floats,
                           FitDestination *out) {
  if (!fake())
    return false;
  FakeSlot &s = g_fake[g_next_ticket % kFitSlots];
  if (s.busy)
    return false;
  s.data.resize(voxels);
  s.mask.resize(voxels);
  s.reference.resize(reference_floats);
  *out = {s.data.data(), s.mask.data(), s.reference.data()};
  return true;
}

int fit_batch_submit(const FitBatch &batch) {
  if (!fake() || batch.boxes.empty())
    return -1;
  const int ticket = g_next_ticket;
  FakeSlot &s = g_fake[ticket % kFitSlots];
  if (s.busy)
    return -1;
  s.fits = fit_batch_emulated(batch);
  s.busy = true;
  ++g_next_ticket;
  return ticket;
}

bool fit_batch_collect(int ticket, std::vector<fitdev::FitF> *out) {
  FakeSlot &s = g_fake[ticket % kFitSlots];
  if (!s.busy)
    return false;
  s.busy = false;
  *out = std::move(s.fits);
  return true;
}

bool fit_batch_device(const FitBatch &batch, std::vector<fitdev::FitF> *out) {
  const int ticket = fit_batch_submit(batch);
  return ticket >= 0 && fit_batch_collect(ticket, out);
}

const char *fit_device_name() {
  return fake() ? "a pretended device (testing)" : nullptr;
}
void fit_device_times(double out[3]) { out[0] = out[1] = out[2] = 0.0; }

} // namespace mxi
