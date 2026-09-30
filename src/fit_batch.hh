#pragma once

// Profile fitting in single precision for many boxes at once, for --gpu: the
// boxes packed into flat arrays a device can take, fitted there -- or, where
// there is no device, by the same single-precision functions on the CPU
// (fit_device.hh) -- and the fits handed back in order.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "fit_device.hh"
#include "geometry.hh"
#include "reference.hh"
#include "shoebox.hh"

namespace mxi {

struct FitBatch {
  fitdev::Setup setup;
  std::vector<fitdev::PanelF> panels;
  std::vector<fitdev::BoxF> boxes;
  std::vector<float> data;        //: every box's voxels, in turn
  std::vector<std::uint8_t> mask; //: the same
  std::vector<float> reference;   //: every box's local reference, side^3 each
  std::size_t voxels = 0, corner_floats = 0; //: the scratch a batch needs
};

//: The setup and panels for an experiment and grid, which every batch shares.
FitBatch make_fit_batch(const Experiment &e, const GridSpec &spec, double gain);

//: Adds a box: its pixels and mask, its reflection's s1 and phi, its constant
//: background, and its local reference from profile_at.
void add_to_batch(FitBatch *batch, const Shoebox &box, const Vec3 &s1,
                  double phi, double background,
                  const std::vector<double> &local_reference);

//: Empties a batch of its boxes, keeping its setup.
void clear_batch(FitBatch *batch);

//: The fits, in the order the boxes were added, on the CPU in single
//: precision: fit_box_emulated, in parallel over the boxes.
std::vector<fitdev::FitF> fit_batch_emulated(const FitBatch &batch);

//: The fits on the device, if the build has one and the machine too; false
//: and nothing done otherwise.
bool fit_batch_device(const FitBatch &batch, std::vector<fitdev::FitF> *out);

//: The name of the device the fits would run on, or nullptr if none.
const char *fit_device_name();

//: A fit as ProfileFit, to record as the CPU's is recorded.
ProfileFit to_profile_fit(const fitdev::FitF &f);

} // namespace mxi
