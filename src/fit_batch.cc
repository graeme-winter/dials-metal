#include "fit_batch.hh"

#include <algorithm>

#include "parallel.hh"

namespace mxi {

static_assert(fitdev::kValid == shoebox_mask::kValid, "the mask's bits");
static_assert(fitdev::kBackground == shoebox_mask::kBackground,
              "the mask's bits");
static_assert(fitdev::kForeground == shoebox_mask::kForeground,
              "the mask's bits");
// fit_metal.metal declares these structs again, the shading language unable to
// include fit_device.hh: every field four bytes, in the same order, and these
// sizes are what it assumes.
static_assert(sizeof(fitdev::F3) == 12, "F3 as the shader has it");
static_assert(sizeof(fitdev::PanelF) == 72, "PanelF as the shader has it");
static_assert(sizeof(fitdev::Setup) == 68, "Setup as the shader has it");
static_assert(sizeof(fitdev::BoxF) == 60, "BoxF as the shader has it");
static_assert(sizeof(fitdev::FitF) == 24, "FitF as the shader has it");

namespace {

fitdev::F3 f3(const Vec3 &v) {
  return {static_cast<float>(v.x), static_cast<float>(v.y),
          static_cast<float>(v.z)};
}

} // namespace

FitBatch make_fit_batch(const Experiment &e, const GridSpec &spec,
                        double gain) {
  FitBatch b;
  fitdev::Setup &s = b.setup;
  s.s0 = f3(e.beam.s0());
  s.axis_lab = f3(e.goniometer.lab_axis());
  s.osc_start = static_cast<float>(Scan::radians(e.scan.osc_start));
  s.osc = static_cast<float>(Scan::radians(e.scan.osc_width));
  s.z_offset = static_cast<float>(e.scan.z_offset);
  s.side = spec.side();
  s.sub = spec.subdivisions;
  const double span_d = spec.half_width * spec.sigma_d,
               span_m = spec.half_width * spec.sigma_m;
  s.span_d = static_cast<float>(span_d);
  s.step_d = static_cast<float>(2.0 * span_d / static_cast<double>(s.side));
  s.span_m = static_cast<float>(span_m);
  s.step_m = static_cast<float>(2.0 * span_m / static_cast<double>(s.side));
  s.gain = static_cast<float>(gain);
  s.iterations = 3;
  for (std::size_t k = 0; k < e.detector.size(); ++k) {
    const Panel &p = e.detector[k];
    fitdev::PanelF f;
    f.origin = f3(p.origin);
    f.fast = f3(p.fast);
    f.slow = f3(p.slow);
    f.normal = f3(p.normal());
    f.pixel_fast = static_cast<float>(p.pixel_size[0]);
    f.pixel_slow = static_cast<float>(p.pixel_size[1]);
    f.mu = static_cast<float>(p.mu);
    f.thickness = static_cast<float>(p.thickness);
    f.parallax = p.parallax ? 1 : 0;
    f.conditional = p.parallax_conditional ? 1 : 0;
    b.panels.push_back(f);
  }
  return b;
}

void add_to_batch(FitBatch *batch, const std::vector<BatchEntry> &entries,
                  bool to_device) {
  std::size_t voxels = 0, corners = 0, reference = 0;
  for (const BatchEntry &x : entries) {
    const Shoebox &box = *x.box;
    fitdev::BoxF f;
    f.x0 = box.bbox[0];
    f.y0 = box.bbox[2];
    f.z0 = box.bbox[4];
    f.nx = box.nx();
    f.ny = box.ny();
    f.nz = box.nz();
    f.panel = box.panel;
    f.s1 = f3(x.s1);
    f.phi = static_cast<float>(x.phi);
    f.background = static_cast<float>(x.background);
    f.voxels_at = static_cast<std::uint32_t>(voxels);
    f.reference_at = static_cast<std::uint32_t>(reference);
    f.corners_at = static_cast<std::uint32_t>(corners);
    voxels += box.size();
    reference += x.reference->size();
    corners += 2 * static_cast<std::size_t>(f.nx + 1) *
               static_cast<std::size_t>(f.ny + 1);
    batch->boxes.push_back(f);
  }
  batch->voxels = voxels;
  batch->corner_floats = corners;
  batch->reference_floats = reference;
  FitDestination into;
  batch->in_place =
      to_device && fit_batch_destination(voxels, reference, &into);
  if (!batch->in_place) {
    batch->data.resize(voxels);
    batch->mask.resize(voxels);
    batch->reference.resize(reference);
    into = {batch->data.data(), batch->mask.data(), batch->reference.data()};
  }
  batch->data_at = into.data;
  batch->mask_at = into.mask;
  batch->reference_at = into.reference;
  for_each_index(entries.size(), [&](std::size_t k) {
    const BatchEntry &x = entries[k];
    const fitdev::BoxF &f = batch->boxes[k];
    std::copy(x.box->data.begin(), x.box->data.end(), into.data + f.voxels_at);
    std::copy(x.box->mask.begin(), x.box->mask.end(), into.mask + f.voxels_at);
    for (std::size_t i = 0; i < x.reference->size(); ++i)
      into.reference[f.reference_at + i] =
          static_cast<float>((*x.reference)[i]);
  });
}

void clear_batch(FitBatch *batch) {
  batch->boxes.clear();
  batch->data.clear();
  batch->mask.clear();
  batch->reference.clear();
  batch->voxels = batch->corner_floats = batch->reference_floats = 0;
  batch->data_at = nullptr;
  batch->mask_at = nullptr;
  batch->reference_at = nullptr;
  batch->in_place = false;
}

std::vector<fitdev::FitF> fit_batch_emulated(const FitBatch &batch) {
  std::vector<fitdev::FitF> out(batch.boxes.size());
  for_each_index(batch.boxes.size(), [&](std::size_t k) {
    const fitdev::BoxF &box = batch.boxes[k];
    std::vector<float> corners(2 * static_cast<std::size_t>(box.nx + 1) *
                               static_cast<std::size_t>(box.ny + 1));
    std::vector<float> profile(static_cast<std::size_t>(box.nx) * box.ny *
                               box.nz);
    std::vector<float> blend(
        static_cast<std::size_t>(box.nz) *
        static_cast<std::size_t>(batch.setup.side * batch.setup.side));
    out[k] = fitdev::fit_box_emulated(
        batch.setup, batch.panels.data(), box, batch.data_at + box.voxels_at,
        batch.mask_at + box.voxels_at, batch.reference_at + box.reference_at,
        corners.data(), profile.data(), blend.data());
  });
  return out;
}

ProfileFit to_profile_fit(const fitdev::FitF &f) {
  ProfileFit p;
  p.valid = f.valid != 0;
  p.intensity = static_cast<double>(f.intensity);
  p.variance = static_cast<double>(f.variance);
  p.correlation = static_cast<double>(f.correlation);
  p.iterations = f.iterations;
  p.measured = static_cast<double>(f.measured);
  return p;
}

// fit_batch_device and fit_device_name are in fit_device_none.cc, or with CUDA
// fit_device_cuda.cc and fit_cuda.cu: a program links one.

} // namespace mxi
