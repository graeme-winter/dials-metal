#pragma once

// An experiment list from an NXmx HDF5 master file, as dials.import writes it:
// the beam, the detector's panels, the goniometer, the scan and the image set,
// read from the NeXus groups and their NXtransformations chains, in DIALS's
// (imgCIF) frame. What the file gets wrong or leaves out, the overrides say.

#include <array>
#include <optional>
#include <string>
#include <vector>

#include "json.hh"

namespace mxi {

struct ImportOverrides {
  std::optional<double> wavelength; //: A
  std::optional<double> distance;   //: mm, along the panel's normal
  std::optional<std::array<double, 2>> beam_centre; //: pixels, fast and slow
  std::optional<double> mu;          //: 1/mm, the sensor's attenuation
  std::optional<double> trusted_max; //: the top of the trusted range
  std::optional<std::array<int, 2>> image_range; //: first and last, from 1
};

//: The experiment list, and notes on what was assumed or overridden -- each a
//: line a person should read.
json::Value import_nxmx(const std::string &master,
                        const ImportOverrides &overrides,
                        std::vector<std::string> *notes);

//: A sensor's linear attenuation coefficient, 1/mm, at a wavelength in A: the
//: NIST (Hubbell and Seltzer) mass attenuation coefficients, interpolated
//: log-log in energy, times the material's density. 0 for a material not
//: tabulated.
double attenuation_coefficient(const std::string &material, double wavelength);

} // namespace mxi
