// The experimental geometry, and the two maps that matter: a spot on the
// detector to a point in the crystal's reciprocal space, and back again.
//
// CONVENTIONS -- VALIDATED
// -----------------------
// These were beliefs until they met a real .expt and indexed.refl, at which
// point one of them turned out to be backwards and a whole correction turned
// out to be missing. `docs/conventions.md` records which, and what is still
// unvalidated. The closed-loop test could not have caught either, because an
// error present in both the forward and reverse map cancels exactly.
//
// As implemented:
//
//   s0        = -direction * (1 / wavelength). NEGATIVE. dxtbx stores the beam
//               direction pointing from the sample back towards the source, so
//               s0, which points along propagation, is its negation. Getting
//               this wrong puts the Ewald sphere on the wrong side of the
//               origin and nothing indexes at all: with the sign flipped the
//               residual against A h is 2.1 per reciprocal Angstrom against a
//               reciprocal cell edge of 0.0128, so it fails loudly.
//   px <-> mm = parallax corrected, not a multiplication by the pixel size.
//               See Panel below; it is worth 1.6 pixels at the panel edge.
//   s1        = unit vector from sample to the corrected position, over the
//               wavelength.
//   q         = s1 - s0, the scattering vector in the laboratory frame.
//   R(phi)    = setting * rotation(axis, phi) * fixed, the full goniometer
//               rotation taking crystal frame to laboratory frame.
//   r0        = R(phi)^-1 * q, the reciprocal lattice point in the crystal
//               frame, which is what indexing works on.
//   A         = the setting matrix with a*, b*, c* as its COLUMNS, so that
//               r0 = A * h with h a column of Miller indices.
//   phi(z)    = osc_start + (z - z_offset) * osc_width.
//
// The last one has a known wrinkle. DIALS anchors the scan z coordinate to
// image number one rather than to the start of the array, so z_offset is zero
// and the two agree only when first_image is one. That is the common case and
// the disagreement is a constant otherwise. z_offset is left settable rather
// than hard-coded so the choice is visible.

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "linalg.h"

namespace mxi {

struct Beam {
  Vec3 direction{0.0, 0.0, 1.0};
  double wavelength = 1.0;

  // Note the sign: see the convention block at the top of this file.
  Vec3 s0() const { return -direction.normalized() / wavelength; }
};

struct Panel {
  std::string name;
  Vec3 fast{1.0, 0.0, 0.0};
  Vec3 slow{0.0, -1.0, 0.0};
  Vec3 origin{0.0, 0.0, 100.0};
  double pixel_size[2] = {0.1, 0.1};
  std::int64_t image_size[2] = {1000, 1000};
  double trusted_min = 0.0;
  double trusted_max = 65535.0;

  // Parallax correction. An X-ray entering a thick silicon sensor at an angle
  // travels some way through it before being absorbed, so the pixel that fires
  // is displaced outwards from where the ray met the front face. On a 0.45 mm
  // sensor that reaches 1.6 pixels at the edge of the panel, which is twenty
  // times the centroid precision and impossible to absorb into the geometry.
  //
  // Ignoring it does not make the model wrong by a constant; it makes it wrong
  // by a smooth function of the scattering angle, which refinement then partly
  // absorbs into the detector distance, leaving a radial residual behind.
  bool parallax = false;
  double mu = 0.0;         // attenuation coefficient, per mm
  double thickness = 0.0;  // sensor thickness, mm

  Vec3 normal() const { return fast.cross(slow).normalized(); }

  // The lateral displacement caused by absorption at depth, evaluated at the
  // position given. Both directions of the conversion evaluate it at their own
  // input, which is not self-inverse to better than about a thousandth of a
  // millimetre -- and is exactly what DIALS does, checked against real
  // xyzobs.mm and xyzcal.px to the last bit.
  std::pair<double, double> parallax_offset(double mm_fast,
                                            double mm_slow) const;

  // Pixels to millimetres on the panel face, and back.
  std::pair<double, double> px_to_mm(double px_fast, double px_slow) const;
  std::pair<double, double> mm_to_px(double mm_fast, double mm_slow) const;

  // Laboratory position of a point given in millimetres on the panel face.
  Vec3 lab_coord_mm(double mm_fast, double mm_slow) const {
    return origin + fast * mm_fast + slow * mm_slow;
  }
  // Laboratory position of a point given in pixels, parallax included.
  Vec3 lab_coord_px(double px_fast, double px_slow) const {
    const auto mm = px_to_mm(px_fast, px_slow);
    return lab_coord_mm(mm.first, mm.second);
  }

  // Where a diffracted ray meets this panel, in pixels. Returns nothing if the
  // ray is parallel to the panel, hits it from behind, or lands off the edge.
  // Rays travelling away from the detector are rejected by the sign of the
  // path length, not by checking the panel bounds -- a ray going backwards can
  // still intersect the infinite plane inside the bounds.
  std::optional<std::pair<double, double>> intersect(const Vec3 &s1) const;
};

struct Detector {
  std::vector<Panel> panels;

  std::size_t size() const { return panels.size(); }
  const Panel &operator[](std::size_t i) const { return panels[i]; }

  // The panel a ray hits, and where, searching in order. Multi-panel detectors
  // can overlap in projection, so the first hit wins and that is deliberate:
  // it matches the order panels are serialised in.
  std::optional<std::tuple<std::size_t, double, double>> intersect(
      const Vec3 &s1) const;
};

struct Goniometer {
  Vec3 axis{1.0, 0.0, 0.0};
  Mat3 fixed = Mat3::identity();
  Mat3 setting = Mat3::identity();

  // A multi-axis goniometer, as dxtbx serialises one: a list of axes, the
  // angle each is set to, and which of them the scan turns. Everything below
  // the scan axis is carried by the sample and becomes the fixed rotation;
  // everything above it moves the scan axis itself and becomes the setting
  // rotation.
  //
  // With every angle at zero -- the common case, and the case in the insulin
  // data -- both come out as the identity, which means a single-axis dataset
  // cannot test this decomposition at all.
  static Goniometer from_axes(const std::vector<Vec3> &axes,
                              const std::vector<double> &angles_deg,
                              std::size_t scan_axis);

  Mat3 rotation_at(double phi) const {
    return setting * rotation(axis, phi) * fixed;
  }
  // The rotation axis as it appears in the laboratory, which is what
  // prediction needs -- the sample turns about the setting-rotated axis.
  Vec3 lab_axis() const { return (setting * axis).normalized(); }
};

struct Scan {
  std::int64_t first_image = 1;
  std::int64_t last_image = 1;
  double osc_start = 0.0;    // degrees
  double osc_width = 0.0;    // degrees per image
  std::int64_t batch_offset = 0;
  // See the convention note at the top of this file. Zero reproduces DIALS.
  double z_offset = 0.0;

  std::int64_t num_images() const { return last_image - first_image + 1; }

  double phi_from_z(double z) const {
    return radians(osc_start + (z - z_offset) * osc_width);
  }
  double z_from_phi(double phi) const {
    if (osc_width == 0.0) return z_offset;
    return (degrees(phi) - osc_start) / osc_width + z_offset;
  }
  double phi_start() const { return phi_from_z(0.0); }
  double phi_end() const { return phi_from_z(static_cast<double>(num_images())); }

  static double radians(double d) { return d * 3.14159265358979323846 / 180.0; }
  static double degrees(double r) { return r * 180.0 / 3.14159265358979323846; }
};

struct UnitCell {
  double a = 1, b = 1, c = 1, alpha = 90, beta = 90, gamma = 90;
  double volume() const;
};

struct Crystal {
  // Setting matrix: a*, b*, c* as columns, so r0 = A * (h, k, l).
  Mat3 A = Mat3::identity();
  std::string space_group_hall = " P 1";

  static Crystal from_real_space(const Vec3 &a, const Vec3 &b, const Vec3 &c);
  // Real-space basis vectors, the rows of A inverse.
  Vec3 real_a() const;
  Vec3 real_b() const;
  Vec3 real_c() const;
  UnitCell cell() const;
  double d_spacing(int h, int k, int l) const;
};

struct Experiment {
  Beam beam;
  Detector detector;
  Goniometer goniometer;
  Scan scan;
  std::optional<Crystal> crystal;
  std::string identifier;
};

// --------------------------------------------------------------------------
// The two maps
// --------------------------------------------------------------------------

// Laboratory scattering vector of an observed spot, before the goniometer
// rotation is undone.
Vec3 lab_scattering_vector(const Experiment &e, std::size_t panel, double px_fast,
                           double px_slow);

// Reciprocal lattice point in the crystal frame: what indexing consumes.
// `z` is the observed centroid's scan coordinate in images.
Vec3 reciprocal_lattice_point(const Experiment &e, std::size_t panel,
                              double px_fast, double px_slow, double z);

// Convenience: map a whole list of observations at once.
struct Observation {
  std::size_t panel = 0;
  double px_fast = 0.0, px_slow = 0.0, z = 0.0;
};
std::vector<Vec3> reciprocal_lattice_points(const Experiment &e,
                                            const std::vector<Observation> &obs);

}  // namespace mxi
