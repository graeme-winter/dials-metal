// Reference profiles and profile fitting: Kabsch (2010a) sections 3.3 and 3.4,
// Leslie (1999) section 6.
//
// WHY FIT AT ALL
// --------------
// Summation adds every foreground pixel with equal weight, so a pixel at the
// edge of a spot contributes its share of the background noise and almost none
// of the signal. Leslie section 6.6: for a weak reflection, weighting pixels by
// the profile reduces the variance by a factor of
//
//     sum(P^2) * m / (sum P)^2
//
// which is about two for a typical spot, so the standard error falls by root
// two. It gains nothing for strong reflections, where the background is
// negligible, and Leslie section 6.4 shows the fitted intensity reduces exactly
// to the summed one in that limit -- which is a useful check on an
// implementation rather than a mere remark.
//
// ON THE GRID, NOT ON THE DETECTOR
// --------------------------------
// Both learning and fitting happen in the Kabsch frame, on a cube of
// (2n+1) points spanning plus and minus `half_width` sigmas. That is the whole
// point of the frame: a reflection's shape is the same there wherever it sits
// on the detector and whatever its rocking curve, so profiles from different
// reflections can be averaged, which on the detector they cannot.
//
// A pixel is not a point in that frame -- it maps to a region -- so each pixel
// is subdivided and its counts shared between the grid points its pieces land
// in. Kabsch uses five subdivisions a side.
//
// TWO PASSES OVER THE IMAGES
// --------------------------
// The reference profiles have to exist before anything can be fitted, and they
// are learned from the reflections themselves. So the images are read twice:
// once to learn, once to fit. This is the slow way and it is the right way
// round: a profile learned from part of the scan and applied to the rest would
// be a different algorithm, and one whose errors would be hard to attribute.

#pragma once

#include <cstddef>
#include <vector>

#include "geometry.hh"
#include "shoebox.hh"

namespace mxi {

struct GridSpec {
  int n = 4;               //: so the side is 2n + 1
  double sigma_d = 0.0;    //: degrees
  double sigma_m = 0.0;    //: degrees
  double half_width = 3.0; //: in sigmas
  int subdivisions = 5;    //: per detector axis, as Kabsch uses

  int side() const { return 2 * n + 1; }
  std::size_t size() const {
    const std::size_t s = static_cast<std::size_t>(side());
    return s * s * s;
  }
  std::size_t at(int i1, int i2, int i3) const {
    const std::size_t s = static_cast<std::size_t>(side());
    return (static_cast<std::size_t>(i3) * s + static_cast<std::size_t>(i2)) *
               s +
           static_cast<std::size_t>(i1);
  }
};

//: One shoebox in the Kabsch frame: counts and background on the grid.
//:
//: Both, because fitting needs the background at each grid point and not only
//: the total -- the fit is of a profile plus a background to the raw counts,
//: and a grid point's background depends on how many pixels landed in it.
struct Transformed {
  bool valid = false;
  std::vector<double> data;       //: counts
  std::vector<double> background; //: the fitted background, same mapping
  //: How much of each pixel's area reached the grid, per point. A grid point
  //: no pixel reached is not a zero measurement, it is no measurement.
  std::vector<double> coverage;
  //: The fraction of the shoebox's counts that fell outside the grid. A
  //: reflection that mostly misses is not one to learn from.
  double outside = 0.0;
};

Transformed transform_shoebox(const Experiment &e, const Shoebox &box,
                              const Vec3 &s1, double phi_calculated,
                              const GridSpec &spec);

//: The same thing done the obvious way: every subdivision of every pixel
//: through epsilon_of, every valid voxel through every subdivision and plane.
//: transform_shoebox is tested against this, which is its specification.
Transformed transform_shoebox_direct(const Experiment &e, const Shoebox &box,
                                     const Vec3 &s1, double phi_calculated,
                                     const GridSpec &spec);

//: Reference profiles, over the detector AND over the scan.
//:
//: Kabsch section 3.3 and Leslie section 6.1 both divide the detector: the
//: profile changes across its face through obliquity of incidence, the
//: projected diffracting volume and absorption in the sensor. MOSFLM uses nine
//: or twenty-five regions.
//:
//: The profile also changes along the SCAN, because the crystal does -- it is
//: refined scan-varying for exactly that reason, and on a long sweep the
//: sample itself changes. Measured here: with one profile per region for a
//: whole scan, the fitted intensities drift from 0.986 of DIALS' at the start
//: to 0.951 at the end, which is the shape of a profile that fits the middle
//: of a scan and neither end. So the scan is divided too.
//:
//: And the profile used at a reflection is a WEIGHTED AVERAGE of the nearby
//: ones rather than the nearest, with weights falling linearly with distance,
//: as Leslie section 6.1 describes. Taking the nearest makes the model jump at
//: a region boundary, so two neighbouring reflections either side of one are
//: fitted with different profiles and their intensities differ by more than
//: their positions warrant.
struct ReferenceProfiles {
  GridSpec spec;
  int divisions = 3; //: `divisions * divisions` regions across a panel
  int blocks = 1;    //: divisions along the scan
  std::size_t panels = 1;
  //: The range of images the blocks divide, which is the range being
  //: INTEGRATED and not the whole scan.
  //:
  //: Dividing the scan instead leaves every block outside the range empty: a
  //: slice of 180 frames of an 1800 frame scan, cut into eighteen blocks, put
  //: everything in two of them and had the other sixteen borrow the detector
  //: average -- 144 of 162 cells, which is not a scan-varying profile model at
  //: all.
  double first_image = 0.0;
  double last_image = 1.0;
  //: `profile[region]` is one normalised profile, summing to one.
  std::vector<std::vector<double>> profile;
  std::vector<std::size_t> spots;
  bool finalised = false;

  std::size_t regions_per_panel() const {
    return static_cast<std::size_t>(divisions) *
           static_cast<std::size_t>(divisions) *
           static_cast<std::size_t>(blocks);
  }
  std::size_t region_count() const { return panels * regions_per_panel(); }
  //: The index of one cell, by its position in the three divisions.
  std::size_t index_of(std::size_t which_panel, int block, int j, int i) const;
  //: Which cell a reflection falls in, for LEARNING: a contribution goes to
  //: one cell, so that the profiles stay independent estimates.
  std::size_t region_of(const Panel &panel, std::size_t which_panel,
                        double px_fast, double px_slow, double z) const;
};

struct Neighbour {
  std::size_t region;
  double weight;
};

//: The cells near a reflection and how much each counts, for FITTING.
//:
//: Trilinear in the three divisions, so the profile varies smoothly across the
//: detector and along the scan instead of jumping at a boundary. At an edge
//: the weights fall back onto the cells that exist rather than reaching for
//: ones that do not.
std::vector<Neighbour> neighbours_of(const ReferenceProfiles &reference,
                                     const Panel &panel,
                                     std::size_t which_panel, double px_fast,
                                     double px_slow, double z);

//: The interpolated profile at a reflection, normalised to unit sum.
std::vector<double> profile_at(const ReferenceProfiles &reference,
                               const Panel &panel, std::size_t which_panel,
                               double px_fast, double px_slow, double z);

ReferenceProfiles make_reference(const GridSpec &spec, int divisions,
                                 int blocks, std::size_t panels,
                                 double first_image, double last_image);

//: Add one reflection's transformed shoebox to its region's profile.
//:
//: The contribution is normalised to unit sum first, so that a strong
//: reflection does not simply outvote everything else: the profile wanted is
//: the average SHAPE, not the average spot.
bool add_reference(ReferenceProfiles *reference, std::size_t region,
                   const Transformed &t);

//: Normalise every region's profile. A region with too few spots borrows the
//: average of its OWN SCAN BLOCK across the whole detector, rather than being
//: left empty: the profile drifts along the scan, which is why there are blocks
//: at all, so the whole scan's average is the wrong stand-in -- and a block's
//: average is known as soon as the block's reflections have all been seen,
//: which lets integration fit in one pass over the images. A block with no
//: spots at all borrows the nearest earlier block's average, or failing that
//: the nearest later one's.
void finalise_reference(ReferenceProfiles *reference, std::size_t least = 10);

//: The scan block a cell is in.
int block_of_cell(const ReferenceProfiles &reference, std::size_t cell);

//: The pieces finalise_reference is made of, for finalising a block at a time:
//: the average over one block's cells across the detector, normalised, before
//: any of them is normalised (false if the block saw no spots); and one
//: block's cells finalised against the average it borrows from.
bool block_average(const ReferenceProfiles &reference, int block,
                   std::vector<double> *average);
void finalise_block(ReferenceProfiles *reference, int block, std::size_t least,
                    const std::vector<double> &average);

//: The reference profile evaluated at each voxel of a shoebox.
//:
//: The transpose of `transform_shoebox`: that one carries counts from pixels
//: onto the grid, this one carries the profile from the grid back onto the
//: pixels, through the same subdivisions and the same eps3 overlaps.
//:
//: It exists because the variance of a fit is only right if it is computed
//: over the independent measurements, and those are the pixels. Fitting on the
//: grid treats its points as independent when one pixel's counts are spread
//: over several of them, which overcounts the information: measured against
//: DIALS, fitting on the grid claimed a variance 0.35 of the summed one at
//: high resolution where DIALS has 0.84, so the errors were too small by a
//: factor of 1.6 and everything downstream weighted by them was wrong.
std::vector<double> profile_on_pixels(const Experiment &e, const Shoebox &box,
                                      const Vec3 &s1, double phi_calculated,
                                      const GridSpec &spec,
                                      const std::vector<double> &reference);

//: The same thing done the obvious way: every subdivision of every pixel
//: through epsilon_of, every voxel through every subdivision and plane.
//: profile_on_pixels is tested against this, which is its specification, and
//: is about an order of magnitude faster.
std::vector<double>
profile_on_pixels_direct(const Experiment &e, const Shoebox &box,
                         const Vec3 &s1, double phi_calculated,
                         const GridSpec &spec,
                         const std::vector<double> &reference);

struct ProfileFit {
  bool valid = false;
  double intensity = 0.0;
  double variance = 0.0;
  //: Pearson correlation between the reference profile and this reflection's
  //: background-subtracted grid. DIALS writes it as `profile.correlation` and
  //: it is the honest way to spot a reflection the profile does not describe.
  double correlation = 0.0;
  int iterations = 0;
  //: The fraction of the profile that was actually measured, one when nothing
  //: was masked. A fit is an extrapolation below that and the caller decides
  //: how far it will go: a reflection with a tenth of itself visible has an
  //: intensity, and not one anybody should merge.
  double measured = 1.0;
};

//: Fit `reference` to `t`, returning the scaled intensity and its variance.
//:
//: Minimises sum_j w_j (K P_j + B_j - D_j)^2 over the scale K alone, with the
//: background held at what the GLM found. The weights are Poisson,
//: w_j = 1 / v_j, and v_j depends on K through v = B + K P -- so it is
//: iterated, as Kabsch does. Two or three rounds is plenty.
ProfileFit fit_profile(const std::vector<double> &reference,
                       const Transformed &t, double gain = 1.0,
                       int iterations = 3);

//: Fit `pixel_profile` to the shoebox's own pixels.
//:
//: The same weighted least squares as `fit_profile`, over the foreground
//: voxels rather than the grid, so the variance counts each measurement once.
ProfileFit fit_on_pixels(const Shoebox &box,
                         const std::vector<double> &pixel_profile,
                         double gain = 1.0, int iterations = 3);

} // namespace mxi
