// Indexing by three-dimensional Fourier transform.
//
// The idea, which goes back to Bricogne (1986): a set of reciprocal lattice
// points is a lattice plus noise, so the function
//
//     f(x) = sum_j exp(2 pi i r_j . x)
//
// is large wherever x is a real-space lattice vector, because every term is
// then in phase. Sampling f on a grid is a Fourier transform of the delta
// functions at the r_j, and its peaks are candidate real-space basis vectors.
// Three of them that are not coplanar define a cell; the right three are the
// ones under which most reflections get integer indices.
//
// Grid sizing is the part worth getting right, and it is fixed by two
// requirements rather than chosen:
//
//   * the real-space grid must reach beyond the largest cell edge, which caps
//     the reciprocal spacing at 1 / (2 * max_cell);
//   * the reciprocal grid must hold every point out to 1 / d_min without
//     aliasing, which floors the grid size at 4 * max_cell / d_min.
//
// For insulin at 3 Angstrom with a 100 Angstrom cell that is 134, so 256 is
// the next power of two and the real-space sampling comes out at 0.78
// Angstrom. Peak positions are then refined by centroid, because 0.78 Angstrom
// is far too coarse to be a cell edge.
//
// Multiple sweeps are indexed TOGETHER. Their reciprocal lattice points are
// pooled into one crystal frame -- each sweep keeps its own goniometer, scan
// and detector, so the pooling is exact, not an approximation. That is what
// guarantees a common basis across sweeps instead of several independently
// chosen lattices that then have to be reindexed onto one another. The price
// is a single UB and perfect goniometry, which is not true of any real
// multi-sweep experiment; refinement's first job is to break that constraint.

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "expt.h"
#include "geometry.h"
#include "refl.h"

namespace mxi {

//: Where the time went, in seconds. Filled in whether or not anyone asks,
//: because the cost of a clock is nothing beside what it settles: this was
//: about to be ported to a device on the assumption that the FFT dominated.
struct IndexTiming {
  double reciprocal_points = 0.0;
  double max_cell = 0.0;
  double candidate_vectors = 0.0;  //: the whole of it
  double fft = 0.0;                //: the transform alone
  double peak_search = 0.0;        //: the modulus, the peaks, the sort
  double choose_basis = 0.0;
  double fit_and_reduce = 0.0;
  double macrocycles = 0.0;   //: all of the below together
  double subset_copy = 0.0;   //: copying the table and picking the strong half
  double refine = 0.0;        //: refinement proper
  double reassign = 0.0;      //: indices reassigned from the refined model
  //: Triples of candidate vectors actually scored, and how many were skipped
  //: as too nearly degenerate. The cost of choosing a basis is this count
  //: times the number of reflections, and the count depends on the data: the
  //: same code was 3.8 per cent of one run and 30.8 per cent of another.
  std::size_t triples_scored = 0;
  std::size_t triples_skipped = 0;
  double total = 0.0;
};

struct IndexOptions {
  // Zero means "work it out": d_min from the reflections' own resolution
  // range, max_cell from the nearest-neighbour spacing of the reciprocal
  // lattice points.
  double d_min = 0.0;
  double max_cell = 0.0;
  std::size_t grid = 0;
  //: How far from an integer a Miller index may fall and still count.
  double tolerance = 0.3;
  //: Candidate basis vectors taken from the peak list.
  std::size_t n_candidates = 30;
  //: Shortest real-space vector treated as a candidate, to keep the search
  //: away from the origin peak and its immediate neighbourhood.
  double min_cell = 3.0;
  //: Macrocycles of assign, refine, re-assign. One means assign once and
  //: stop, which is what this did originally and is not enough.
  int macrocycles = 3;
  //: Refine the model on strong reflections only, then assign to all.
  //:
  //: This is the whole point of the macrocycle. A weak, marginally indexed
  //: population is internally consistent, so it does not produce outliers --
  //: it drags the model until it fits, after which nothing about it looks
  //: anomalous and no amount of outlier rejection will find it. Measured on
  //: insulin: refining on everything puts the detector 0.27 mm further away
  //: than DIALS and the cell 0.15 per cent large; refining on reflections with
  //: at least ten signal pixels reproduces DIALS to ten microns.
  bool refine_on_strong = true;
  bool verbose = false;
};

struct IndexResult {
  Crystal crystal;
  std::size_t n_indexed = 0;
  std::size_t n_total = 0;
  double fraction_indexed() const {
    return n_total ? static_cast<double>(n_indexed) / static_cast<double>(n_total) : 0.0;
  }
  double rmsd_index = 0.0;  // RMS distance of h from the nearest integer
  std::vector<Vec3> candidates;
  double d_min = 0.0;
  double max_cell = 0.0;
  std::size_t grid = 0;
  //: Reflections used to refine the model, as opposed to indexed by it.
  std::size_t n_refined_on = 0;
  int cycles_run = 0;
  IndexTiming timing;
};

// Map every reflection into the crystal frame of its own experiment. The
// `id` column selects the experiment; a reflection whose id is out of range is
// an error rather than something to skip, because silently dropping data is
// how an indexer comes to report a confident answer from a third of the spots.
std::vector<Vec3> reciprocal_lattice_points(const ExperimentList &experiments,
                                            const Table &reflections);

// Estimate the largest real-space cell edge from the spacing of the
// reciprocal lattice points.
double estimate_max_cell(const std::vector<Vec3> &points);

// Candidate real-space basis vectors, strongest first.
std::vector<Vec3> find_candidate_vectors(const std::vector<Vec3> &points,
                                         const IndexOptions &options,
                                         double d_min, double max_cell,
                                         std::size_t grid);

// Choose three of the candidates, reduce the resulting cell, and return it.
bool choose_basis(const std::vector<Vec3> &candidates,
                  const std::vector<Vec3> &points, double tolerance,
                  Crystal *crystal, std::size_t *n_indexed);

// Replace the basis by a reduced one describing the same lattice.
Mat3 reduce_basis(const Mat3 &real_space_rows);

// The whole thing. On success the crystal is set on every experiment in the
// list and `miller_index` is written into the table.
IndexResult index(ExperimentList &experiments, Table &reflections,
                  const IndexOptions &options = {});

}  // namespace mxi
