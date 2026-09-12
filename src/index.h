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
