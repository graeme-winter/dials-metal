#pragma once

// Scaling: grouping observations by symmetry, the merged intensities, and the
// least-squares fit of the scaling model (Beilsten-Edmands et al. 2020, eqns
// 2 and 3).

#include <cstddef>
#include <vector>

#include "expt.hh"
#include "refl.hh"
#include "scale_model.hh"
#include "symmetry.hh"

namespace mxi {

enum class IntensityChoice { profile, summation };

struct ScaleDataOptions {
  IntensityChoice use = IntensityChoice::profile;
  //: dials.scale's default: a partial below 0.4 is too uncertain to scale up.
  double partiality_cutoff = 0.4;
  double d_min = 0.0;
};

//: The observations scaling uses, as parallel arrays.
struct ScaleData {
  //: Corrected as dials.scale corrects them before scaling: I lp / qe /
  //: partiality, and the variance by the square of the same factor.
  std::vector<double> intensity, variance;
  std::vector<ScaleObservation> observation;
  std::vector<std::size_t> group; //: index into `unique`
  std::vector<Miller> unique;     //: the symmetry-unique index of each group
  std::vector<std::size_t> row;   //: in the reflection table
  std::vector<double> d;
  std::vector<bool> outlier;
  //: Both estimates, corrected likewise, for combining them; has_sum is false
  //: where summation failed, as it does across a module gap.
  std::vector<double> prf, prf_variance, sum, sum_variance;
  std::vector<bool> has_sum;
  //: The variances as integration gave them, corrected, before an error model.
  std::vector<double> variance_before;
  std::size_t size() const { return intensity.size(); }
};

//: The observations of one sweep fit to scale: integrated by the chosen
//: method, with a positive variance, partiality at the cutoff or above, a
//: finite d beyond d_min, and not absent in the space group. Each carries its
//: position in the scan, 1/2d^2, and -- if the shape has an absorption term --
//: its harmonics at s1 and the reverse incident beam in the crystal frame.
ScaleData build_scale_data(const ExperimentList &experiments,
                           const Table &reflections, const SpaceGroup &group,
                           const ScaleModelShape &shape,
                           const ScaleDataOptions &options = {});

//: g for every observation.
std::vector<double> inverse_scales(const ScaleModel &model,
                                   const ScaleData &data);

//: <I_h> = sum w g I / sum w g^2 over each group's observations that are not
//: outliers, w = 1 / variance (the paper's eqn 3). A group with none is zero.
std::vector<double> merged_intensities(const ScaleData &data,
                                       const std::vector<double> &g);

//: A random subset of whole groups for fitting, as the paper's section 4.4:
//: groups are drawn until there are at least min_groups of them and
//: min_observations observations among them, or all are used.
std::vector<std::size_t>
select_for_fitting(const ScaleData &data, std::size_t min_groups = 2000,
                   std::size_t min_observations = 50000, unsigned seed = 1);

struct ScaleFitOptions {
  int max_iterations = 50;
  //: Weights of the restraints toward zero, sum B_i^2 and sum P_lm^2.
  double decay_restraint = 0.1;
  double absorption_restraint = 1000.0;
};

struct ScaleFitResult {
  int iterations = 0;
  bool converged = false;
  double target_start = 0.0, target_end = 0.0;
  std::size_t observations = 0;
};

//: Fit the model by Levenberg-Marquardt over `use` (observation indices; all
//: that are not outliers if empty). The residuals are sqrt(w)(I - g <I_h>),
//: and are differentiated with <I_h> held: at the least-squares <I_h> the
//: target's derivative with respect to it is zero, so the gradient held is the
//: whole gradient and the fixed point is the true minimum. After each step the
//: scale is normalised, which a common factor with <I_h> leaves unseen.
ScaleFitResult fit_scale_model(ScaleModel &model, const ScaleData &data,
                               const ScaleFitOptions &options = {},
                               const std::vector<std::size_t> &use = {});

//: Observations whose normalised deviation from their group's weighted mean,
//: excluding themselves, exceeds zmax (Evans 2006), flagged in data.outlier.
//: Every observation is retested, earlier outliers too; within a group the
//: worst is removed and the rest retested, so one cannot hide another. Returns
//: how many are flagged.
std::size_t reject_outliers(ScaleData &data, const std::vector<double> &g,
                            double zmax = 6.0);

struct ErrorModel {
  double a = 1.0, b = 0.02;
  std::size_t used = 0; //: observations it was determined from
};
//: sigma'^2 = a^2 (sigma^2 + (b I)^2), determined as dials.scale does: a from
//: the slope of the central normal probability plot, |x| < 1.5, of the
//: normalised deviations (eqn 12); b by minimising eqn 17 over logarithmically
//: spaced intensity bins; the two alternated to convergence from a = 1,
//: b = 0.02; on groups with <I> above 25 and <I / sigma^2> above 0.85.
ErrorModel refine_error_model(const ScaleData &data,
                              const std::vector<double> &g);
//: The error model applied to every observation's variance, from the variance
//: before any: applying it twice does not compound.
void apply_error_model(ScaleData &data, const ErrorModel &model);

//: I = w I_prf + (1 - w) I_sum with w = 1 / (1 + (I_sum / I_mid)^3), eqns 13
//: and 14; the variance as for fully correlated estimates, which they are,
//: being from the same pixels. I_mid of zero is profile fitting alone, and of
//: infinity summation alone.
void combine_intensities(ScaleData &data, double i_mid);
//: Rmeas with the scales given, over observations that are not outliers.
double rmeas(const ScaleData &data, const std::vector<double> &g);
//: The I_mid of lowest Rmeas, among profile alone, summation alone, and powers
//: of ten between; leaves data combined with it.
double choose_intensity_combination(ScaleData &data,
                                    const std::vector<double> &g);

} // namespace mxi
