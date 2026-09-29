#pragma once

// Scaling: grouping observations by symmetry, the merged intensities, and the
// least-squares fit of the scaling model (Beilsten-Edmands et al. 2020, eqns
// 2 and 3).

#include <cstddef>
#include <string>
#include <utility>
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
  //: A byte each, not std::vector<bool>: outlier rejection writes these from
  //: several threads at once, and neighbouring bits of one word are one write.
  std::vector<std::uint8_t> outlier;
  //: Both estimates, corrected likewise, for combining them; has_sum is false
  //: where summation failed, as it does across a module gap.
  std::vector<double> prf, prf_variance, sum, sum_variance;
  std::vector<bool> has_sum;
  //: Which half of its Friedel pair each observation is (true for I+ and
  //: for centric reflections), and whether each group is centric.
  std::vector<bool> plus;
  std::vector<bool> centric;
  //: The variances as integration gave them, corrected, before an error model.
  std::vector<double> variance_before;
  //: I^2 var(g) / g^2, the uncertainty of the scale carried into the
  //: observation's own units: var(I / g) = (sigma^2 + this) / g^2. Zero until
  //: propagate_scale_variances sets it; the error model is refined and applied
  //: to sigma^2 plus this, so that a and b correct what remains after it.
  std::vector<double> scale_term;
  //: The term for observation i, zero where none was propagated -- including
  //: where the vector is shorter, as it is in data built by hand: reading it
  //: unguarded, the error model's tests ran off the end of an empty one.
  double scale_term_at(std::size_t i) const {
    return i < scale_term.size() ? scale_term[i] : 0.0;
  }
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

//: The covariance of the model's parameters at a fit, from the normal
//: equations of the variable-projection Jacobian and the restraints. Two
//: directions are fixed by the model's normalisation -- the scale's mean at one
//: and the relative B's at zero -- and the covariance is that of the parameters
//: under those constraints, Z (Z^T N Z)^-1 Z^T for Z spanning them: the first
//: of those directions is exactly null in N, and without the constraints
//: nothing can be inverted. Scaled by the goodness of fit, chi^2 over the
//: degrees of freedom, which count every merged <I> as a parameter fitted.
struct ParameterCovariance {
  std::vector<double> matrix; //: n x n, row major
  double goodness_of_fit = 1.0;
  std::size_t degrees_of_freedom = 0;
  bool ok = false;
};
ParameterCovariance parameter_covariance(const ScaleModel &model,
                                         const ScaleData &data,
                                         const ScaleFitOptions &options,
                                         const std::vector<std::size_t> &use);
//: var(g) for every observation: grad_g^T C grad_g.
std::vector<double>
inverse_scale_variances(const ScaleModel &model, const ScaleData &data,
                        const ParameterCovariance &covariance);
//: Set each observation's scale_term, I^2 var(g) / g^2.
void propagate_scale_variances(ScaleData &data, const std::vector<double> &g,
                               const std::vector<double> &g_variance);

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

//: One resolution shell of merging statistics, or all of them.
struct MergingShell {
  double d_max = 0.0, d_min = 0.0;
  std::size_t observations = 0, unique = 0, possible = 0, possible_acentric = 0;
  double multiplicity = 0.0, completeness = 0.0; //: completeness as a fraction
  double mean_i = 0.0;       //: mean over merged reflections of <I>
  double i_over_sigma = 0.0; //: mean over merged reflections of <I>/sigma(<I>)
  //: With Friedel mates merged, and with them apart ("I+/-").
  double rmerge = 0.0, rmeas = 0.0, rpim = 0.0;
  double rmerge_anom = 0.0, rmeas_anom = 0.0, rpim_anom = 0.0;
  double cc_half = 0.0;
  //: As iotbx.merging_statistics, which dials.scale reports: acentric
  //: reflections with both mates measured over the acentric ones possible; the
  //: observations over the groups with mates apart; the correlation of
  //: I+ - I- between random halves; the slope of the normal probability plot
  //: of dI / sigma(dI) over |x| < 0.9; sqrt(2 <(F+ - F-)^2> / <F+^2 + F-^2>);
  //: and mean |dI| over mean sigma(dI).
  double anom_completeness = 0.0, anom_multiplicity = 0.0, cc_anom = 0.0;
  double anom_slope = 0.0, df_over_f = 0.0, di_over_sig_di = 0.0;
  std::size_t anomalous_pairs = 0;
};

//: Everything for the groups given except what needs the cell -- possible
//: reflections, completeness, the resolution limits -- over their observations
//: that are not outliers. The random halves are drawn with a fixed seed.
MergingShell merge_groups(const ScaleData &data, const std::vector<double> &g,
                          const std::vector<std::size_t> &groups);

//: merge_groups for each resolution shell -- equal volume in reciprocal
//: space, highest resolution LAST -- and overall, with completeness: every
//: reflection the crystal's cell allows to d_min, in the asymmetric unit and
//: not absent, and the acentric ones among them.
std::vector<MergingShell> merging_statistics(const ScaleData &data,
                                             const std::vector<double> &g,
                                             const SpaceGroup &group,
                                             const Crystal &crystal, int shells,
                                             MergingShell *overall);

struct ScaleRunOptions {
  bool combine = true;    //: choose profile, summation or a mix by Rmeas
  bool absorption = true; //: when the sweep is wide enough for it
  double d_min = 0.0;
  ScaleFitOptions fit;
};

struct ScaleRun {
  ScaleData data;
  ScaleModel model{ScaleModelShape{}};
  std::vector<double> g;
  ErrorModel error_model;
  double i_mid = 0.0; //: 0 profile alone, infinity summation alone
  std::vector<ScaleFitResult> fits;
  std::vector<double> g_variance; //: var(g) for every observation
  ParameterCovariance covariance;
  std::size_t outliers = 0;
  //: Observations in the subset the model was fitted to.
  std::size_t fitted_on = 0;
  //: Wall seconds of each step, in order, for --timing.
  std::vector<std::pair<std::string, double>> timing;
};

//: Scale one sweep as the paper's figure 2: outliers on the unscaled data; fit;
//: outliers; the intensity combination; outliers; error model; fit; outliers;
//: error model; a final fit, outlier rejection and error model.
ScaleRun scale_sweep(const ExperimentList &experiments,
                     const Table &reflections, const SpaceGroup &group,
                     const ScaleRunOptions &options = {});

namespace flag {
//: dials.scale's, from dials/array_family/reflection_table.h, which is what
//: dials.merge selects by.
constexpr std::int64_t kOutlierInScaling = 1 << 23;
constexpr std::int64_t kExcludedForScaling = 1 << 24;
constexpr std::int64_t kScaled = 1 << 26;
} // namespace flag

//: The scaling written into the table as dials.scale writes it:
//: inverse_scale_factor and its variance, intensity.scale.value and variance
//: (corrected and combined, the error model's variance), and the flags --
//: scaled, an outlier in scaling, or excluded -- ADDED to those already there.
//: Rows not scaled get an inverse scale of 1 and a variance of 0. The
//: intensity's variance already carries the scale's, as dials.scale's does, so
//: nothing downstream should add inverse_scale_factor_variance to it again.
void write_scaling(Table &reflections, const ScaleData &data,
                   const std::vector<double> &g,
                   const std::vector<double> &g_variance);

} // namespace mxi
