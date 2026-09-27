# Symmetry and scaling: a plan

STATUS: a plan, being built. Decided: gemmi as the dependency for space groups
(a submodule at v0.7.5, wrapped in `src/symmetry.hh`), and mxi_scale first.
Built so far: space groups and changes of basis. Program names here are in plain
text, not code, while they do not exist yet.

The pipeline ends at integration. What follows it in DIALS is dials.symmetry,
choosing the point group and space group, and dials.scale, putting every
observation on one scale and refining an error model. This plans mxi_scale and
mxi_symmetry, the scaling after Beilsten-Edmands et al. (2020), Acta Cryst.
D76, 385-399, parameterised differently where noted.

## Reference numbers

The DIALS logs for the 300 image insulin sweep that the integration was compared
on also contain dials.symmetry and dials.scale. They are what the two programs
are checked against first.

* **Symmetry.** The best solution is I m -3, from a P 1 cell of 67.42 67.46
  67.46 A, 109.44 109.47 109.46 degrees, to 77.898 A cubic.
* **Scaling.** 11 parameters -- 6 for scale, 5 for decay, no absorption surface
  for a 30 degree sweep -- and DIALS warns that 64 per cent of them are
  poorly determined. Error model a = 0.531, b = 0.0246: a below one, so the
  integration's sigmas were overestimated for this sweep. Overall
  completeness 71 per cent, multiplicity 2.7, I/sigma 31.4, Rmeas 0.045, CC half
  0.999, resolution suggested 1.63 A at CC half = 0.3.

## mxi_scale

### The model

Inverse scale g = C(phi) . T(t, d) . S(s0, s1), as the paper's physical model.

* **Scale C(phi).** A cubic B-spline in rotation angle, where DIALS uses a
  Gaussian-weighted average of the nearest three parameters. Control points 10
  to 20 degrees apart, as the paper finds, and fewer on a narrow sweep: DIALS
  used six on 30 degrees, and most were poorly determined. The clamped uniform
  B-spline the crystal model uses (`src/geometry.cc`, `Crystal::A_at`) serves,
  in scalar form.
* **Decay T = exp(B(t) / 2d^2).** B(t) a cubic B-spline in t, image number for
  now as a stand-in for dose, with DIALS' weak restraint of sum B_i^2 toward
  zero. The zero-dose work belongs here later: t as dose, and the extrapolation
  as a use of the fitted B(t).
* **Absorption S = 1 + sum P_lm [Y_lm(s1) + Y_lm(s0)] / 2** in the crystal frame,
  as the paper's equation 7: l_max = 4 by default, 24 parameters, odd terms
  included since they absorb miscentring, with a restraint of sum P_lm^2 toward
  zero. Off for narrow sweeps, as DIALS turns it off, since there is not the
  angular coverage to determine it.

One overall scale is degenerate with the merged intensities and must be fixed:
the mean of C, or one control point.

### Fitting

The target is the paper's equation 2, sum w (I - g <I>)^2 plus restraints,
with <I_h> in closed form (equation 3) for the current g. At about 70
parameters the normal matrix is tiny, so full-matrix Levenberg-Marquardt with
analytic derivatives can run throughout, where DIALS uses L-BFGS until the last
cycle. The cost is building J^T J, linear in observations; the final inversion
gives parameter uncertainties, propagated to the scale factors' variances.

Around the fit, as the paper's figure 2:

1. **Reflection selection.** A random subset of symmetry-equivalent groups, at
   least 2000 groups and 50000 reflections (section 4.4), for optimisation; the
   model is applied to everything.
2. **Outlier rejection.** A reflection's normalised deviation from the weighted
   mean of its group, excluding itself, above 6 sigma (Evans 2006). Repeated as
   the model improves, retesting everything.
3. **Intensity combination.** I = w I_prf + (1 - w) I_sum with w = 1 / (1 +
   (I_sum / I_mid)^3), I_mid chosen from logarithmically spaced values by Rmeas.
4. **Error model.** sigma'^2 = a^2 [sigma^2 + (b I)^2]. a from the slope of the
   central normal probability plot (|x| < 1.5) of the normalised deviations; b
   by minimising the paper's equation 17 over logarithmically spaced intensity
   bins with a fixed; alternated to convergence from a = 1, b = 0.02, on groups
   with <I> above 25 and <I / sigma^2> above 0.85.
5. **Free-set validation.** Rmeas on a tenth of the groups held out, as the
   check for overfitting that DIALS' own warning above is a symptom of.

### Output

A scaled, unmerged reflection table -- inverse_scale_factor and its variance,
the adjusted intensity variances -- that dials.merge and dials.export can take
on to MTZ. And, in the manner of the other programs' reports, merging
statistics by resolution shell: completeness, multiplicity, I/sigma, Rmeas,
Rpim, CC half.

### Checking it

* Against dials.scale on the same integrated input: scale factors image by
  image, the error model's a and b, merging statistics by shell.
* A plant: synthetic intensities with known C, B and S applied, recovered.
* Free-set Rmeas, work against free.
* Whether `mxi_integrate --postrefine` changes the scaled result, which was left
  open when it was measured.

## mxi_symmetry

As dials.symmetry and POINTLESS do it:

1. **The lattice's metric symmetry**, from the reduced cell, by a search for
   two-fold axes as Le Page's. Every candidate point group is a subgroup of it.
2. **Normalised intensities**, E^2 in resolution shells, so that correlations
   compare like with like.
3. **A score for each symmetry element** from the correlation between
   reflections it relates, and **Laue group likelihoods** combining them
   (Evans 2011).
4. **Systematic absences** along screw axes, for the space group.
5. **Reindexing** to the chosen setting.

The multi-crystal case -- indexing ambiguities between crystals, as dials.cosym
resolves them -- comes after, and has its own prior work to build on.

## Decisions

* **Scaling first.** mxi_scale takes the space group from the model, or from
  the command line with a change of basis until mxi_symmetry can choose one:
  for the 300 image sweep, I 2 3 and `b+c,a+c,a+b` as dials.symmetry reported.
  That operator takes this pipeline's own indexing to 77.93 77.90 77.89 A at
  89.99 89.94 89.93 degrees, and leaves none of 20914 reflections forbidden by
  the I centring.
* **gemmi**, pinned at a release, for operators, the reciprocal asymmetric unit,
  absences and centric reflections, behind `src/symmetry.hh`.
* **Reports** of merging statistics, in the manner of dials.scale's HTML, go in
  a later mxi_report. mxi_scale prints its summary tables, as the others do.
* **MTZ and mmCIF** go in a later mxi_export. mxi_scale writes a scaled,
  unmerged reflection table and the models.
