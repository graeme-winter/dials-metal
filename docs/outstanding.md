# What is outstanding

STATUS: the one list of open work, gathered from every document. An item's
evidence is in the document named; when an item is done, it moves to "Closed"
at the end and the document it concerns says what was found. Numbers do not
change when items close. Items are grouped by where the work is, and
marked **correctness** where a number is known to be wrong, or **capability**
where something is missing.

## Integration

* **1.** **correctness -- Partials come out high.** Judged by symmetry equivalents,
   observations with partiality below 0.9 are 10.5 per cent high, and those in
   the first two and last five images 4 to 5 per cent: partiality looks
   underestimated at the ends of the scan, so dividing by it overshoots.
   `docs/scaling.md`, the CC half of the lowest shell.
* **2.** **correctness -- Gap-crossing reflections are about 3 per cent low**, 3.4 by
   symmetry equivalents in scaling; the suspect is profile learning from
   reflections that are themselves cut. `docs/integration.md`, open questions.
* **3.** **correctness -- sigma_m is 0.129 degrees here against DIALS' 0.089** on the
   same 300 image sweep; it sets how large the boxes are. Not yet investigated.
* **4.** **correctness -- The background is about one per cent high at low
   resolution.** A guard ring around the foreground would likely fix it.
   `docs/integration.md`.
* **5.** **correctness -- A dependence of the z centre on strength remains** after
   `--postrefine`, +0.018 images for the weakest to -0.036 for the strongest:
   the integrator's centre on an asymmetric rocking curve.
   `docs/integration.md`.
* **6.** **correctness -- For a thin spot `xyzres.px` in z reports the foreground
   window.** `docs/integration.md`.
* **7.** **capability -- Overlapping reflections are not detected** (Leslie 6.3,
   6.7.1), nor overloads: there is no overload flag. `docs/integration.md`.
* **8.** **decision -- Whether `--postrefine` should be the default**, and what it does
   to scaled intensities: it removes the z offset, at twice the integration
   time. `docs/integration.md`.
* **9.** **correctness -- A reflection predicted just off the panel keeps a
   calculated position of zero** after refinement; `update_predictions`
   reports those rows, and what to write for them is open. `src/refine.hh`.
* **10.** **structure -- Integration is one body inside `run_program`**, so
    `--postrefine` composes two runs of the program. Making it a function --
    one public function per task -- is the refactoring that would remove that.
* **11.** **speed -- Reading frames and profile fitting are a third of the time
    each** on a 16M sweep. `docs/gpu.md`.

## Refinement

* **12.** **correctness -- Outlier rejection is milder than DIALS'**: 11714
    reflections kept against 9580 on the same sweep, RMSD 0.30 px against 0.19.
* **13.** **capability -- No e.s.d.s for the cell**: refinement computes no
    covariance. The scaling covariance shows how.
* **14.** **capability -- The cell is not constrained to the lattice's symmetry**:
    77.93 77.90 77.89 A after `mxi_symmetry`, where DIALS reports 77.898.
    `docs/symmetry.md`.
* **15.** **correctness, untested -- Three conventions with no data to settle them**:
    the order of goniometer composition, `first_image != 1`, and multi-panel
    detectors. `docs/conventions.md`, still open.

## Symmetry

* **16.** **capability -- One sweep only.** With more than one, sweeps must be put on
    the same side of the indexing ambiguity -- dials.symmetry and this already
    choose opposite sides on one cubic sweep -- by a reference or by
    dials.cosym's method. `docs/symmetry.md`.
* **17.** **difference -- Normalisation is by resolution shells**, where DIALS fits an
    anisotropic maximum-likelihood model; the identity CC is 0.858 against
    0.931 on one sweep. The first thing to revisit if a harder case disagrees.
* **18.** **difference -- The resolution limit takes the last shell above the
    threshold**, where DIALS fits a curve: 1.85 A against 2.46 from <I>/<sigma>
    on one sweep.
* **19.** **decision -- Monoclinic candidates are named C 1 2/m 1**, the reference
    setting, where DIALS chooses the setting with beta nearest 90, I 1 2/m 1.
* **20.** **capability -- No obliquity (delta) column** in the table of subgroups.

## Scaling

* **21.** **capability -- One sweep only.**
* **22.** **capability -- No free-set validation**, the check for overfitting.
* **23.** **capability -- No resolution estimate**, and so no "Suggested" column.
* **24.** **capability -- Friedel mates are always merged in scaling**; there is no
    anomalous option.
* **26.** **capability -- The scaling model is not written** to `scaled.expt`, as
    dials.scale writes its `scaling_model`.

## New programs

* **27.** **mxi_report**: HTML reports of merging statistics, as dials.scale's.
* **28.** **mxi_export**: MTZ and mmCIF.

## The spot finder

* **30.** **structure -- What is still duplicated** between the spot finder and the
    rest. `docs/spotfinder.md`.

* **34.** **untested -- Three paths of the CUDA backend**: its 32-bit
  instantiation, which wants 32-bit data, and its `stage0_tile` and
  `stage2_direct` kernels, which
  `SPOTFINDER_GPU_STAGE0=tile SPOTFINDER_GPU_STAGE2=direct` selects on any data.
  The paths it does run agree with Metal exactly. `docs/spots.md`.
* **35.** **speed -- Where the CUDA run's time goes**: 196 images a second on an
  RTX 4060 against Metal's 502. A 16-bit frame of 3108 x 3262 pixels is 20 MB,
  so 196 a second is 4 GB/s across PCIe, which Apple silicon's unified memory
  does not pay; decompression on the host is the other candidate. A run
  without `-gpu` on each machine would say which. `docs/spots.md`.

## Infrastructure

* **31.** **The C++ suite takes about two minutes**, most of it a handful of older
    refinement and indexing tests at 10 to 13 seconds each.
* **32.** **The device port is designed but not written.** `docs/gpu.md`.
* **33.** **Three documents describe the spot finder** -- `docs/spots.md`,
    `docs/spots_notes.md`, `docs/spotfinder.md` -- which could be one.

## Closed

* **25.** dials.merge on `scaled.refl`: it takes it, run by Graeme.
* **29.** The Metal backend of the spot finder: it works. 1800 images of 3108 x
  3262 pixels thresholded in 3.6 seconds, 502 a second, on a Mac with 16
  threads; 146118 spots. `docs/spots.md`.

## Found in DIALS

Not work here, but worth taking upstream, each measured and in the document
named.

* **dials.scale's error model is biased low** by (n-1)/n: its normalised
  deviations are multiplied by sqrt((n-1)/n) (`calc_deltahl`) where the exact
  variance of a difference from a mean the observation is part of is needed. A
  planted a of 1.3 comes back 1.158 in groups of eight and 0.662 in pairs.
  `docs/scaling.md`.
* **Its anomalous slope and dI/s(dI) inherit that bias**: 1.890 and 1.476 on a
  sweep with no measurable anomalous signal, where the same data give 0.955 and
  0.787 with the exact variance. `docs/scaling.md`.
* **dials.scale's scale-factor variance is not a propagation**: the variance is
  multiplied by 1 + sigma_g / g after the error model, rather than
  I^2 var(g) / g^2 added before it. `docs/scaling.md`.
* **Refinement against the spot finder's centres leaves a z offset** that
  depends on strength and on orientation -- about 0.1 images, period 180
  degrees -- because those centres depend on strength. DIALS' integration shows
  it too. `docs/integration.md` and `docs/integration_history.md`.
