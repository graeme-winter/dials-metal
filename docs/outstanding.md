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
* **18.** **difference -- The <I>/<sigma> limit takes the last shell above
  the threshold**, where dials.symmetry fits a curve: 1.85 A against 2.46 on
  one sweep. The CC half limit is now dials.estimate_resolution's tanh fit.
* **19.** **decision -- Monoclinic candidates are named C 1 2/m 1**, the reference
    setting, where DIALS chooses the setting with beta nearest 90, I 1 2/m 1.
* **20.** **capability -- No obliquity (delta) column** in the table of subgroups.

## Scaling

* **21.** **capability -- One sweep only.**
* **22.** **capability -- No free-set validation**, the check for overfitting.
* **24.** **capability -- Friedel mates are always merged in scaling**; there is no
    anomalous option.
* **26.** **capability -- The scaling model is not written** to `scaled.expt`, as
    dials.scale writes its `scaling_model`.

* **41.** **measure -- `mxi_scale --threads`**: byte-identical on any count,
  speed unmeasured -- one core here. Serial on the 16M sweep: 8.0 s.
  `docs/scaling.md`.

## New programs

* **27.** **mxi_report**: HTML reports of merging statistics, as dials.scale's.
* **28.** **mxi_export**: MTZ and mmCIF.

## The spot finder

* **30.** **structure -- What is still duplicated** between the spot finder and the
    rest. `docs/spotfinder.md`.

* **34.** **untested -- Two paths of the CUDA backend**: its 32-bit
  instantiation, which wants 32-bit data, and its `stage2_direct` kernel
  (`SPOTFINDER_GPU_STAGE2=direct`). `stage0_tile` has run on real data on an RTX
  4060, output byte-identical to the direct kernel's, and slower. The paths it
  runs agree with Metal exactly. `docs/spots.md`.
* **35.** **speed -- Where the spot finder's time goes.** `mxi_find --timing`,
  on 3600 frames of 4148 x 4362 pixels, 10 GB compressed, with 16 threads:
  * **Metal, 13.9 s:** reading from HDF5 39 per cent of the threads' time, 24 ms a
    frame from an external drive -- 10 GB in 13.9 s is 0.72 GB/s, in the range of
    a USB 3.2 drive; a second run with the file in the page cache would say
    whether it is the drive or HDF5's lock -- thresholding 32, decompressing 16.
  * **CUDA, RTX 4060, 32.0 s: the kernels are the limit.** 8.9 ms a frame
    against Metal's 3.9. With 8 threads the wall is the same, 32.5 s, and only
    the waiting falls -- reading back 138 thread-seconds to 43, the rest 137 to
    45 -- while the device's work stays about 110 s: the GPU runs some 3.4
    frames' kernels at once and that is its rate. Of that work, at 8 threads:
    stage 2 (11 x 11) 38 per cent, stage 0 (7 x 7) 27, stage 1 (5 x 5) 23, the
    upload 11 (3.2 ms for 36 MB, 11 GB/s). Stage 0 by tile is slower on this
    card, 37.4 thread-seconds against 32.4, so direct stays the default there,
    the reverse of Apple silicon. Faster means the windowed sums cost less:
    running sums independent of the window's size, or fused stages, measured on
    the card. A fused kernel for that is written (`SPOTFINDER_GPU_FUSED=1`,
    `docs/spots.md`), verified against the CPU tile by tile and compiled for the
    4060, not yet run: `ctest -R dext_gpu`, then `cmp` and `--timing` against the
    three kernels on real data. HDF5's lock shows too: reading takes 4.6 thread-seconds on 8
    threads and 22.5 on 16, for the same data. `docs/spots.md`.
* **36.** **speed -- Prediction took 29 per cent of one integration and 1 of
  another**: 6.3 s of a 300 image integration of a scan-varying model, nearly as
  long for 60 images, and 0.2 s of a 3600 image one on the 16M sweep. What makes
  the first slow is not yet known. `mxi_integrate --timing`.
* **37.** **speed -- `mxi_symmetry` names every subgroup**, 36 per cent of its
  run, for the table; the choice needs only the one chosen.

* **42.** **measure -- One pass over the images on the 16M sweep.** On a smaller
  data set on a MacBook it is byte-identical to two passes (same SHA-1) and 7.0 s
  against 8.1, the second pass's decompression gone. On the 16M sweep, where two
  passes decompressed for 80 thread-seconds of a 28.4 s integration, the saving
  and the most shoeboxes held -- estimated at 49000 and 1.1 GB -- are still to be
  measured: `mxi_integrate --timing` against `--two-pass`. `docs/integration.md`.
* **40.** **speed -- `mxi_symmetry` scores its elements by ordered-map lookup**,
  1.3 of its 4.2 s on the 16M sweep; a hash table, or sorted indices, would do.

## Infrastructure

* **31.** **The C++ suite takes about two minutes**, most of it a handful of older
    refinement and indexing tests at 10 to 13 seconds each.
* **32.** **The device port is designed but not written.** `docs/gpu.md`.
* **33.** **Three documents describe the spot finder** -- `docs/spots.md`,
    `docs/spots_notes.md`, `docs/spotfinder.md` -- which could be one.

## Closed

* **25.** dials.merge on `scaled.refl`: it takes it, run by Graeme.
* **38.** Integration decompressing every frame twice: one pass over the
  images now, each reflection fitted once its scan blocks are final, the table
  byte-identical to two passes'. A sparse profile cell borrows its own scan
  block's average rather than the whole scan's, which one pass needs, and a
  scan block is 10 degrees by default. `docs/integration.md`.
* **39.** Writing a reflection table at about a second a million rows: a
  column is copied, not appended a byte at a time, and written straight to the
  file rather than built into one string first. 630930 rows, 261 MB: 0.875 s to
  0.14-0.19, of which the disk is 0.06; reading 0.29 to 0.24. Every table
  byte-identical. `src/refl.cc`.
* **23.** A resolution estimate: dials.estimate_resolution's tanh fit through
  CC half, at 0.3, and its significance limit, with the "Suggested" column.
  `docs/scaling.md`.
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
