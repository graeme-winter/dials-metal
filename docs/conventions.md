# Conventions: what was checked, and what it cost

These were beliefs. They have now met a real `imported.expt`, `indexed.expt`
and `indexed.refl` from `dials.index` on a 300-image insulin sweep.

**Two were wrong.** Neither could have been caught by the test suite that
existed, because the closed-loop test predicts reflections and maps them back,
applying each convention once in each direction so that an error cancels
itself. That is the same trap that let the `.refl` reader in `mxeq` pass
seventy tests against a format wrong in two ways.

`tests/test_real_geometry.cc` now checks all of this against embedded real
numbers on every build.

---

## Wrong: the sign of s0

**Was.** `s0 = direction / wavelength`.

**Is.** `s0 = -direction / wavelength`. dxtbx stores the beam direction
pointing from the sample back towards the source, so s0, which points along
propagation, is its negation. The detector origin in the real file has
`z = -167` with a beam direction of `+z`, which is the visible tell.

**What it cost.** Residual against `A h` of 2.097 per reciprocal Angstrom,
against a reciprocal cell edge of 0.0148 — a hundred and forty lattice
spacings. This fails loudly, which is the one mercy.

**Knock-on.** The synthetic fixture in `test_geometry.cc` had the panel at
`+z`, which was consistent only with the wrong sign. Correcting the sign made
it predict nothing at all — and a geometry fixture that cannot diffract is a
silent test, not a failing one. The prediction tests now assert counts so that
this cannot recur.

## Missing entirely: the parallax correction

**Was.** `mm = px * pixel_size`.

**Is.** Winter *et al.* (2018), *Acta Cryst.* **D74**, 85-97, Appendix A. An
X-ray entering a 0.45 mm silicon sensor at an angle travels into it before
being absorbed, so the pixel that fires is displaced outward from where the ray
met the front face. The implementation here is eqn (7) for the attenuation
length and eqn (8) for the displacement, arrived at by fitting against real
files before the paper was to hand and identical to it. With `mu = 3.663 /mm`:

```
depth  = 1/mu - (t/cos_t + 1/mu) * exp(-mu * t / cos_t)
offset = depth * (u . fast), depth * (u . slow)
px -> mm:  mm = px * pixel_size - offset(px * pixel_size)
mm -> px:  px = (mm + offset(mm)) / pixel_size
```

`u` is the unit vector from the sample to the point, `cos_t` its angle to the
panel normal. Each direction evaluates the offset at *its own input*, which is
not exactly self-inverse — and is exactly what DIALS does. Both directions
reproduce DIALS to the last bit of a double: `px_to_mm` against
`xyzobs.mm.value`, `mm_to_px` against `xyzcal.px`.

**What it cost.** 0.048 mm median, 0.118 mm maximum — **1.58 pixels**, twenty
times the centroid precision. Omitting it does not shift the model by a
constant but by a smooth function of scattering angle, which refinement partly
absorbs into the detector distance and leaves as a radial residual.

**Which way round.** Eqn (8) describes a *predicted ray* arriving at the
detector, so adding the offset is millimetres to pixels and subtracting it is
pixels to millimetres. Fitting alone could not have settled that; it only said
the sign was consistent.

**`xyzobs.mm` is not a usable quantity.** It is computed once at import and
never recomputed, so it carries the inverse parallax correction under whatever
geometry was current then, and after refinement it describes a detector that no
longer exists. Reproducing it requires the **imported** panel; `xyzcal.px`
requires the **refined** one, and mixing them leaves 1.7e-3 mm that looks like
a flaw in the parallax model and is not.

Nothing in `src/` reads it, and nothing should: work from `xyzobs.px` through
the current detector model every time. The same conclusion has been applied to
`mxeq`, where it was a fallback join key -- two pipelines that refined to
slightly different detectors would have been compared in two different
millimetre frames.

The third component is the exception and stays usable: it is the rotation
angle, which comes from the scan and has no dependence on the detector model,
so refinement cannot make it stale.

## Wrong shape: the goniometer

**Was.** `rotation_axis`, `fixed_rotation`, `setting_rotation`.

**Is.** A multi-axis goniometer: `axes`, `angles`, `scan_axis`. Axes below the
scan axis are carried by the sample and compose into the fixed rotation; axes
above it move the scan axis itself and compose into the setting rotation.

The entry in `angles` for the scan axis itself is ignored, and must be: that
axis does not have one setting, it has a different one on every frame, and the
scan supplies it.

**Still untested.** Every other angle in this dataset is zero, so both compose
to the identity. **A single-axis insulin sweep cannot test the decomposition at
all.** It needs a dataset with a non-zero chi or phi setting -- the chemical
crystallography case in the DIALS paper, four sweeps on a fixed-chi goniometer
at 57.74 degrees, is exactly the shape of data that would.

## Wrong shape: the scan

**Was.** `oscillation: [start, width]`.

**Is.** `properties.oscillation` is a per-image array of start angles.

Treating the width as constant is fair, and now checked rather than assumed:
`Scan::from_oscillation` records the worst departure from a constant width as a
fraction of the width. On this data it is 3.6e-14, which is accumulation
round-off and nothing else.

The width is taken from the **endpoints**, not from the first two elements. The
array is generated by repeated addition, so adjacent differences subtract two
nearby doubles and lose most of their significance; spanning the scan divides
that error by the number of images. One bad value near the start would
otherwise set the width for the whole sweep.

---

## Confirmed correct

**The scan angle.** `phi = radians(osc_start + z * osc_width)`, exact to 1e-12
against `xyzobs.mm.value[2]`, which is the rotation angle in radians.

**`entering`.** `s0 . (m2 x r) < 0`. Agrees with DIALS on 100% of 13072 indexed
reflections. Worth noting how quietly this one would have failed: with the flag
inverted, prediction still lands in exactly the right place and only the keyed
joins in `mxeq` break, pairing each reflection's entering observation with its
exiting one.

**The setting matrix.** `A` has a\*, b\*, c\* as columns; real-space vectors
are the rows of `A^-1`, so `A` is the plain inverse. Confirmed against DIALS'
`real_space_a/b/c` and its own `rlp` column.

**The whole chain.** Real `xyzobs.px.value` mapped into the crystal frame lands
on `A h` with median residual 1.9e-4, and agrees with DIALS' own `rlp` column
to 4.3e-5 — forty times closer than either sits to the ideal lattice. What is
left is the indexing residual, not an error here.

---

## Still open

1. **The goniometer decomposition.** Needs a non-zero setting angle.
2. **`first_image != 1`.** The z anchor wrinkle is still untested; insulin
   starts at image 1, where the two conventions agree by construction.
3. **Multi-panel detectors.** One panel here, and the `hierarchy` block is
   ignored entirely.
4. **The provenance of the `s1` column.** Normalised to exactly
   `1/wavelength`, and 9.5e-5 from the value computed from `xyzobs.mm` under
   the refined geometry -- closer than any other candidate tried, but not
   equal. Given that `xyzobs.mm` is itself import-time, a column derived from
   it under a refined detector is mixing two geometries, which would explain a
   residual of about this size. Nothing here depends on `s1`; like `xyzobs.mm`
   it is best not consumed.
