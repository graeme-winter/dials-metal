# Conventions, and how to falsify them

Everything in this file is a belief about how DIALS defines its geometry. None
of it has been checked against an `.expt` written by DIALS, because none was
available when it was written.

This is not a theoretical worry. The `.refl` reader in `mxeq` was wrong about
its format in two ways and passed seventy tests, because the only thing
checking it was a round trip against its own writer. The same trap is open
here and it is wider: the closed-loop test predicts reflections and maps them
back, so **any convention error that appears in both the forward and the
reverse map cancels exactly and the tests stay green**. A green suite here
means self-consistency and nothing more.

Each convention below has the test that would falsify it. All of them need one
real `.expt` and its `indexed.refl`.

---

## 1. The sign of s0

**As implemented.** `s0 = direction / wavelength`, with `direction` the unit
vector along beam propagation, source towards sample.

**If wrong.** Every scattering vector `q = s1 - s0` has the wrong origin. The
Ewald sphere sits on the wrong side and prediction produces reflections at
angles that do not occur.

**Falsifying test.** Read a real `.expt`, take its crystal and its scan, and
predict. Compare against the `xyzcal.px` in the matching `indexed.refl`. If s0
is inverted, nothing will match at all — this is a loud failure, not a subtle
one.

## 2. The goniometer composition

**As implemented.** `R(phi) = setting * rotation(axis, phi) * fixed`, taking a
vector from the crystal frame to the laboratory frame. The inverse is applied
as `fixed^T * rotation(axis, -phi) * setting^T`, which assumes both matrices
are orthogonal — true for goniometer matrices, and not checked.

**If wrong.** With `fixed` and `setting` both identity, which is the common
single-axis case, no composition error is visible at all. It only appears on a
multi-axis goniometer or a non-trivial setting rotation.

**Falsifying test.** Needs an `.expt` with a non-identity `fixed_rotation` or
`setting_rotation`. A single-axis insulin dataset cannot test this, and a green
suite on insulin says nothing about it.

## 3. The rotation axis and the sense of phi

**As implemented.** Right-handed rotation about `rotation_axis` by increasing
phi with increasing image number.

**If wrong.** Predicted reflections appear at `-phi`, so a reflection predicted
near the start of the scan is observed near the end. On a symmetric sweep this
can look almost plausible, which makes it the most dangerous of the set.

**Falsifying test.** Compare predicted `z` against observed `z` for indexed
reflections on a real dataset. A sign error shows as a reflection of the
`z_pred` against `z_obs` line about the scan midpoint, not as scatter.

## 4. The detector pixel-to-lab map

**As implemented.**
`lab = origin + px_fast * pixel_size[0] * fast + px_slow * pixel_size[1] * slow`,
with `px_fast` the first component of `xyzobs.px` and `px_slow` the second.

**If wrong.** The likeliest error is fast and slow being swapped, which on a
square detector is a transpose that survives a beam-centre check and shows up
only as a systematically wrong cell.

**Falsifying test.** Take real `xyzobs.px.value` from an indexed `.refl`, map
to reciprocal space, and require the result to be close to `A h` with `A` and
`h` from the same files. Millipixel agreement confirms the whole chain. This is
the single most valuable test to run, because it checks 1, 3, 4 and 6 at once.

## 5. Whether `entering` means entering

**As implemented.** `entering` is true where `s0' . (m2 x r) < 0`, that is,
while the reciprocal lattice point is approaching the Ewald sphere.

**If wrong.** Every flag is inverted. Nothing about the geometry breaks, and
prediction still lands in the right place. What breaks is every keyed join
against DIALS output, because `(hkl, entering, id)` then pairs the entering
observation of a reflection with its exiting one — the two are at different
angles on the scan, so the comparison reports a large disagreement that has
nothing to do with intensity.

**Falsifying test.** Join a real `indexed.refl` on `(hkl, entering)` and check
that `xyzcal.px` z agrees. If the flag is inverted the z values pair up
crossed, which is unmistakable.

## 6. The scan z anchor

**As implemented.** `phi = osc_start + (z - z_offset) * osc_width`, with
`z_offset` defaulting to zero, which reproduces DIALS.

**Known wrinkle.** DIALS anchors z to image number one rather than to the start
of the array. The two agree when `first_image == 1`, which is nearly always,
and differ by a constant otherwise. `z_offset` exists so that the choice is a
visible parameter rather than a baked-in constant, and so that the day someone
processes a scan starting at image 0 there is one number to change.

**Falsifying test.** A dataset with `first_image != 1`. Insulin will not do it.

## 7. The setting matrix

**As implemented.** `A` has a\*, b\*, c\* as its **columns**, so `r0 = A h`.
The real-space vectors are the **rows** of `A^-1`, so
`A = from_rows(a, b, c)^-1` — the plain inverse, not the inverse transpose.

**Status.** This one is already tested properly, in
`setting_matrix_is_the_inverse_not_the_inverse_transpose`, using a triclinic
cell. On a cubic cell both forms are diagonal and the wrong one passes, which
is why the test cell is deliberately awkward. It is listed here because it is
the same class of error and because reading a real `.expt` will exercise it
against DIALS' own serialisation of `real_space_a/b/c`.

---

## What to do with one real `.expt`

In rough order of how much each buys:

1. **`indexed.expt` + `indexed.refl` from the same run.** Map `xyzobs.px.value`
   to reciprocal space, compare against `A * miller_index`. Agreement to a
   millipixel validates conventions 1, 3, 4, 6 and 7 in one go, and the residual
   distribution says which one is wrong if it fails.
2. **The same pair, checking `entering`.** Validates 5.
3. **Predict from the real `.expt` and compare to real `xyzcal.px`.** The end to
   end check, and the one that would be the acceptance test for prediction.
4. **Anything with a non-identity setting or fixed rotation.** Validates 2,
   which nothing else can.
