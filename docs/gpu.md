# Computing the refinement target on a device

An exploration, with the measurements that motivate it. Nothing here is
implemented yet.

## It is the only thing worth moving

Measured on insulin, 13072 indexed reflections, single-threaded:

    one target evaluation        4.0 ms        0.31 us per reflection
    static refinement            0.62 s        15 parameters, 10 steps
    scan-varying, 9 points       5.22 s        87 parameters, 13 steps

A numerical Jacobian needs one evaluation per parameter plus one for the base,
so the scan-varying refinement does 88 x 13 = 1144 evaluations, about 15
million reflection-evaluations. Multiplying that out gives 4.58 s of the 5.22 s
measured, **88 per cent**; for the static case the same estimate comes to 106
per cent, which is to say everything and a little measurement noise.

There is no point accelerating anything else. The normal equations are 87 x 87,
the outlier rejection is a median, and the I/O happens once.

## The shape of the work

One reflection-evaluation is: build the reciprocal lattice point from the
setting matrix, solve the Ewald condition for the rotation angle, rotate,
intersect the panel plane, apply the parallax correction, subtract the
observation. Around a hundred flops and five to ten transcendentals -- `atan2`,
`acos`, `exp`, a `sin`/`cos` pair.

Every one of them is independent. The natural decomposition is two-dimensional,
one thread per (reflection, parameter): 13072 x 88 = 1.15 million for insulin,
and 11000 x 200 = 2.2 million for a four-sweep scan-varying l-cysteine. That is
enough parallelism to saturate anything.

**Do not materialise the Jacobian.** At 13072 x 3 x 87 doubles it is 27 MB per
iteration, and the only thing it is used for is `J^T W J` and `J^T W r`. Each
thread should accumulate its reflection's contribution into the normal matrix
directly, leaving a reduction over reflections of an 87 x 87 symmetric matrix --
3828 upper-triangle entries, small enough to hold in shared memory per
threadgroup and combine at the end.

## The B-spline makes the Jacobian banded

This is the part that connects to the choice of interpolation. A cubic
B-spline has local support: the setting matrix at scan position t depends on
exactly four control points. So a reflection's residual depends on 4 x 9 = 36
crystal parameters, not on all 9N of them, plus the 6 detector parameters.

At 9 control points that is 42 of 87 parameters per reflection, a saving of
half. At 30 control points it would be 42 of 276, a saving of six sevenths. An
interpolating spline would not have this property -- its global solve couples
every control point to every observation -- which is a concrete reason to
prefer the B-spline beyond its smoothness.

The banding is by scan position, so reflections sort naturally into bands by
image number, and a threadgroup covering one band touches a contiguous slice of
the parameter vector.

## Precision: measured, and the earlier estimate was wrong

Apple GPUs have no double precision at all. `src/target.h` is the whole target
written once and templated on the scalar type, so the same code compiles at
both precisions -- which matters, because a separate float implementation
disagreeing with the double one could be the precision or could be a
transcription error, and the measurement could not tell them apart. Compiled
with `double` it reproduces `centroid_residual` to 5e-13 px; only then does the
`float` result mean anything.

**The residual survives float32.** Median difference 2.2e-4 px against
residuals of 0.32 px -- under a tenth of a per cent, random per reflection, and
averaging away over thirteen thousand of them.

**A finite-difference derivative does not.** Compared against the analytical
derivative, on real insulin:

    double, step 1e-6 relative          median error 4.4e-07     6.4 digits
    float,  step 1e-6 relative          median error 1.00        none at all
    float,  step 3e-4 (near sqrt eps)   median error 1.4e-02     1.8 digits

At the step the refinement uses, a float finite difference is entirely noise.
At the best step available to float it does not reach two digits, and its
ninety-ninth percentile is above five, meaning some entries have the wrong
sign.

### The estimate that was wrong, and why

An earlier entry in this file claimed float32 would leave 4.1 digits in a
numerical derivative. The reasoning was: a 1e-6 relative parameter step changes
the residual by 1.5e-3 of its own size, float epsilon is 1.2e-7, so four digits
survive.

That compares the change against the *residual*. But the residual is a
difference of detector positions of order two thousand pixels, so its absolute
error in float32 is epsilon times the position, about 2.4e-4 px -- not epsilon
times the residual. The change being measured is 1.5e-3 x 0.32 px, about
5e-4 px. Signal and noise are the same size, which is exactly the 1.00 relative
error measured.

The general form of the mistake: **relative precision belongs to the quantity
the arithmetic is carried in, not to the quantity you are interested in.**

### What follows

Analytical derivatives are not a nicety for a device port, they are a
precondition. That reverses the earlier plan, which had them fourth on the list
as an optimisation.

`tests/test_precision.cc` asserts the failure as well as the success, so that
nobody later assumes numerical differentiation would port as it stands.

### The analytical derivative in float32: measured

`src/derivatives_t.h` is the same treatment applied to the derivatives.
Compiled with `double` it reproduces `crystal_derivatives`,
`detector_derivatives` and `beam_derivatives` **bit for bit** -- the same
operations in the same order -- so the float number is not confounded by a
transcription difference.

    finite difference, float, step 1e-6    median relative error  1.00
    finite difference, float, step 3e-4    median relative error  1.4e-02
    analytical,        float               median relative error  1.3e-07

Seven orders of magnitude, and structural rather than lucky: an analytical
derivative never forms the difference of two nearly equal positions, so there
is no cancellation to spend the significance on. The ninety-ninth percentile is
1.2e-5 and the worst case 3.5e-2, the latter on the same near-tangential
reflections that trouble everything else.

So the whole target and its Jacobian can be computed in single precision. The
device port is not blocked on precision, provided the derivatives are
analytical -- which is the conclusion the earlier estimate had exactly
backwards.

## Analytical derivatives: written, and validated

`src/derivatives.h` implements Appendix A of Waterman et al. (2016) for the
crystal parameters. They matter more for a device than for a CPU, and for a
reason that is not speed: a finite difference is a difference of two nearly
equal residuals, and on a float32 device some of the significance is spent on
the cancellation however carefully the step is chosen. Measured earlier, float32
leaves about four digits in a numerical derivative; an analytical one leaves
seven. They also remove the per-parameter factor entirely -- one evaluation per
reflection instead of forty-two.

The chain, with this code's parameterisation:

    dphi/dp   = -(R_phi dr0/dp . s1) / ((e x r_phi) . s0)      eqn (40)
    dr_phi/dp = (e x r_phi) dphi/dp + R_phi dr0/dp             eqn (46)
    dv/dp     = D dr_phi/dp                                    eqn (45)
    dX/dp     = (w du/dp - u dw/dp) / w^2                      eqn (43)

Refining the nine elements of A directly pays for itself here: since r0 = A h,
the derivative with respect to element A(i, j) is just h_j sitting in row i and
zero elsewhere. Through U and B it would be a chain through the metrical
matrix. And the derivative with respect to a B-spline control point is that
same vector times the control point's weight -- so `spline_weights` returns the
four indices and weights, and the banding falls out with no extra work.

**They are checked against central finite differences on the real refined
insulin geometry**, element by element, in `tests/test_derivatives.cc`: median
relative agreement below 1e-8 and the 99th percentile below 1e-5. That test is
the reason they are allowed to exist. A wrong analytical derivative does not
crash; it converges smoothly to the wrong answer and reports a small residual
doing it, so the numerical version stays in the tree as the oracle.

Detector and beam follow Appendix B. Neither the detector nor the beam appears
in r0, and the detector does not appear in s0 either, so `dphi` for a detector
parameter is exactly zero rather than merely small -- asserted, because a
nonzero value there would mean the chain rule had picked up a term that does
not exist.

One piece is not in the paper. DIALS measures its residual in millimetres and
radians; this code measures it in pixels and images, and the parallax
correction sits between the two, so the conversion is a 2x2 Jacobian rather
than a division by the pixel size. It is worth 7e-4 relative on the diagonal at
the corner of an Eiger2 panel, which is small and is a hundred times the
tolerance the derivatives are held to. `Panel::mm_to_px_jacobian` computes it
analytically, for the same reason as everything else here.

**Available as `mxi_refine --analytic`.** On insulin it reaches the same model
to five decimal places in detector distance and four in cell, six times faster;
on four sweeps of l-cysteine with about two hundred parameters, 12 s against
33 s.

### Where the two Jacobians disagree, and why it is not the analytical one

Compared entry by entry over the whole Jacobian -- crystal, detector and beam,
static and scan-varying:

    median               1.5e-9 to 3.6e-8
    99th percentile      1.4e-2 to 1.5e-1
    grossly different    0.07 to 0.3 per cent of entries

The tail is not a gradual loss of accuracy. It is a small set of entries where
the two disagree completely, often in sign. Those are reflections where the
perturbation moves which Ewald root lies nearest the observation, so the finite
difference compares two different branches of the prediction and its value is
meaningless. **The analytical derivative is the correct one there.** That is
worth knowing beyond precision: it means the numerical path has been feeding
refinement a fraction of a per cent of invalid derivatives all along.

Still to do:

1. Port, accumulating the normal matrix on the device, exploiting the banding.
   `target.h` and `derivatives_t.h` are written to be compiled as they stand;
   what is missing is the dispatch, the reduction and the host side.
2. Check the normal matrix itself in float32. Every measurement so far is of a
   single reflection's contribution; accumulating thirteen thousand of them is
   a sum of positive quantities of widely differing size, which is a different
   question and deserves its own measurement rather than an assumption.

## The volume cutoff, and a guess it did not support

Eqn (40) divides by the volume of the parallelepiped formed by the rotation
axis, the reciprocal lattice vector and the beam, which vanishes for
reflections near the rotation axis. DIALS discards any below 0.05. Those are
the reflections with large Lorentz factors and genuinely ill determined phi,
and the analytical derivative makes the reason explicit rather than empirical:
the measured rotation-angle derivative per unit |h| is 3.9 times larger in the
smallest-volume quartile than the largest.

It is **not**, however, the four per cent of reflections whose forward and
reverse maps disagree under a scan-varying model, which is what this code had
previously guessed. Tested: at a cutoff of 0.05 the volume criterion removes
5.4 per cent of reflections and only 12 per cent of the disagreements, leaving
the rate essentially unchanged at 3.97 per cent. Nor is it the iteration count
in the forward map -- three, six and twelve passes all give 4.27 per cent -- nor
reflections whose two Ewald roots are close, which show the same 4.2 per cent as
those whose roots are ninety degrees apart. **That population is unexplained.**
Outlier rejection removes it and refinement then works, which is a workaround
rather than an answer.
