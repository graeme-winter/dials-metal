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

## Precision: measured, and less of an obstacle than expected

Apple GPUs have no double precision at all, and the geometry here is carried in
double throughout. The worry is finite differences: a derivative is the
difference of two nearly equal residuals, and catastrophic cancellation in
float32 could leave nothing.

Measured, on a step of 1e-6 relative -- what the refinement currently uses --
the residual changes by a median of 1.5e-3 of its own size. Against float32's
epsilon of 1.2e-7 that leaves **4.1 significant digits** in the derivative. With
a step of 3e-4 relative, which is near sqrt(epsilon) and is what a float32
implementation should use, the change is 0.46 relative and **6.6 digits**
remain.

So numerical differentiation survives float32, which was not obvious. A
Gauss-Newton step only has to be a descent direction, and it is accepted only if
the residual actually falls, so four digits is ample.

The residual itself is a weaker constraint than it looks. It is a difference of
detector positions of order 2000 px, so float32 rounding contributes about
2.4e-4 px, against residuals of 0.3 px: one part in a thousand, random per
reflection, averaging away over thirteen thousand of them.

**The experiment to run first needs no GPU.** Compute the whole target in
`float` on the CPU and compare against the `double` version, reflection by
reflection, and then run a full refinement that way and compare the refined
model. That is cheap, decisive, and tests the arithmetic rather than an
estimate of it. If it holds, the port is mechanical; if it does not, the
measurement says exactly where.

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

1. Run the whole target in `float` on the CPU and compare against `double`.
2. Port, accumulating the normal matrix on the device, exploiting the banding.

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
