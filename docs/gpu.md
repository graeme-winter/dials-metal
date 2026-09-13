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

## What would still be wrong

Analytical derivatives would remove the per-parameter factor entirely -- one
evaluation per reflection instead of 42 -- and are perhaps twenty times faster
again. They were deliberately not written, because a wrong analytical
derivative does not crash: it converges smoothly to the wrong answer and
reports a small residual doing it. If they are written, the numerical version
is the oracle they must be checked against, and it should stay in the tree for
that purpose.

The order in which to do this, then:

1. Run the whole target in `float` on the CPU. Compare residuals and the
   refined model against `double`.
2. Port the target evaluation, keeping numerical differentiation, accumulating
   the normal matrix on the device.
3. Exploit the banding.
4. Only then consider analytical derivatives, with the numerical version as
   the acceptance test.
