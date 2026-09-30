# Where the time goes, and what a device would take

STATUS: measurements of where integration's, indexing's and refinement's time
goes, and what moving each to a GPU would and would not buy -- a device port of
the refinement target is designed here and not written. The spot finder's GPU
threshold, which is written and runs on Metal and CUDA, is in `docs/spots.md`.
Within each part the sections are as they were written: a later one can
overturn an earlier one, and says so.

## Integration

Profile fitting and reading frames, and whether a device is the way to make
them faster: measured each time before anything moved.

### What a device would buy, measured first

Profile fitting was 53 per cent of integrating a 3600 image Eiger 16M sweep --
the obvious thing to move. It was not the arithmetic that cost: 97 per cent of
it was carrying the reference profile onto the pixels, and nearly all of that
was one call to `epsilon_of` for every subdivision of every pixel, 8100 a box.
Computing it at the pixel corners and interpolating, with the subdivisions near
a cell boundary done exactly, was 4.8 times faster on the CPU and gave the same
answer. Fitting went from 35.8 seconds to 12.2.

Reading frames had the same shape: 5.62 reads a frame where two passes need 2,
from windows that discarded their boxes. Reading each frame once a pass took it
from 19.6 seconds to 11.0.

The run went from 67.7 seconds to 36.9 with no device. Where it stands now:

    reading frames (wall)    11.0 s   30 per cent   decompression 98 thread-s
    profile fitting          12.2 s   33 per cent   independent per reflection
    opening shoeboxes         4.3 s   12 per cent

Profile learning, which uses the same geometry, has since had the same change:
4.3 times faster a box, about two seconds less on that run.

Profile fitting is still the natural device workload -- one independent fit per
reflection, a few thousand voxels each, and a reference profile shared by
thousands -- but it is now a third of a much shorter run, and the case for
moving it should be made from a profile of this version, not the old one.
Reading is decompression on the CPU and would not move with it.

#### Measured again, with one pass over the images

On 1800 frames on a MacBook, 16 threads, one pass: 6.6 s, of which profile
fitting is 2.4 (36.5 per cent), reading frames 1.1, writing 0.7, the profile
model, background and summation, and opening shoeboxes about 0.45 each. Fitting,
split by `mxi_integrate --timing` on the 300 image sweep: carrying the reference
profile onto the pixels 86 per cent of its thread-seconds, the least squares 14,
interpolating the reference 1. And of the carrying, a throwaway probe put 90 per
cent in `pixel_cells`: for each subdivision of each pixel -- 5 x 5 a pixel, some
10000 a box -- where it falls in the grid, from e1 and e2 interpolated between
the pixel's corners and computed exactly within a thousandth of a cell boundary.
Geometry alone: no pixel value enters it.

That is a device's kind of work -- independent across boxes and within them, its
inputs a box's extent, s1 and phi -- with two obstacles. Metal has no double
precision and the RTX 4060 runs it at a sixty-fourth of single, while the
geometry and its boundary test are written in double; in single a subdivision
near a boundary can land in the other cell, so a device's fit would differ from
the CPU's by that, and the byte-for-byte comparisons this code has been held to
would become tolerances. And fitting is a third of the run, so a device can take
at most that. Before a device: `pixel_cells` does twice some work it could do
once -- (e + span) / step and its floor, for the boundary test and again for the
cell -- recomputes the four bilinear weights for every pixel where they depend
only on the subdivision, and grows each pixel's cell list with vector pushes;
and a strong reflection's cells are computed twice, for learning and for
fitting. Those can be removed with the arithmetic unchanged, on the CPU, for
every platform.

Done, with the arithmetic unchanged and integrated.refl byte-identical: the
weights formed once a box, (e + span) / step once, the cells gathered in a
local array, `std::floor` replaced by an exact truncation (on baseline x86-64 it
is a libm call, and there were 540 million), and the arithmetic split from the
branches so that the compiler vectorises it. On one thread on the 300 image
sweep, carrying the profile onto the pixels went from 3.98 s to 2.16, profile
learning, which uses the same geometry, from 2.18 to 1.20, and the fit from 4.71
thread-seconds to 2.83. The corners' exact geometry, one call a pixel corner,
is now a quarter of what is left.

And a strong reflection's cells were computed twice, for learning and again for
fitting, from the same box, s1, phi and grid. In one pass the box learned from
is the box fitted, so its cells are kept with it: carrying the profile onto the
pixels went to 1.54 s, and the fit to 2.27 thread-seconds, 52 per cent less than
at the start, for 0.01 GB more held on the 300 image sweep. Byte-identical, in one
pass and in two, which computes them afresh.

## Indexing

The transform, the peak search and the refinement inside the macrocycles,
measured phase by phase with `mxi_index --timing`.

### Where indexing's time actually goes

Measured with `mxi_index --timing` on 78618 reflections of insulin, a 256^3
grid, one core:

    reciprocal points        0.008 s    0.1%
    max cell                 0.059 s    0.4%
    candidate vectors        6.676 s   50.2%
      the transform          3.213 s   24.2%
      the peak search        3.246 s   24.4%
      the rest of it         0.217 s    1.6%
    choose basis             4.097 s   30.8%
    fit and reduce           0.013 s    0.1%
    macrocycles              2.446 s   18.4%
    indexing total          13.304 s

**Both tables below were taken with a broken instrument** and the peak search
figures in them are about four times too large: the timer's closing assignment
had been placed inside the grid accessor, which the peak search calls some four
hundred and fifty million times, so every call read the clock. Corrected, on
the slower machine, the peak search is 0.598 s rather than 3.2 and the
transform and the basis search dominate. The tables are kept because the
reasoning built on them is instructive and because the correction is the point.

On another machine, same code and the same number of reflections:

    candidate vectors        3.843 s   78.8%
      the transform          2.783 s   57.1%
      the peak search        1.045 s   21.4%
    choose basis             0.183 s    3.8%
    macrocycles              0.788 s   16.2%
    indexing total           4.875 s

**The two runs disagree about which phase is largest**, and the disagreement is
not noise: `choose basis` is 3.8 per cent of one and 30.8 per cent of the
other, a factor of twenty-two, while the transform differs by fifteen per cent.

The reason is in the code. Choosing a basis scores every triple of candidate
vectors whose volume clears a degeneracy filter, so its cost is the number of
triples that clear it times the number of reflections. `--timing` reports the
count: 3743 of 4060 scored in the slow run, and something near 170 in the fast
one. How many survive depends on how nearly parallel the candidate vectors are,
which depends on the data. **The phase is data-dependent by more than an order
of magnitude and no single measurement of it means anything.**

The transform and the peak search are not: they are fixed work for a given grid
size, 16.7 million points either way, and they are 78.8 and 44.6 per cent of
the two runs. Those are the phases to attack.

What each is worth on the faster machine if it cost nothing:

    without the transform     2.33x
    without the peak search   1.27x
    without choose basis      1.04x
    without macrocycles       1.19x

#### What each phase is, as work

* **The transform**, a 256^3 complex FFT: 16.7 million points. A device library
  call -- cuFFT, or vDSP and MPS on Apple -- and milliseconds there. The one
  phase where the device version is someone else's code.
* **The peak search**: the modulus of 16.7 million voxels, then each compared
  with its twenty-six neighbours, then a sort. One thread per voxel, no
  communication, a reduction at the end. As good a fit for a device as exists.
* **Choosing the basis**: every triple of thirty candidate vectors scored
  against every reflection. Four thousand triples by seventy-eight thousand
  reflections, each independent. The same shape as the peak search and,
  measured here, the biggest single piece.
* **The macrocycles**: assignment, which is a pass over the reflections, and
  refinement, whose device port is designed in the rest of this document.

#### Threads before devices, and the transform is the one that matters

All of this is one core. From 4.875 s, with the peak search, the basis search
and the macrocycles threaded -- all three are embarrassingly parallel -- and
the transform treated two ways:

    4 threads, transform threads perfectly   1.28 s  3.8x
    4 threads, transform threads 2x only     1.97 s  2.5x
    8 threads, transform threads perfectly   0.68 s  7.2x
    8 threads, transform threads 2x only     1.72 s  2.8x

Threading everything except the transform hits a floor at about 1.7 s however
many cores are thrown at it, because the transform is 57 per cent of the run.
**Under a second needs the transform.** A threaded radix-2 will not get there
on its own; the realistic options are a library on the host -- vDSP, FFTW,
MKL -- or the device.

All of which is `std::thread` and no device, no memory transfers and no second
implementation to keep in step. The device is worth doing after that rather
than instead of it, and the honest comparison for any device backend is against
the threaded host version and not against this one.



An exploration, with the measurements that motivate it. Nothing here is
implemented yet.

### The transform, and the half of it that is not needed

With analytical derivatives in use the transform is the largest phase of
indexing, 42.5 per cent of a 1.069 second run. Two things are left in it.

**FFTW's own threading**, which `--fft-threads` now turns on where FFTW was
built with it. Debian ships `libfftw3_omp` and Homebrew `libfftw3_threads`;
either will do and the build reports which it found, or says it found neither
and runs on one thread.

**A real-to-complex transform.** The grid is filled by adding 1.0 at each
reciprocal lattice point and nothing else, so its imaginary part is zero
everywhere. A complex-to-complex transform of real data does twice the
arithmetic and holds twice the memory for an output that is Hermitian
symmetric: `F(-k)` is the conjugate of `F(k)`, so half of the 16.7 million
points computed are a reflection of the other half.

`fftw_plan_dft_r2c_3d` computes the `n * n * (n/2 + 1)` that are independent.
The obstacle is not the transform, it is the peak search, which walks the full
cube and compares each point with its twenty-six neighbours. On the half grid
some of those neighbours are the conjugates of points on the other side, and
the wrapping that the search already does for periodicity would have to become
a wrapping that also conjugates. That is a change to the part of this code
where an error would be least visible -- a peak list that is subtly wrong still
indexes, as the FFTW sign convention showed -- so it wants the agreement test
extended to the peak list itself before it is attempted, not afterwards.

### Where the refinement inside indexing spends its time

With the transform handed to FFTW, the rounding instruction enabled and the
clock out of the peak search, refinement is what is left. On a fast machine it
is 47.7 per cent of indexing; measured here, split three ways:

    macrocycles                 2.462 s   37.5%
      copy and select           0.181 s    2.8%
      refinement                2.239 s   34.1%
      reassignment              0.035 s    0.5%
        the jacobian            1.888 s   28.7%
        the normal equations    0.111 s    1.7%

**The Jacobian is 84 per cent of refinement**, which is what the rest of this
document assumed and had not shown. The normal equations are five per cent, and
the remaining eleven is residuals, outlier rejection and the solve.

So the thing to move is the Jacobian, and it has the right shape: one thread
per reflection and parameter, no communication, the analytical derivatives
already written and already validated against finite differences.

Two things about its current form that a port should not inherit:

* it is allocated fresh every iteration -- ten vectors of 118797 doubles on
  this data, 9.5 MB freed and reallocated per iteration -- where one buffer
  reused across iterations would do, and on a device must;
* it is parameter-major, a vector per parameter over all reflections. That is
  the wrong way round for a device, where the reflection is the thread and the
  parameters of one reflection want to be adjacent.

Neither is worth changing on the host for its own sake without measuring what
the allocation costs. Both are worth knowing before writing the device version,
because the layout is the part that is expensive to change afterwards.

## Refinement

The refinement target is where a device port was first designed, and the
design is here -- why the target is the thing to move, the shape of its work,
its precision in float32, the analytical derivatives it needs -- followed by
what scan-varying refinement and file input then turned out to cost on the
CPU. The port is designed and not written.

### It is the only thing worth moving

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

### The shape of the work

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

### The B-spline makes the Jacobian banded

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

### Precision: measured, and the earlier estimate was wrong

Apple GPUs have no double precision at all. `src/target.hh` is the whole target
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

#### The estimate that was wrong, and why

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

#### What follows

Analytical derivatives are not a nicety for a device port, they are a
precondition. That reverses the earlier plan, which had them fourth on the list
as an optimisation.

`tests/test_precision.cc` asserts the failure as well as the success, so that
nobody later assumes numerical differentiation would port as it stands.

#### The analytical derivative in float32: measured

`src/derivatives_t.hh` is the same treatment applied to the derivatives.
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

### Analytical derivatives: written, and validated

`src/derivatives.hh` implements Appendix A of Waterman et al. (2016) for the
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

#### Where the two Jacobians disagree, and why it is not the analytical one

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

### The volume cutoff, and a guess it did not support

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

### Scan-varying refinement, where the normal equations were the cost

The static case measured above, under "It is the only thing worth moving", is
not the command anyone actually runs. This is:

    mxi_refine indexed.expt indexed.refl --scan-varying 18 --beam --analytic

Eighteen control points takes the parameter count from about ten to about a
hundred and seventy, and the normal equations are quadratic in it. Measured:

    the jacobian               1.791 s   23.5%
    the normal equations       4.848 s   63.6%
    total                      7.629 s

A scan-varying crystal is a cubic B-spline, so a reflection touches four
control points and no more: about forty of those hundred and seventy
parameters have a nonzero derivative and the rest are structurally zero. The
accumulation skipped the zeros in its outer index and not its inner one, so
each nonzero was multiplied against every parameter below it, and
`jacobian[b][row]` walked across a hundred and seventy separate heap
allocations -- one load per parameter, which is the worse half of the cost.

Gathering the nonzero entries of a row once and using them against each other:

    the normal equations       4.848 -> 1.517 s
    total                      7.629 -> 4.204 s

The order of accumulation is unchanged, ascending in both indices, so every sum
is formed from the same terms in the same sequence. The refined crystal,
detector and beam are identical, compared as written files rather than to a
tolerance.

That the sparsity was there to be used is a property of the model rather than
of the data, so it holds for any scan-varying refinement and the saving grows
with the number of control points.

#### The search for the nonzeros was itself linear in the parameters

Using the sparsity meant finding it, and the finding was done by reading every
entry of every column of the Jacobian for each row:

    for (a = 0; a < n; ++a) { ja = jacobian[a][row]; if (ja == 0) continue; ... }

That is O(n) per row however sparse the row is, and each probe lands in a
different heap allocation. At a hundred and seventy parameters it did not
matter. On ten full rotations of a crystal, where the scan-varying model runs
to thousands of parameters, it was 267 of the 291 seconds the refinement took.

The pattern does not need finding. `build_analytic_jacobian` knows which
parameters it is about to write -- four control points from the spline, six
detector, two beam -- so it records the span as it goes. The accumulation then
iterates thirty-six parameters instead of several thousand.

    control points      before      after
              18       1.748 s     1.134 s
              60       5.180 s     1.204 s
             120      10.602 s     1.573 s

Linear in the control points before, nearly flat after. The refined crystal and
detector are identical.

The finite-difference path records nothing, because it genuinely does not know
which parameters it touched, and falls back to the search. It is slower than
the analytical path by a much larger factor anyway.

#### The Jacobian threads at 1.29x, and that is an allocation

On sixteen cores, a scan-varying refinement of ten rotations:

    --jacobian-threads 1     the jacobian  17.714 s
    --jacobian-threads 16    the jacobian  13.722 s

The inner loop is one reflection per thread with nothing shared, so 1.29x is
not a threading problem. It is what surrounds the loop:

    jacobian->assign(n, std::vector<double>(observations.size() * 3, 0.0));

Parameters by reflections by three doubles, allocated and zeroed every
iteration. At a million reflections that is gigabytes of serial memory traffic
wrapped around a parallel inner loop, and no thread count touches it.

The fix is not to allocate it. Every entry is written once and read once, by
the accumulation that immediately follows, so the two can be fused: compute a
reflection's derivatives and accumulate them into the normal equations there
and then, and the array never exists. That removes the allocation, the traffic
in both directions, and the parameter-major layout already noted as wrong for a
device. It is what a device port has to do anyway, since gigabytes an iteration
is not a thing to move across a bus.

#### And then they were the only serial phase left

On a sixteen-core machine, with the Jacobian threaded and the sparsity used:

    read                       0.261 s   17.3%
    the jacobian               0.085 s    5.6%
    the normal equations       0.884 s   58.4%
    total                      1.513 s

The Jacobian is five per cent there and forty per cent on a single core, which
is not a statement about the Jacobian: it is threaded and the normal equations
were not. `--normal-threads` threads them.

It is a reduction, so this one is not free. Each thread accumulates into its
own `n` by `n` matrix -- 231 kB apiece at a hundred and seventy parameters,
which is why the thread count is capped by what the partials cost -- and the
partials are summed at the end. **A threaded run and a serial one do not agree
bit for bit**, because floating point addition is not associative and the terms
are summed in a different order. They agree to about 1.6e-12 relative, seven
orders below the convergence tolerance, and the test pins 1e-10.

Two threaded runs at the same setting do agree exactly: the chunk boundaries
come from the thread count, not from how the threads happen to be scheduled. A
parallel reduction that answered differently run to run would make every
comparison downstream meaningless, so that is tested as well.

`--normal-threads 1` keeps the serial sum where the last bits matter.

### Refinement on its own is dominated by file I/O

`mxi_refine --timing`, on 78618 reflections whose table carries shoeboxes:

    read                       0.293 s   28.7%
    build the target rows      0.032 s    3.1%
    the jacobian               0.251 s   24.5%
    the normal equations       0.111 s   10.9%
    the solve                  0.000 s    0.0%
    the trial residuals        0.110 s   10.7%
    outlier rejection          0.008 s    0.8%
    write                      0.130 s   12.7%
    total                      1.023 s

**Reading and writing are 41 per cent of it**, more than the Jacobian, and
nine tenths of that is pixel data refinement never looks at. The file is
100 MB, of which the shoebox column is 79. The same run on a table with the
shoeboxes stripped:

    read     0.293 -> 0.049 s
    write    0.130 -> 0.042 s
    total    1.023 -> 0.615 s

A third of the run was carrying shoeboxes from one file to another, and almost
none of it was the file. Measured against the floor:

    fread of the same 100 MB       0.055 s   at 1.8 GB/s
    memcpy of 100 MB               0.012 s
    read_reflections               0.310 s

The reader had

    std::string raw((std::istreambuf_iterator<char>(in)), {});

which goes through the stream one character at a time, regrowing the string as
it goes. Sized, resized and read in one go it is 0.148 s. The writer was
growing its output from empty to a hundred megabytes, doubling and copying
everything it had each time; reserved, 0.138 s.

    read      0.310 -> 0.148 s
    write     0.163 -> 0.138 s

A table read and written back is byte for byte what it was, which is the only
check worth making on a change to framing: values that survive a round trip
would survive most ways of getting the framing wrong.

Which is worth saying plainly: refinement's arithmetic is now a smaller part of
`mxi_refine` than its file handling, and no amount of threading or device work
on the Jacobian will change that. The cheapest remaining second is in the
reader.

The solve is 0.000 s and that is not a broken timer: it is a Cholesky of a ten
by ten matrix, done a handful of times.
