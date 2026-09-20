# Why each join is the way it is

Every check in this package is a join followed by a diff. The diff is easy.
The join is where a comparison silently stops measuring what it claims to, and
each boundary gets it wrong in its own way.

## strong

**There is no key.** Two spot finders produce two lists of positions in
whatever order their grouping emitted them, and they need not be the same
length. So the join is spatial: mutual nearest neighbour within a radius.

Mutual matters. Without it, a cluster of three spots in A all claim the same
row of B and the matched count exceeds what actually corresponds.

The radius default is generous (2 px). A tight radius does not report
disagreement as disagreement; it reports it as absence, moving spots out of the
matched population and into the unmatched count where their offsets can no
longer be seen.

**The interesting output is the unmatched population, not the matched one.**
Two finders that agree on 95% of spots with sub-pixel centroids, differing only
on the weakest 5%, are equivalent for every practical purpose. Two that match
100% but with a systematic quarter-pixel shift are not. Only the split into
matched, only-in-A and only-in-B separates those, which is why the check
describes the intensity distribution of each separately: if the unmatched spots
are the weak ones, the two finders agree and differ in where they cut.

**Spot depth is reported because a three-dimensional grouping bug hides here.**
A reflection swept through a degree lands on several images. Grouped per image
instead of six-connected in x, y and z, it becomes three spots with three
centroids, and nothing downstream reassembles them. That failure shows up as a
depth distribution pinned at one image, and otherwise not until indexing
quietly does worse.

## indexed

**The order of operations is the whole content of this check.** Match
spatially, then find the reindexing operator over the matched pairs, then
compare Miller indices. Compare the indices first and two perfectly good
solutions related by a change of basis report near-total disagreement.

The operator is found from the data, not assumed: both tables describe the same
observed spots, so a spatial match gives a correspondence, and over that
correspondence the right operator is the one under which nearly all the indices
agree. The candidate pool is the cubic and hexagonal lattice groups, explicitly
closed under transposition and inversion.

That closure is not decoration. The cubic group is closed under transposition
already, because its metric is the identity and so `M^T = M^-1`. The hexagonal
group is not, for the same reason in reverse. Without the explicit closure the
search would depend on whether the operator acts on Miller indices as a row or
a column vector, and would do so **for hexagonal and rhombohedral lattices
only** -- which is how a convention bug survives every test case someone thinks
to write down.

`(0, 0, 0)` rows are excluded from the search. An unindexed reflection agrees
under every operator, so a mostly-unindexed table would otherwise score well
for whichever candidate happened to come first.

The check also reports whether the winning operator is a symmetry of the cell
metric. One that is not cannot be a reindexing, and a high agreement for it is a
result to look at rather than to use.

## refined

**Refined parameters are not comparable pointwise and this check does not try.**
The same lattice can be described in a different basis, the same orientation
reached by a different rotation, and a scan-varying model can have a different
number of intervals depending on what the refinement chose. Comparing parameter
vectors would compare parameterisations.

What is compared is the cell, the misorientation, and the geometry.

The misorientation comes from the polar decomposition of `A_b M A_a^-1`,
minimised over the candidate lattice operators. The minimisation is the point:
without it, two models related by a reindexing report ninety degrees and the
check is useless on exactly the case it exists for. The polar decomposition
rather than the trace directly, because when the two cells differ slightly the
matrix is not orthogonal and `arccos` of its trace means nothing.

Panel axis angles are reported alongside the origin shift because a rotation of
a panel about its own centre leaves the distance and very nearly the origin
alone.

A differing image range raises a warning rather than an error: it is legitimate,
and it offsets every z comparison downstream by a constant that is not a bug.

## integrated

**The key exists, so the difficulty moves to the metric.**

A percentage difference in intensity is meaningless averaged over a dataset
spanning five orders of magnitude. A single correlation coefficient over
unbinned intensities is dominated by the hundred strongest reflections and can
sit at 0.9999 while the outer shell is systematically out by four percent --
which is not hypothetical, it is what the fixture in `tests/` is built to
demonstrate. So everything is reported in equal-volume resolution shells, and
Spearman is printed next to Pearson.

Shells are equal in reciprocal volume rather than equal in count, so that a
shell means the same thing between two runs that selected different numbers of
reflections.

**The pull needs its warning read.** `(I_a - I_b) / sqrt(s_a^2 + s_b^2)` is the
right statistic for two *independent* measurements. These are the same photons
processed twice; the errors are correlated, and the width will be far below one
even when the disagreement is real. A width of 0.2 is not a pass. Use the pull
for its shape -- symmetry, and tail weight -- and the relative difference for
its size.

**Duplicate keys are resolved by frame, not by intensity.** On a scan that
turns through more than 360 degrees the same reflection is recorded twice under
the same key. Pairing them by intensity would pair whichever two happened to
agree and hide the disagreement it was called to find.

Resolution is taken from the experiment model where one is given, in preference
to the table's `d` column, because `d` is itself a pipeline output and binning
by it would bin the comparison by one of the things being compared.

## scaled

**Equivalence stops being about individual reflections here.** Two scaling runs
can disagree about every inverse scale factor and produce merged data that
agrees to within the noise, because an overall scale and a smooth function of
the scan are not observable in the merged result. Insisting the scale factors
match would be insisting the two models be the same model.

So the comparison is on what survives that freedom: the internal consistency of
each dataset alone (CC-half, Rmeas), and the correlation between the two merged
sets after both are mapped into the same asymmetric unit.

**Read `cc_merged` against the two CC-half values, not against 1.** Two
datasets cannot correlate with each other better than their own internal
consistency permits. A `cc_merged` close to the CC-half values is agreement to
within the noise, and that is the result being looked for.

The half-dataset split for CC-half is random with a fixed seed rather than by
parity of the observation index. Observation order is a pipeline-dependent
thing, and splitting on it would make CC-half compare two pipelines' orderings
instead of their data.

The scale factors are still reported -- as a ratio, with the note that a
constant ratio is an unobservable overall scale and it is the *spread* that
says the two models differ.
