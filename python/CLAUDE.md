# mxeq -- working notes

Read this before changing anything. It records constraints that were arrived at
rather than chosen, and the reasons are not recoverable from the code.

## What this is for

Checking that a GPU-accelerated MX pipeline produces output **equivalent** to
DIALS', boundary by boundary. DIALS is the acceptance oracle. `mx` is the
development oracle. Neither is imported here.

## Hard constraints

**No cctbx, no DIALS, no dxtbx imports. Ever.** The whole point is that this
runs in a container with no crystallographic software in it, against files
written by a DIALS build somewhere else. Reaching for `dials.array_family` to
read a `.refl` would be more convenient and would destroy the tool.

**No thresholds in version one.** Every check reports distributions. A
tolerance chosen before the distribution has been seen on real data is a
tolerance chosen by guessing, and it will then be tuned until it passes, which
is worse than having none. Thresholds go in a config file against the `--json`
output once `ins10_1` has been through it.

**Equivalence, not identity.** Identity is the right target for the threshold
kernel and is achievable there. It is not achievable at integration and
demanding it would report every difference in summation order as a failure.

**Exit status is not a verdict.** Zero means the comparison ran. If that ever
changes, the change must be explicit and documented here, because scripts will
have been written against the current meaning.

## Format -- validated

The `.refl` layout is now pinned against real files (`dials.find_spots` and
`mxi_find` output, 13766 rows, ten columns including a shoebox):

```
[ "dials::af::reflection_table", 2, { identifiers, nrows, data } ]
data[column] = [ type_name, [ nrows, blob ] ]
```

Text is msgpack `str`, payloads are `bin`, blobs are little-endian with the
components of a compound type adjacent. Element widths: `int` 4, `std::size_t`
8, `double` 8, `vec3<double>` 24, `int6` 24. The count inside a payload is the
number of **rows**, not scalars, and is checked rather than trusted.

`tests/data/` holds a 48-row cut of real output, sliced by
`tests/make_real_fixture.py` rather than written by `refl.dumps` -- see below
for why that distinction is the whole point. Its first 31 bytes are identical
to the real file's.

The identifiers key is read under both `identifiers` and
`experiment_identifiers` because DIALS has spelled it more than one way. Both
the two- and three-element top level, and wrapped or bare payloads, are
accepted: it costs nothing, and a reader that understands only the one file it
was tested against is not much of a reader.

## Never join on a millimetre column

`xyzobs.mm.value` is computed once at import, carries the inverse parallax
correction under the geometry current at that moment, and is **never
recomputed**. After refinement it describes a detector that no longer exists.

Two pipelines that refined to slightly different detectors would then be
compared in two different millimetre frames, and the difference would read as a
centroid disagreement. On the insulin data the stale-versus-current gap is
1.7e-3 mm. Pixels are the raw measurement and cannot go stale, so
`POSITION_COLUMNS` holds only pixel columns and `position_column` raises with a
pointed message if a table has nothing but millimetres.

The third component is the exception: it is the rotation angle, which comes
from the scan and has no dependence on the detector model, so it cannot go
stale. That is why it is usable for pinning the scan convention and the first
two are not.

## Things that were got wrong once

**The format was wrong in two ways and the whole suite passed anyway.** The
first reader assumed a two-element top level (it is three, with a version) and
a bare payload blob (it is a `[count, blob]` pair). Seventy tests passed,
because the only thing checking the format was a round-trip against this
package's own writer -- and a writer built from the same wrong assumption is
perfectly consistent with it.

A round-trip proves self-consistency. **Only a file written by something else
proves a format.** That is why `tests/data` exists, why the fixture is cut by
slicing a real document's blobs rather than by calling `refl.dumps`, and why
`tests/test_format.py` asserts on raw header bytes and element widths instead
of on decode-then-encode.

**The escape hatch failed on the first real file it saw.** `refl.structure()`
was documented as never raising, and raised `UnicodeDecodeError` -- it decoded
every `bin` to test whether it was printable text, and a column of packed
doubles is not valid UTF-8. A diagnostic that only works on data that was
already readable is not a diagnostic. It now has a `printable()` helper that
returns `None` instead of raising, and a blanket `except` around the walk.

**The A matrix is the plain inverse of the real-space matrix, not the inverse
transpose.** The real-space matrix has a, b, c as its *rows*, so `M A = I` by
the definition of the reciprocal basis. The inverse transpose passes every test
on an orthogonal cell, because both are diagonal there. `test_readers.py` uses a
triclinic cell for this and nothing else.

**The candidate operator pool must be closed under transposition explicitly.**
It was originally claimed that the pool was closed already, because it is built
from lattice groups. The cubic group is -- its metric is the identity, so
`M^T = M^-1`. The hexagonal group is not, for the same reason in reverse. A
test caught it. Left unclosed, the reindexing search would have been
convention-dependent for hexagonal and rhombohedral lattices only.

**`fraction_matched` divides by the larger table, not the smaller.** Otherwise a
table containing one reflection that happens to be in the other scores 1.0.

**The `.expt` scan shape changed and the old reader failed silently.** Current
DIALS writes `properties.oscillation` as a per-image array of start angles, not
a top-level `[start, width]`. Reading only the old key against a current file
returned `(0.0, 0.0)` — for *both* sides of a comparison, which then agreed
perfectly and said nothing. A silent pass is worse than a failure. Both forms
are now read, pinned against a real scan block in `tests/data/real_scan.json`,
and the width comes from the endpoints rather than the first two elements
because the array is built by repeated addition.

## Design decisions worth not relitigating

**One module per boundary, not one parameterised function.** The join, the
metrics and the failure modes differ at every boundary. Sharing a body between
them would put the differences in flags.

**The reindexing operator is found from the data, never assumed, and never
applied silently.** It is found at the `indexed` boundary and passed explicitly
to later ones with `--operator`. A later boundary given no operator when one is
needed reports a low match fraction and says so; it does not go looking.

**Duplicate keys are resolved by frame order.** Never by intensity: that would
pair whichever two happened to agree.

**CC-half splits randomly with a fixed seed, not on observation parity.**
Observation order is pipeline-dependent, and splitting on it would make CC-half
compare two pipelines' orderings rather than their data.

**Resolution comes from the `.expt` in preference to the `d` column.** `d` is a
pipeline output; binning by it bins the comparison by one of the things being
compared.

## Testing

`tests/fixtures.py` generates synthetic data with *planted* answers -- a known
reindexing operator, a known intensity ratio, a known misorientation. A test
that asserts against a number nobody chose is a test that passes for reasons
nobody knows.

Fixtures are imported by name via `tests/conftest.py`, not shipped in the
package. They generate test data and have no business being installed.

## Style

`black` formatting. Four-space indent, double quotes. No `from x import *`.
Comments explain why, not what; a comment restating the line above it is
deleted, not improved.

## Say which file was wrong, not which byte

`check refined` compares models and the rest compare reflections. Handing the
wrong kind failed inside a UTF-8 decoder:

    UnicodeDecodeError: 'utf-8' codec can't decode byte 0x93 in position 0

which names the byte and not the mistake. 0x93 is msgpack's header for a
three-element array, so a reflection table announces itself in its first byte
and the message can say so.

A pair of reflection tables passed to `check refined` is now routed to the
reflection comparison rather than refused, because comparing two refined
pipelines usually does mean comparing their reflections. Anything genuinely
mismatched is refused with what it is and what was wanted.

## Two residuals over different reflections do not compare

DIALS refined 1800 images of insulin on 61679 of 78618 reflections; this
package on 71454. The headline rmsds are not answering the same question, and
the difference between them is partly a difference in which reflections were
averaged.

`check indexed` now reports both sides over the reflections both predicted. On
insulin, over the same 75343:

    DIALS   median 0.3307 px   99% 1.517
    ours    median 0.3997 px   99% 1.785

So the gap is real and about a fifth, rather than an artefact of the set.
