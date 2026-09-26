"""The reflections two integrations disagree about, as a table to look at.

A trend says a disagreement exists and a table of worst offenders says which
reflections they are.  Neither shows the pixels, and the pixels are where the
answer is: whether a strong neighbour is leaking into the box, whether the spot
is somewhere other than where it was predicted, whether the box has a module
gap through it.

So this writes a reflection table holding only the ones that disagree, which
`dials.image_viewer` will draw on the images:

    dials.image_viewer imported.expt disagree.refl

The rows are OURS, not the reference's, because the question is what our boxes
did.  Their value is carried along in `reference.intensity` and the ratio in
`disagreement.ratio`, so the viewer's table shows both.
"""

from __future__ import annotations

import numpy as np

from . import match, refl


def _key(table: refl.ReflectionTable) -> list[np.ndarray]:
    hkl = table.columns["miller_index"]
    entering = table.columns["entering"].ravel().astype(int)
    return [hkl[:, 0], hkl[:, 1], hkl[:, 2], entering]


def select(
    ours: refl.ReflectionTable,
    theirs: refl.ReflectionTable,
    value: str = "intensity.sum.value",
    value_b: str | None = None,
    factor: float | None = None,
    difference: float | None = None,
    sigmas: float | None = None,
    floor: float = 0.0,
    radius: float = 5.0,
    limit: int = 0,
) -> tuple[refl.ReflectionTable, str]:
    """Rows of `ours` whose `value` disagrees with `theirs`.

    Three ways of asking, because they find different things:

    * `factor` -- a RELATIVE difference, `a/b` outside [1/factor, factor].
      Finds what is proportionally wrong, which is most things, and is
      dominated by weak reflections where a ratio means least.
    * `difference` -- an ABSOLUTE difference, `|a - b|` in counts.  Finds what
      is wrong by a lot, which on a constant offset is the strong reflections
      and on a scale error the strong ones too.  A background biased by a
      fraction of a count costs every reflection the same number of counts, and
      only this sees it as one thing.
    * `sigmas` -- the difference in units of the two variances added, which is
      the only one of the three that knows whether a disagreement is larger
      than the measurement.  Needs the matching variance column.

    Each given criterion must be exceeded, so passing two narrows rather than
    widens.  With none given, `factor` defaults to 2.

    `floor` skips pairs where both values are below it: a ratio between two
    numbers near zero means nothing and there are a great many of them, a
    reflection of 0.6 against 0.2 being a factor of three and nothing at all.
    It does not apply to `difference` or `sigmas`, which are not confused by
    small numbers in the first place.
    """
    if factor is None and difference is None and sigmas is None:
        factor = 2.0
    # A different column on each side, so one file can be compared with itself:
    # our summed intensity against our fitted one is the comparison that says
    # which of the two is misbehaving, and it needs no second program.
    other = value_b or value
    for name, table, column_name in (
        ("ours", ours, value),
        ("theirs", theirs, other),
    ):
        for column in ("miller_index", "entering", "xyzcal.px", column_name):
            if column not in table.columns:
                raise ValueError(f"{name} has no {column}")

    from .trends import paired

    ia, ib, unpartnered = paired(ours, theirs, radius)

    va = ours.columns[value].ravel()[ia]
    vb = theirs.columns[other].ravel()[ib]
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = va / vb

    apart = np.ones(va.shape, dtype=bool)
    asked = []
    if factor is not None:
        big_enough = (np.abs(va) > floor) | (np.abs(vb) > floor)
        # Either way round, and a sign flip counts: a reflection one program
        # calls positive and the other negative disagrees however small the
        # numbers.
        relative = (
            ~np.isfinite(ratio)
            | (ratio > factor)
            | (ratio < 1.0 / factor)
            | (np.sign(va) != np.sign(vb))
        )
        apart &= relative & big_enough
        asked.append(f"a factor of {factor:g}")
    if difference is not None:
        apart &= np.abs(va - vb) > difference
        asked.append(f"{difference:g} counts")
    if sigmas is not None:
        variance_a = value.replace(".value", ".variance")
        variance_b = other.replace(".value", ".variance")
        if variance_a not in ours.columns or variance_b not in theirs.columns:
            raise ValueError(
                f"--sigma needs {variance_a} in ours and {variance_b} in "
                "theirs, and one of them is missing"
            )
        combined = np.sqrt(
            np.maximum(ours.columns[variance_a].ravel()[ia], 0.0)
            + np.maximum(theirs.columns[variance_b].ravel()[ib], 0.0)
        )
        with np.errstate(invalid="ignore", divide="ignore"):
            pull = np.abs(va - vb) / combined
        apart &= np.isfinite(pull) & (pull > sigmas)
        asked.append(f"{sigmas:g} sigma")

    chosen = np.where(apart)[0]
    if limit and len(chosen) > limit:
        # The worst, not the first: an arbitrary truncation would hide whatever
        # is at the far end of the scan.
        order = np.argsort(-np.abs(np.log10(np.abs(ratio[chosen]) + 1e-30)))
        chosen = chosen[order[:limit]]
        chosen.sort()

    rows = ia[chosen]
    out = refl.ReflectionTable(nrows=len(rows))
    for name, column in ours.columns.items():
        out.columns[name] = column[rows]
        if name in ours.types:
            out.types[name] = ours.types[name]
    out.identifiers = dict(ours.identifiers)
    # What the other program made of the same reflection, so the viewer's table
    # shows both without another file being opened.
    out.columns["reference.intensity"] = vb[chosen].astype(float)
    out.types["reference.intensity"] = "double"
    out.columns["disagreement.ratio"] = ratio[chosen].astype(float)
    out.types["disagreement.ratio"] = "double"

    report = (
        f"{ours.nrows} rows against {theirs.nrows}, matched {len(ia)}"
        + (f", {unpartnered} unpartnered" if unpartnered else "")
        + f"\n{len(rows)} disagree on {value}"
        + (f" against {other}" if other != value else "")
        + " by more than "
        + " and ".join(asked)
        + (f" (of {int((apart & big_enough).sum())} found)" if limit else "")
    )
    if len(rows):
        px = out.columns["xyzcal.px"]
        report += (
            f"\n  frames {px[:, 2].min():.0f} to {px[:, 2].max():.0f}"
            f", median ratio {np.nanmedian(ratio[chosen]):.3g}"
        )
    return out, report
