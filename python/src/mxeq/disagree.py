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
    factor: float = 2.0,
    absolute: float = 0.0,
    radius: float = 0.5,
    limit: int = 0,
) -> tuple[refl.ReflectionTable, str]:
    """Rows of `ours` whose `value` disagrees with `theirs` by more than
    `factor`, either way round.

    `absolute` skips pairs where both are small, since a ratio between two
    numbers near zero is meaningless and there are a great many of them: a
    reflection of 0.6 against 0.2 is a factor of three and is nothing at all.
    """
    for name, table in (("ours", ours), ("theirs", theirs)):
        for column in ("miller_index", "entering", "xyzcal.px", value):
            if column not in table.columns:
                raise ValueError(f"{name} has no {column}")

    matching, duplicates = match.match_keys(
        _key(ours),
        _key(theirs),
        tie_break_a=ours.columns["xyzcal.px"][:, 2],
        tie_break_b=theirs.columns["xyzcal.px"][:, 2],
    )
    ia, ib = matching.index_a, matching.index_b
    close = (
        np.abs(ours.columns["xyzcal.px"][ia, 2] - theirs.columns["xyzcal.px"][ib, 2])
        <= radius
    )
    ia, ib = ia[close], ib[close]
    if len(ia) == 0:
        raise ValueError(
            "nothing matched: if one file was written before a reindexing step "
            "and the other after, the Miller indices are not comparable"
        )

    va = ours.columns[value].ravel()[ia]
    vb = theirs.columns[value].ravel()[ib]
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = va / vb
    big_enough = (np.abs(va) > absolute) | (np.abs(vb) > absolute)
    # Either way round, and a sign flip counts: a reflection one program calls
    # positive and the other negative disagrees however small the numbers.
    apart = (
        ~np.isfinite(ratio)
        | (ratio > factor)
        | (ratio < 1.0 / factor)
        | (np.sign(va) != np.sign(vb))
    )
    chosen = np.where(apart & big_enough)[0]
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
        + (f", {duplicates} duplicate keys" if duplicates else "")
        + f"\n{len(rows)} disagree on {value} by more than {factor}x"
        + (f" (of {int((apart & big_enough).sum())} found)" if limit else "")
    )
    if len(rows):
        px = out.columns["xyzcal.px"]
        report += (
            f"\n  frames {px[:, 2].min():.0f} to {px[:, 2].max():.0f}"
            f", median ratio {np.nanmedian(ratio[chosen]):.3g}"
        )
    return out, report
