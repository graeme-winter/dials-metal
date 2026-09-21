"""Where two integrations disagree, and against what.

A single correlation says two columns disagree.  It does not say whether the
disagreement is with resolution, with intensity, with position on the detector,
with how close a reflection sits to the rotation axis, or with how well the
profile describes it -- and those point at different causes.  This bins a ratio
against each of those in turn, so the shape of the disagreement is visible.

Read the tables rather than the medians.  A ratio that is flat everywhere and
merely noisy is a different problem from one that trends with resolution, and
the latter is usually geometry while the former is usually weighting.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import match, refl, stats


@dataclass
class Explanatory:
    """One thing to bin against."""

    name: str
    values: np.ndarray
    unit: str = ""
    #: Highest first, as resolution shells are conventionally shown.
    descending: bool = False


def _beam_centre(table: refl.ReflectionTable) -> tuple[float, float]:
    """The centre the reflections sit around, for a radius with no experiment
    list to hand.  The median is used rather than the mean because a few
    reflections at the edge of a detector would otherwise drag it."""
    px = table.columns["xyzcal.px"]
    return float(np.median(px[:, 0])), float(np.median(px[:, 1]))


def explanatory_variables(
    a: refl.ReflectionTable, index: np.ndarray, reference: refl.ReflectionTable, ref_index: np.ndarray
) -> list[Explanatory]:
    """The variables worth binning a disagreement against.

    Taken from the REFERENCE where both have the column, so that the bins do
    not move when the thing being investigated changes.
    """
    out: list[Explanatory] = []

    def column(table: refl.ReflectionTable, idx: np.ndarray, name: str):
        if name not in table.columns:
            return None
        return table.columns[name][idx]

    d = column(reference, ref_index, "d")
    if d is not None:
        out.append(Explanatory("resolution", d.ravel(), "A", descending=True))

    signal = column(reference, ref_index, "intensity.sum.value")
    variance = column(reference, ref_index, "intensity.sum.variance")
    if signal is not None and variance is not None:
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = signal.ravel() / np.sqrt(np.maximum(variance.ravel(), 1e-12))
        out.append(Explanatory("I/sigma (reference)", ratio))

    zeta = column(reference, ref_index, "zeta")
    if zeta is not None:
        out.append(Explanatory("|zeta|", np.abs(zeta.ravel())))

    part = column(reference, ref_index, "partiality")
    if part is not None:
        out.append(Explanatory("partiality", part.ravel()))

    cc = column(a, index, "profile.correlation")
    if cc is not None:
        out.append(Explanatory("profile.correlation (ours)", cc.ravel()))

    fg = column(a, index, "num_pixels.foreground")
    if fg is not None:
        out.append(Explanatory("foreground pixels (ours)", fg.ravel().astype(float)))

    back = column(a, index, "background.mean")
    if back is not None:
        out.append(Explanatory("background (ours)", back.ravel(), "counts/pixel"))

    px = column(reference, ref_index, "xyzcal.px")
    if px is not None:
        cx, cy = _beam_centre(reference)
        radius = np.hypot(px[:, 0] - cx, px[:, 1] - cy)
        out.append(Explanatory("detector radius", radius, "px"))
        out.append(Explanatory("frame", px[:, 2], ""))

    return out


def bin_edges(values: np.ndarray, n_bins: int) -> np.ndarray:
    """Equal-population bins, so every row of the table carries the same
    weight.  Equal-width bins on a quantity like I/sigma put almost everything
    in the first bin and tell you nothing."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.array([0.0, 1.0])
    quantiles = np.linspace(0.0, 100.0, n_bins + 1)
    edges = np.unique(np.percentile(finite, quantiles))
    if edges.size < 2:
        edges = np.array([finite.min(), finite.max() + 1e-9])
    return edges


@dataclass
class Row:
    low: float
    high: float
    count: int
    median_ratio: float
    correlation: float
    median_a: float
    median_b: float


def trend(
    value_a: np.ndarray,
    value_b: np.ndarray,
    against: Explanatory,
    n_bins: int = 10,
) -> list[Row]:
    """Bin the ratio of two quantities against a third."""
    ok = (
        np.isfinite(value_a)
        & np.isfinite(value_b)
        & np.isfinite(against.values)
        & (np.abs(value_b) > 0)
    )
    a, b, x = value_a[ok], value_b[ok], against.values[ok]
    edges = bin_edges(x, n_bins)
    rows: list[Row] = []
    for i in range(len(edges) - 1):
        low, high = edges[i], edges[i + 1]
        last = i == len(edges) - 2
        inside = (x >= low) & ((x <= high) if last else (x < high))
        if inside.sum() < 2:
            continue
        ratio = a[inside] / b[inside]
        r, _ = stats.correlation(a[inside], b[inside])
        rows.append(
            Row(
                low=float(low),
                high=float(high),
                count=int(inside.sum()),
                median_ratio=float(np.median(ratio)),
                correlation=float(r),
                median_a=float(np.median(a[inside])),
                median_b=float(np.median(b[inside])),
            )
        )
    if against.descending:
        rows.reverse()
    return rows


def format_trend(against: Explanatory, rows: list[Row], value: str) -> str:
    unit = f" ({against.unit})" if against.unit else ""
    lines = [f"{value} against {against.name}{unit}"]
    lines.append(
        f"  {'from':>10} {'to':>10} {'n':>8} {'ours':>10} {'theirs':>10} "
        f"{'ratio':>8} {'corr':>7}"
    )
    for r in rows:
        lines.append(
            f"  {r.low:10.4g} {r.high:10.4g} {r.count:8d} "
            f"{r.median_a:10.4g} {r.median_b:10.4g} "
            f"{r.median_ratio:8.4f} {r.correlation:7.4f}"
        )
    return "\n".join(lines)


def compare(
    a: refl.ReflectionTable,
    b: refl.ReflectionTable,
    values: list[str],
    n_bins: int = 10,
    radius: float = 0.5,
) -> str:
    """Match two integrated tables and report every trend for every value."""
    needed = ("miller_index", "entering", "xyzcal.px")
    for name, table in (("A", a), ("B", b)):
        missing = [c for c in needed if c not in table.columns]
        if missing:
            raise ValueError(f"{name} has no {', '.join(missing)}")

    def key(table: refl.ReflectionTable) -> list[np.ndarray]:
        hkl = table.columns["miller_index"]
        entering = table.columns["entering"].ravel().astype(int)
        return [hkl[:, 0], hkl[:, 1], hkl[:, 2], entering]

    matching, duplicates = match.match_keys(
        key(a),
        key(b),
        tie_break_a=a.columns["xyzcal.px"][:, 2],
        tie_break_b=b.columns["xyzcal.px"][:, 2],
    )
    ia, ib = matching.index_a, matching.index_b
    # A key can match while the two rows are different observations of it, on
    # a scan that goes round more than once.  The frame settles that.
    close = np.abs(a.columns["xyzcal.px"][ia, 2] - b.columns["xyzcal.px"][ib, 2]) <= radius
    ia, ib = ia[close], ib[close]

    out = [
        f"{a.nrows} rows against {b.nrows}, matched {len(ia)}"
        + (f", {duplicates} duplicate keys" if duplicates else "")
    ]
    if len(ia) == 0:
        out.append("")
        out.append(
            "Nothing matched.  If one file was written before a reindexing step "
            "and the other after, the Miller indices are not comparable and the "
            "geometry has to be pinned to the same crystal."
        )
        return "\n".join(out)

    variables = explanatory_variables(a, ia, b, ib)
    for value in values:
        if value not in a.columns or value not in b.columns:
            out.append("")
            out.append(f"{value}: not in both files")
            continue
        va = a.columns[value].ravel()[ia]
        vb = b.columns[value].ravel()[ib]
        overall, _ = stats.correlation(va, vb)
        finite = np.isfinite(va) & np.isfinite(vb) & (np.abs(vb) > 0)
        out.append("")
        out.append(
            f"=== {value}: correlation {overall:.4f}, "
            f"median ratio {np.median(va[finite] / vb[finite]):.4f} "
            f"over {int(finite.sum())} ==="
        )
        for against in variables:
            rows = trend(va, vb, against, n_bins)
            if rows:
                out.append("")
                out.append(format_trend(against, rows, value))
    return "\n".join(out)
