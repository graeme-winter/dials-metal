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


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Correlation of the RANKS.

    Pearson is dominated by the largest values, so a handful of gross outliers
    in a bin of five thousand drags it down while the other four thousand nine
    hundred agree perfectly.  Spearman cannot be moved that way: an outlier is
    one rank whatever its size.  Where the two disagree, the disagreement is
    the useful reading -- Pearson low and Spearman high means a few disasters,
    both low means real scatter.
    """
    if a.size < 3:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    denominator = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    if denominator <= 0:
        return float("nan")
    return float((ra * rb).sum() / denominator)


def robust_spread(ratio: np.ndarray) -> float:
    """Median absolute deviation of the ratio, scaled to a standard deviation.

    The spread the bulk of the data has, as against the one a few outliers
    impose on it.
    """
    if ratio.size == 0:
        return float("nan")
    middle = np.median(ratio)
    return float(1.4826 * np.median(np.abs(ratio - middle)))


@dataclass
class Row:
    low: float
    high: float
    count: int
    median_ratio: float
    correlation: float
    rank_correlation: float
    spread: float
    outliers: int
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
        spread = robust_spread(ratio)
        middle = float(np.median(ratio))
        # Beyond five robust standard deviations of the bulk, or five per cent
        # away from it, whichever is the wider.
        #
        # The floor matters. Where the bulk agrees exactly the median absolute
        # deviation is zero, and five times zero excludes nothing however wrong
        # a value is -- so a bin of a thousand with twenty catastrophes in it
        # reported none. A disagreement of five per cent is worth counting
        # whatever the rest of the bin is doing.
        limit = max(5.0 * spread, 0.05 * abs(middle))
        outliers = (
            int(np.sum(np.abs(ratio - middle) > limit)) if limit > 0 else 0
        )
        rows.append(
            Row(
                low=float(low),
                high=float(high),
                count=int(inside.sum()),
                median_ratio=middle,
                correlation=float(r),
                rank_correlation=spearman(a[inside], b[inside]),
                spread=spread,
                outliers=outliers,
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
        f"  {'from':>10} {'to':>10} {'n':>7} {'ours':>9} {'theirs':>9} "
        f"{'ratio':>7} {'spread':>7} {'corr':>7} {'rho':>7} {'out':>5}"
    )
    for r in rows:
        lines.append(
            f"  {r.low:10.4g} {r.high:10.4g} {r.count:7d} "
            f"{r.median_a:9.4g} {r.median_b:9.4g} "
            f"{r.median_ratio:7.4f} {r.spread:7.4f} "
            f"{r.correlation:7.4f} {r.rank_correlation:7.4f} {r.outliers:5d}"
        )
    return "\n".join(lines)


def worst_offenders(
    a, ia: np.ndarray, b, ib: np.ndarray, value: str, n: int
) -> str:
    """The reflections that disagree most, so they can be looked at.

    A count of outliers says how many there are.  Only their Miller indices and
    positions say what they have in common, and whether they are a class of
    reflection or scattered through the data.
    """
    va = a.columns[value].ravel()[ia]
    vb = b.columns[value].ravel()[ib]
    ok = np.isfinite(va) & np.isfinite(vb) & (np.abs(vb) > 0)
    ratio = np.full(va.shape, np.nan)
    ratio[ok] = va[ok] / vb[ok]
    middle = np.nanmedian(ratio)
    order = np.argsort(-np.abs(ratio - middle))
    lines = [f"the {n} worst {value}, against a median ratio of {middle:.4f}"]
    lines.append(
        f"  {'h':>5} {'k':>5} {'l':>5} {'fast':>8} {'slow':>8} {'frame':>8} "
        f"{'ours':>11} {'theirs':>11} {'ratio':>9}"
    )
    hkl = a.columns["miller_index"][ia]
    px = a.columns["xyzcal.px"][ia]
    shown = 0
    for i in order:
        if not np.isfinite(ratio[i]):
            continue
        lines.append(
            f"  {hkl[i,0]:5d} {hkl[i,1]:5d} {hkl[i,2]:5d} "
            f"{px[i,0]:8.1f} {px[i,1]:8.1f} {px[i,2]:8.1f} "
            f"{va[i]:11.4g} {vb[i]:11.4g} {ratio[i]:9.3f}"
        )
        shown += 1
        if shown >= n:
            break
    return "\n".join(lines)


def compare(
    a: refl.ReflectionTable,
    b: refl.ReflectionTable,
    values: list[str],
    n_bins: int = 10,
    radius: float = 0.5,
    worst: int = 0,
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
        if worst:
            out.append("")
            out.append(worst_offenders(a, ia, b, ib, value, worst))
    return "\n".join(out)
