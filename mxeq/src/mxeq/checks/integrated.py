"""The integration boundary: same reflections, different intensities.

Here there is a key -- (Miller index, entering, experiment) -- and the join is
cheap.  The difficulty moves to the metric.  A percentage difference in
intensity is meaningless when averaged over a dataset whose intensities span
five orders of magnitude, and a single correlation coefficient is dominated by
the hundred strongest reflections.  So everything is reported in resolution
shells, and the headline quantities are the relative difference and the shape
of the pull distribution rather than either on its own.

If the two solutions are in different bases, the caller must supply the
reindexing operator found at the indexing boundary.  Nothing is assumed: with
no operator given, a low match fraction here says the join failed, and the
indexing check says why.
"""

from __future__ import annotations

import numpy as np

from .. import match
from ..expt import ExperimentList
from ..refl import ReflectionTable
from ..report import Report
from ..stats import (
    correlation,
    d_spacing,
    describe,
    pull,
    relative_difference,
)
from .common import crystal_setting, note_row_counts, shell_table

#: Intensity columns compared in pairs, value with its variance.
INTENSITY_PAIRS = (
    ("intensity.sum.value", "intensity.sum.variance"),
    ("intensity.prf.value", "intensity.prf.variance"),
    ("intensity.scale.value", "intensity.scale.variance"),
)

#: Scalar columns worth a distribution of their own.  These are where a
#: disagreement in intensity usually comes from, so reporting them next to it
#: turns "the intensities differ" into a cause.
SCALAR_COLUMNS = (
    "partiality",
    "lp",
    "qe",
    "dqe",
    "background.mean",
    "background.sum.value",
    "profile.correlation",
    "num_pixels.foreground",
    "num_pixels.background",
    "d",
)


def _keys(table: ReflectionTable, operator: np.ndarray | None) -> list[np.ndarray]:
    hkl = table["miller_index"]
    if operator is not None:
        hkl = (hkl.astype(np.int64) @ operator).astype(np.int64)
    columns: list[np.ndarray] = [hkl]
    if "entering" in table:
        columns.append(table["entering"].astype(np.int64))
    if "id" in table:
        columns.append(table["id"].astype(np.int64))
    return columns


def _frames(table: ReflectionTable) -> np.ndarray | None:
    for name in ("xyzcal.px", "xyzobs.px.value"):
        if name in table:
            return table[name][:, 2]
    return None


def check(
    a: ReflectionTable,
    b: ReflectionTable,
    file_a: str = "A",
    file_b: str = "B",
    experiments: ExperimentList | None = None,
    operator: np.ndarray | None = None,
    n_bins: int = 10,
) -> Report:
    report = Report("integrated", file_a, file_b)
    note_row_counts(report, a, b)

    if "miller_index" not in a or "miller_index" not in b:
        report.warn("one of the tables has no miller_index column")
        return report

    m, n_duplicate = match.match_keys(
        _keys(a, operator), _keys(b, None), _frames(a), _frames(b)
    )

    s = report.section("keyed match")
    s.note("  key                          miller_index, entering, id")
    if operator is not None:
        s.scalar("operator", np.asarray(operator).tolist())
    s.scalar("n_matched", m.n_matched)
    s.scalar("n_only_a", m.n_only_a)
    s.scalar("n_only_b", m.n_only_b)
    s.scalar("fraction_matched", m.fraction_matched)
    s.scalar("n_duplicate_keys", n_duplicate)
    if n_duplicate:
        s.note("  Duplicate keys were paired in frame order, not by intensity.")
    if m.n_matched == 0:
        report.warn(
            "nothing matched on the key. If the indexing check found a "
            "non-identity operator, pass it with --operator."
        )
        return report

    if m.fraction_matched < 0.9:
        report.warn(
            f"only {m.fraction_matched:.1%} of reflections matched: the two runs "
            "predicted substantially different reflection lists, which is a "
            "result in its own right and not only a nuisance for this join"
        )

    ia, ib = m.index_a, m.index_b

    # Resolution comes from the experiment model where one is available, and
    # from the table's own d column otherwise. Computing it from the model is
    # preferred because a d column is itself a pipeline output and would then
    # be binning the comparison by one of the things being compared.
    setting = crystal_setting(experiments)
    if setting is not None:
        d = d_spacing(a["miller_index"][ia], setting)
    elif "d" in a:
        d = np.asarray(a["d"], dtype=float)[ia]
        report.warn("binning by A's own d column; pass an .expt for an independent one")
    else:
        d = None
        report.warn("no resolution available: give an .expt to get shell tables")

    for value_column, variance_column in INTENSITY_PAIRS:
        if value_column not in a or value_column not in b:
            continue
        va = np.asarray(a[value_column], dtype=float)[ia]
        vb = np.asarray(b[value_column], dtype=float)[ib]
        section = report.section(value_column)
        pearson, spearman = correlation(va, vb)
        section.scalar("cc", pearson)
        section.scalar("cc_spearman", spearman)
        section.summary("A", describe(va))
        section.summary("B", describe(vb))
        rel = relative_difference(va, vb)
        section.summary("relative difference", describe(rel))

        have_sigma = variance_column in a and variance_column in b
        if have_sigma:
            sa = np.sqrt(np.maximum(np.asarray(a[variance_column], dtype=float)[ia], 0))
            sb = np.sqrt(np.maximum(np.asarray(b[variance_column], dtype=float)[ib], 0))
            p = pull(va, sa, vb, sb)
            section.summary("pull", describe(p))
            section.note(
                "  The pull width is expected to be well below one: these are the"
            )
            section.note(
                "  same photons twice, so the errors are correlated. Read its shape."
            )
            section.summary(
                "sigma relative difference", describe(relative_difference(sa, sb))
            )
            with np.errstate(invalid="ignore", divide="ignore"):
                section.summary(
                    "I/sigma A", describe(np.where(sa > 0, va / sa, np.nan))
                )
                section.summary(
                    "I/sigma B", describe(np.where(sb > 0, vb / sb, np.nan))
                )

        if d is not None:
            quantities = {"cc": np.full(len(ia), np.nan), "rel_diff": rel}
            # A per-shell correlation has to be computed per shell, not
            # median-reduced from a per-reflection quantity, so it is filled in
            # by hand rather than passed through the median machinery.
            from ..stats import resolution_bins

            which, _ = resolution_bins(d, n_bins)
            per_shell_cc = np.full(len(ia), np.nan)
            for bin_index in range(n_bins):
                pick = which == bin_index
                if pick.sum() >= 3:
                    per_shell_cc[pick] = correlation(va[pick], vb[pick])[0]
            quantities["cc"] = per_shell_cc
            if have_sigma:
                quantities["pull"] = p
            shell_table(
                section, d, quantities, n_bins=n_bins, name=f"{value_column}_shells"
            )

    flags = report.section("scalar columns")
    for column in SCALAR_COLUMNS:
        if column not in a or column not in b:
            continue
        va = np.asarray(a[column], dtype=float)[ia]
        vb = np.asarray(b[column], dtype=float)[ib]
        flags.summary(f"{column} A", describe(va))
        flags.summary(f"{column} B-A", describe(vb - va))

    if "xyzcal.px" in a and "xyzcal.px" in b:
        from .common import report_axis_offsets

        # Only where both sides predicted something. A reflection with no
        # prediction carries xyzcal.px of exactly zero, the corner of the
        # detector, and differencing real predictions against that corner
        # reports an offset of hundreds of pixels that means nothing.
        from ..refl import has_prediction

        both = has_prediction(a, ia) & has_prediction(b, ib)
        section = report.section("predicted position offset, A minus B")
        section.scalar("n_predicted_by_both", int(both.sum()))
        if both.sum() >= 10:
            report_axis_offsets(
                section, a["xyzcal.px"][ia][both], b["xyzcal.px"][ib][both]
            )

    return report
