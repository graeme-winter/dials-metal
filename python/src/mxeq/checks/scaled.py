"""The scaling and merging boundary.

At this point equivalence stops being about individual reflections.  Two
scaling runs can disagree about every inverse scale factor and produce merged
data that is the same to within the noise, because an overall scale and a
smooth function of the scan are not observable in the merged result.  So the
comparison is made on quantities that survive that freedom: the internal
consistency of each dataset (CC-half, Rmeas), and the correlation between the
two merged sets after both have been put into the same asymmetric unit.

The unmerged scale factors are still reported, but as a description of what the
two models did, not as a pass condition.
"""

from __future__ import annotations

import numpy as np

from .. import match
from ..expt import ExperimentList
from ..refl import ReflectionTable
from ..report import Report
from ..stats import (
    cc_half,
    correlation,
    d_spacing,
    describe,
    merge,
    relative_difference,
    resolution_bins,
    to_asu,
)
from .common import crystal_setting, note_row_counts
from .integrated import _frames, _keys


def _intensity(table: ReflectionTable) -> tuple[str, str] | None:
    for value, variance in (
        ("intensity.scale.value", "intensity.scale.variance"),
        ("intensity.prf.value", "intensity.prf.variance"),
        ("intensity.sum.value", "intensity.sum.variance"),
    ):
        if value in table and variance in table:
            return value, variance
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
    report = Report("scaled", file_a, file_b)
    note_row_counts(report, a, b)

    columns_a, columns_b = _intensity(a), _intensity(b)
    if columns_a is None or columns_b is None:
        report.warn("no intensity/variance pair found in one of the tables")
        return report
    if columns_a != columns_b:
        report.warn(f"comparing {columns_a[0]} in A against {columns_b[0]} in B")

    # Unmerged, reflection by reflection.
    m, n_duplicate = match.match_keys(
        _keys(a, operator), _keys(b, None), _frames(a), _frames(b)
    )
    s = report.section("unmerged match")
    s.scalar("n_matched", m.n_matched)
    s.scalar("fraction_matched", m.fraction_matched)
    s.scalar("n_duplicate_keys", n_duplicate)

    if m.n_matched:
        ia, ib = m.index_a, m.index_b
        va = np.asarray(a[columns_a[0]], dtype=float)[ia]
        vb = np.asarray(b[columns_b[0]], dtype=float)[ib]
        s.scalar("cc_unmerged", correlation(va, vb)[0])
        s.summary("relative difference", describe(relative_difference(va, vb)))
        for column in ("inverse_scale_factor", "inverse_scale_factor_variance"):
            if column in a and column in b:
                fa = np.asarray(a[column], dtype=float)[ia]
                fb = np.asarray(b[column], dtype=float)[ib]
                s.summary(f"{column} A", describe(fa))
                s.summary(
                    f"{column} ratio B/A", describe(np.where(fa != 0, fb / fa, np.nan))
                )
                s.note("  A constant ratio is an overall scale and is not observable;")
                s.note("  its spread is what says the two models differ.")

    # Merged, which is what anyone downstream actually receives.
    hall = None
    if experiments is not None and len(experiments) and experiments[0].crystal:
        hall = experiments[0].crystal.hall
    if hall is None:
        report.warn(
            "no space group available: pass a scaled .expt to get merging "
            "statistics. Without it only the unmerged comparison above is possible."
        )
        return report

    merged = {}
    for label, table, (value, variance) in (
        ("A", a, columns_a),
        ("B", b, columns_b),
    ):
        hkl = table["miller_index"]
        if label == "A" and operator is not None:
            hkl = (hkl.astype(np.int64) @ operator).astype(np.int64)
        asu = to_asu(hkl, hall)
        if asu is None:
            report.warn(f"could not map {label} into the asymmetric unit")
            return report
        intensity = np.asarray(table[value], dtype=float)
        sigma = np.sqrt(np.maximum(np.asarray(table[variance], dtype=float), 0.0))
        keep = np.isfinite(intensity) & (sigma > 0)
        merged[label] = merge(asu[keep], intensity[keep], sigma[keep])

    q = report.section("merging statistics, each dataset alone")
    rows = []
    for label in ("A", "B"):
        mg = merged[label]
        rows.append(
            [
                label,
                str(len(mg.hkl)),
                f"{float(mg.multiplicity.mean()):.2f}",
                f"{cc_half(mg):.4f}",
                f"{mg.r_meas:.4f}",
            ]
        )
    q.table(["", "unique", "mult", "CC-half", "Rmeas"], rows, name="merging")

    # And the comparison between them, which is the equivalence question.
    common, index_a, index_b = np.intersect1d(
        merged["A"]
        .hkl.view([("h", np.int64), ("k", np.int64), ("l", np.int64)])
        .ravel(),
        merged["B"]
        .hkl.view([("h", np.int64), ("k", np.int64), ("l", np.int64)])
        .ravel(),
        return_indices=True,
    )
    c = report.section("merged A against merged B")
    c.scalar("n_common_unique", len(common))
    c.scalar("n_only_a", len(merged["A"].hkl) - len(common))
    c.scalar("n_only_b", len(merged["B"].hkl) - len(common))
    if len(common) < 3:
        report.warn("too few reflections in common to correlate")
        return report

    ma = merged["A"].intensity[index_a]
    mb = merged["B"].intensity[index_b]
    pearson, spearman = correlation(ma, mb)
    c.scalar("cc_merged", pearson)
    c.scalar("cc_merged_spearman", spearman)
    c.summary("relative difference", describe(relative_difference(ma, mb)))
    c.note("  Compare cc_merged against the two CC-half values above: two datasets")
    c.note("  cannot correlate with each other better than their own internal")
    c.note("  consistency allows, so a cc_merged close to those is agreement to")
    c.note("  within the noise, and that is the result being looked for.")

    setting = crystal_setting(experiments)
    if setting is not None:
        d = d_spacing(merged["A"].hkl[index_a], setting)
        which, edges = resolution_bins(d, n_bins)
        rows = []
        for bin_index in range(n_bins):
            pick = which == bin_index
            if pick.sum() < 3:
                continue
            rows.append(
                [
                    f"{edges[bin_index]:.2f}",
                    f"{edges[bin_index + 1]:.2f}",
                    str(int(pick.sum())),
                    f"{correlation(ma[pick], mb[pick])[0]:.4f}",
                    f"{np.median(relative_difference(ma[pick], mb[pick])): .3g}",
                ]
            )
        c.table(["d_max", "d_min", "n", "cc", "med_rel"], rows, name="merged_shells")

    return report
