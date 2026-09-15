"""The indexing boundary: same spots, possibly a different basis.

The order here matters and is the whole content of this check.  Match
spatially, *then* find the reindexing operator over the matched pairs, *then*
compare Miller indices.  Compare the indices first and a perfectly good pair of
solutions related by a basis change reports near-total disagreement.
"""

from __future__ import annotations

import numpy as np

from .. import match, reindex
from ..expt import ExperimentList
from ..refl import ReflectionTable
from ..report import Report
from ..stats import describe
from .common import note_row_counts, position_column, report_axis_offsets


def _indexed_mask(table: ReflectionTable) -> np.ndarray:
    if "miller_index" not in table:
        return np.zeros(table.nrows, dtype=bool)
    return np.any(table["miller_index"] != 0, axis=1)


def check(
    a: ReflectionTable,
    b: ReflectionTable,
    file_a: str = "A",
    file_b: str = "B",
    experiments_a: ExperimentList | None = None,
    radius: float = 2.0,
) -> Report:
    report = Report("indexed", file_a, file_b)
    note_row_counts(report, a, b)

    if "miller_index" not in a or "miller_index" not in b:
        report.warn("one of the tables has no miller_index column")
        return report

    counts = report.section("how much was indexed")
    mask_a, mask_b = _indexed_mask(a), _indexed_mask(b)
    counts.scalar("n_indexed_a", int(mask_a.sum()))
    counts.scalar("n_indexed_b", int(mask_b.sum()))
    counts.scalar("fraction_indexed_a", float(mask_a.mean()) if a.nrows else 0.0)
    counts.scalar("fraction_indexed_b", float(mask_b.mean()) if b.nrows else 0.0)

    name_a, name_b = position_column(a), position_column(b)
    xyz_a, xyz_b = a[name_a], b[name_b]
    m = match.match_positions(xyz_a, xyz_b, radius=radius)

    s = report.section("spatial match")
    s.scalar("n_matched", m.n_matched)
    s.scalar("fraction_matched", m.fraction_matched)
    if m.n_matched == 0:
        report.warn("no spots matched spatially; nothing further can be compared")
        return report
    s.summary("separation", describe(m.distance))

    hkl_a = a["miller_index"][m.index_a]
    hkl_b = b["miller_index"][m.index_b]

    crystal = experiments_a[0].crystal if experiments_a and len(experiments_a) else None
    found = reindex.find_operator(
        hkl_a, hkl_b, real_space_a=crystal.a if crystal is not None else None
    )

    r = report.section("reindexing")
    r.scalar("n_pairs_both_indexed", found.n_tested)
    r.scalar("operator", found.operator.tolist())
    r.note(
        "  operator rows                "
        + " | ".join(",".join(f"{v:>2d}" for v in row) for row in found.operator)
    )
    r.scalar("is_identity", found.is_identity)
    r.scalar("agreement_best", found.fraction)
    r.scalar("agreement_identity", found.fraction_identity)
    if found.metric_compatible is not None:
        r.scalar("is_a_lattice_symmetry", found.metric_compatible)
        if not found.metric_compatible and found.fraction < 0.9:
            report.warn(
                "the best operator neither preserves the cell metric nor "
                "explains most pairs; it is unlikely to be a change of basis"
            )
        elif not found.metric_compatible:
            r.note("  Not a symmetry of the cell metric, which is normal: two")
            r.note("  different reduced cells of one lattice are related by a")
            r.note("  general unimodular matrix, not by a point group operator.")
    if not found.is_identity:
        r.note("  B is indexed in a different basis from A. Everything downstream")
        r.note("  must be compared through this operator, not directly.")
    if found.fraction < 0.5:
        report.warn(
            f"the best operator only accounts for {found.fraction:.1%} of pairs: "
            "the two solutions are not the same lattice"
        )

    # Disagreement after reindexing, which is the number that matters.
    both = np.any(hkl_a != 0, axis=1) & np.any(hkl_b != 0, axis=1)
    transformed = found.apply(hkl_a[both])
    agree = np.all(transformed == hkl_b[both], axis=1)
    d = report.section("index agreement after reindexing")
    d.scalar("n_compared", int(both.sum()))
    d.scalar("n_agree", int(agree.sum()))
    d.scalar("n_disagree", int((~agree).sum()))
    if (~agree).any():
        delta = np.abs(transformed - hkl_b[both])[~agree]
        d.summary("disagreement |dh| max", describe(delta.max(axis=1)), fmt="{: .3g}")
        d.note(
            "  A disagreement of one unit in a single index is a near-miss in "
            "the assignment;"
        )
        d.note("  a large one means the two are indexing different lattices.")

    # One side indexed, the other not: this is where the two differ in how much
    # of the observed data they can account for, and it is separate from
    # disagreeing about what a spot is.
    e = report.section("indexed on one side only")
    only_a = np.any(hkl_a != 0, axis=1) & np.all(hkl_b == 0, axis=1)
    only_b = np.all(hkl_a == 0, axis=1) & np.any(hkl_b != 0, axis=1)
    e.scalar("n_indexed_only_a", int(only_a.sum()))
    e.scalar("n_indexed_only_b", int(only_b.sum()))

    for table, index, label in ((a, m.index_a, "A"), (b, m.index_b, "B")):
        if "xyzcal.px" not in table or "xyzobs.px.value" not in table:
            continue
        # Only where there IS a prediction. An unindexed reflection carries
        # xyzcal.px of exactly zero, which is the corner of the detector, so
        # including those rows reported an rmsd of 755 px on a dataset whose
        # real value is a third of one.
        from ..refl import has_prediction

        predicted = has_prediction(table, index)
        if predicted.sum() < 10:
            continue
        residual = (
            table["xyzcal.px"][index][predicted]
            - table["xyzobs.px.value"][index][predicted]
        )
        # A distribution, not an rms. On real DIALS output 39 of 13072 indexed
        # reflections are predicted more than 5 px from where they were seen,
        # one of them 3908 px away, and an rms over that reports 34 px for a
        # dataset whose median is 0.28. The rest of this package reports
        # distributions for exactly this reason; this line was the exception.
        e.summary(
            f"|xyzcal - xyzobs| px {label}",
            describe(np.linalg.norm(residual[:, :2], axis=1)),
        )

    if "xyzcal.px" in a and "xyzcal.px" in b:
        # Only where BOTH sides actually predicted something. An unindexed
        # reflection carries xyzcal.px of exactly zero, which is the corner of
        # the detector, so comparing every matched row differences real
        # predictions against that corner and reports a median of 1.5e-320 --
        # a denormal, the signature of a column that was never written -- with
        # a mean of -25 px and a first percentile of -1576.
        #
        # The numbers were not wrong so much as meaningless, and they looked
        # like a catastrophic disagreement between two pipelines that in fact
        # agree. Predicted where nothing was predicted is not a disagreement.
        from ..refl import has_prediction

        pa = has_prediction(a, m.index_a)
        pb = has_prediction(b, m.index_b)
        both = pa & pb
        section = report.section("predicted position offset, A minus B")
        section.scalar("n_predicted_by_both", int(both.sum()))
        section.scalar("n_predicted_by_a_only", int((pa & ~pb).sum()))
        section.scalar("n_predicted_by_b_only", int((pb & ~pa).sum()))
        if both.sum() >= 10:
            report_axis_offsets(
                section,
                a["xyzcal.px"][m.index_a][both],
                b["xyzcal.px"][m.index_b][both],
            )
        else:
            report.warn(
                "fewer than ten reflections were predicted by both, so there "
                "is nothing to compare"
            )

    return report
