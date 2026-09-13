"""The spot-finding boundary: two ``strong.refl`` with no key between them.

This is the only boundary where the two tables need not describe the same set
of objects at all, so the interesting output is not how well the matched spots
agree -- it is how many did not match, and what those look like.  A pipeline
that finds 95% of the same spots with sub-pixel centroids and drops 5% of the
weak ones is a different thing from one that finds all of them slightly wrong,
and only the unmatched population separates the two.
"""

from __future__ import annotations

import numpy as np

from .. import match
from ..refl import ReflectionTable
from ..report import Report
from ..stats import correlation, describe, relative_difference
from .common import (
    note_row_counts,
    position_column,
    report_axis_offsets,
    report_column_agreement,
)


def check(
    a: ReflectionTable,
    b: ReflectionTable,
    file_a: str = "A",
    file_b: str = "B",
    radius: float = 2.0,
    z_scale: float = 1.0,
) -> Report:
    """Compare two strong-spot tables.

    ``radius`` is in pixels and ``z_scale`` converts images to the same units,
    so ``z_scale=1`` treats one image as one pixel.  The default is generous:
    the point of the first pass is to find the population that matches at all,
    and a tight radius reports disagreement as absence.
    """
    report = Report("strong", file_a, file_b)
    note_row_counts(report, a, b)

    name_a, name_b = position_column(a), position_column(b)
    if name_a != name_b:
        report.warn(f"matching {name_a} in A against {name_b} in B")
    xyz_a, xyz_b = a[name_a], b[name_b]

    m = match.match_positions(
        xyz_a, xyz_b, radius=radius, scale=(1.0, 1.0, 1.0 / z_scale if z_scale else 1.0)
    )

    s = report.section("spatial match")
    s.scalar("radius_px", radius)
    s.scalar("n_matched", m.n_matched)
    s.scalar("n_only_a", m.n_only_a)
    s.scalar("n_only_b", m.n_only_b)
    s.scalar("fraction_matched", m.fraction_matched)
    if m.n_matched == 0:
        report.warn(
            "nothing matched. Either the two tables are of different data, or "
            "the positions are in different units, or z is offset by more than "
            "the radius -- check the z offset summary in a wider-radius run."
        )
        return report
    s.summary("separation", describe(m.distance))

    report_axis_offsets(
        report.section("centroid offset, A minus B"), xyz_a[m.index_a], xyz_b[m.index_b]
    )

    # The unmatched population, described by whatever is available to describe
    # it with.  If these are simply the weakest spots then the two finders
    # agree and differ only in where they cut; if they are as strong as the
    # matched ones then something is wrong with the finding itself.
    u = report.section("unmatched, by strength")
    for table, m_index, label in ((a, m.index_a, "A"), (b, m.index_b, "B")):
        mask = np.ones(table.nrows, dtype=bool)
        mask[m_index] = False
        for column in ("intensity.sum.value", "n_signal"):
            if column not in table:
                continue
            u.summary(f"{label} {column} matched", describe(table[column][m_index]))
            u.summary(f"{label} {column} unmatched", describe(table[column][mask]))

    shared = report.section("shared columns over matched pairs")
    report_column_agreement(shared, a, b, m.index_a, m.index_b)
    for column in ("intensity.sum.value", "intensity.sum.variance", "n_signal"):
        if column not in a or column not in b:
            continue
        va = np.asarray(a[column], dtype=float)[m.index_a]
        vb = np.asarray(b[column], dtype=float)[m.index_b]
        shared.summary(f"{column} rel. diff", describe(relative_difference(va, vb)))
        pearson, spearman = correlation(va, vb)
        shared.scalar(f"{column} cc", pearson)
        shared.scalar(f"{column} cc_spearman", spearman)

    # Spot depth: a reflection swept over a degree lands on several images, and
    # grouping it per image instead of in three dimensions splits it.  That
    # failure shows up here as a depth distribution pinned at one, and nowhere
    # else until indexing quietly does worse.
    for table, label in ((a, "A"), (b, "B")):
        if "bbox" not in table:
            continue
        bbox = table["bbox"]
        depth = bbox[:, 5] - bbox[:, 4]
        d = report.section(f"spot depth in images, {label}")
        d.summary("depth", describe(depth), fmt="{: .3g}")
        d.scalar("fraction_single_image", float(np.mean(depth <= 1)))

    z = report.section("counts along the scan")
    n_blocks = 10
    all_z = np.concatenate([xyz_a[:, 2], xyz_b[:, 2]])
    edges = np.linspace(all_z.min(), all_z.max() + 1e-9, n_blocks + 1)
    rows = []
    for k in range(n_blocks):
        lo, hi = edges[k], edges[k + 1]
        in_a = int(np.sum((xyz_a[:, 2] >= lo) & (xyz_a[:, 2] < hi)))
        in_b = int(np.sum((xyz_b[:, 2] >= lo) & (xyz_b[:, 2] < hi)))
        matched = int(np.sum((xyz_a[m.index_a, 2] >= lo) & (xyz_a[m.index_a, 2] < hi)))
        rows.append([f"{lo:.0f}", f"{hi:.0f}", str(in_a), str(in_b), str(matched)])
    z.table(["z_from", "z_to", "n_A", "n_B", "matched"], rows, name="scan_blocks")

    return report
