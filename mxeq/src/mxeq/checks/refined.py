"""The refinement boundary: two experiment models.

Refined parameters are not comparable pointwise.  The same lattice can be
described in a different basis, the same orientation reached by a different
rotation, and a scan-varying model has a different number of parameters
depending on how many intervals the refinement chose.  So nothing here
compares parameter vectors.  What it compares is the cell, the misorientation
between the two settings, and the geometry -- quantities that mean the same
thing regardless of how either model was parameterised.

The misorientation is taken from the polar decomposition of ``A_b A_a^-1``,
minimised over the candidate lattice operators.  Without that minimisation, two
models related by a reindexing report a misorientation of ninety degrees and
the check is useless on exactly the cases it exists for.
"""

from __future__ import annotations

import math

import numpy as np

from ..expt import Crystal, ExperimentList
from ..reindex import candidate_operators
from ..report import Report
from ..stats import describe

CELL_LABELS = ("a", "b", "c", "alpha", "beta", "gamma")


def rotation_angle(m: np.ndarray) -> float:
    """The rotation angle in degrees of the rotational part of ``m``.

    The rotation is extracted by polar decomposition rather than assumed,
    because two settings with slightly different cells give an ``m`` that is
    not orthogonal and ``arccos`` of its trace would then be meaningless.
    """
    u, _, vt = np.linalg.svd(np.asarray(m, dtype=float))
    r = u @ vt
    if np.linalg.det(r) < 0:
        # A reflection is not a misorientation; flip the smallest singular
        # direction to get the nearest proper rotation.
        u[:, -1] *= -1
        r = u @ vt
    cos = (np.trace(r) - 1.0) / 2.0
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def misorientation(a: np.ndarray, b: np.ndarray) -> tuple[float, np.ndarray]:
    """Smallest misorientation between two A matrices, over lattice operators."""
    best_angle = float("inf")
    best_operator = np.eye(3, dtype=int)
    inverse = np.linalg.inv(np.asarray(a, dtype=float))
    for m in candidate_operators():
        angle = rotation_angle(np.asarray(b, dtype=float) @ m.astype(float) @ inverse)
        if angle < best_angle:
            best_angle, best_operator = angle, m
    return best_angle, best_operator


def _report_cells(report: Report, ca: Crystal, cb: Crystal) -> None:
    s = report.section("unit cell")
    cell_a, cell_b = ca.cell, cb.cell
    rows = []
    for k, label in enumerate(CELL_LABELS):
        delta = cell_b[k] - cell_a[k]
        relative = delta / cell_a[k] if cell_a[k] else float("nan")
        rows.append(
            [
                label,
                f"{cell_a[k]:.4f}",
                f"{cell_b[k]:.4f}",
                f"{delta:+.5f}",
                f"{relative:+.2e}",
            ]
        )
    s.table(["", "A", "B", "B-A", "rel"], rows, name="cell")
    s.scalar("volume_a", ca.volume)
    s.scalar("volume_b", cb.volume)
    s.scalar("volume_relative_difference", (cb.volume - ca.volume) / ca.volume)


def check(
    a: ExperimentList,
    b: ExperimentList,
    file_a: str = "A",
    file_b: str = "B",
) -> Report:
    report = Report("refined", file_a, file_b)

    s = report.section("experiments")
    s.scalar("n_experiments_a", len(a))
    s.scalar("n_experiments_b", len(b))
    if len(a) != len(b):
        report.warn(
            "the two files hold different numbers of experiments; comparing the "
            "first of each only"
        )
    if len(a) == 0 or len(b) == 0:
        return report

    ea, eb = a[0], b[0]

    if ea.crystal is not None and eb.crystal is not None:
        _report_cells(report, ea.crystal, eb.crystal)

        o = report.section("orientation")
        angle, operator = misorientation(ea.crystal.setting, eb.crystal.setting)
        o.scalar("misorientation_deg", angle)
        o.scalar("operator", operator.tolist())
        if not np.array_equal(operator, np.eye(3, dtype=int)):
            o.note("  The minimum was found through a non-identity lattice operator:")
            o.note("  the two models describe the same lattice in different bases.")

        v = report.section("scan-varying model")
        v.scalar("scan_varying_a", ea.crystal.scan_varying)
        v.scalar("scan_varying_b", eb.crystal.scan_varying)
        if ea.crystal.scan_varying and eb.crystal.scan_varying:
            cells_a = ea.crystal.cells_at_scan_points()
            cells_b = eb.crystal.cells_at_scan_points()
            v.scalar("n_scan_points_a", len(cells_a))
            v.scalar("n_scan_points_b", len(cells_b))
            if len(cells_a) == len(cells_b):
                for k, label in enumerate(CELL_LABELS):
                    v.summary(
                        f"scan-point {label} B-A",
                        describe(cells_b[:, k] - cells_a[:, k]),
                    )
                angles = [
                    misorientation(
                        ea.crystal.a_at_scan_points[i],
                        eb.crystal.a_at_scan_points[i],
                    )[0]
                    for i in range(len(cells_a))
                ]
                v.summary("scan-point misorientation deg", describe(angles))
            else:
                v.note("  Different numbers of scan points: the per-point cells are")
                v.note("  not comparable. The drift over the scan is still comparable")
                v.note("  through its range, reported below.")
                for label, cells in (("A", cells_a), ("B", cells_b)):
                    span = cells.max(axis=0) - cells.min(axis=0)
                    v.note(
                        f"  {label} cell range over scan     "
                        + " ".join(f"{x:.4f}" for x in span)
                    )
    else:
        report.warn("one of the experiment lists has no crystal model")

    if ea.detector is not None and eb.detector is not None:
        g = report.section("detector")
        g.scalar("n_panels_a", ea.detector.num_panels)
        g.scalar("n_panels_b", eb.detector.num_panels)
        g.scalar("distance_a_mm", ea.detector.normal_distance())
        g.scalar("distance_b_mm", eb.detector.normal_distance())
        g.scalar(
            "distance_difference_mm",
            eb.detector.normal_distance() - ea.detector.normal_distance(),
        )
        if ea.detector.num_panels == eb.detector.num_panels:
            shift = np.linalg.norm(eb.detector.origins - ea.detector.origins, axis=1)
            g.summary("panel origin shift mm", describe(shift))
            # An axis difference is a rotation of the panel, which a distance
            # and an origin shift together can hide entirely.
            for name, axes_a, axes_b in (
                ("fast", ea.detector.fast_axes, eb.detector.fast_axes),
                ("slow", ea.detector.slow_axes, eb.detector.slow_axes),
            ):
                dots = np.clip(np.sum(axes_a * axes_b, axis=1), -1.0, 1.0)
                g.summary(
                    f"panel {name} axis angle deg",
                    describe(np.degrees(np.arccos(dots))),
                )

    if ea.scan is not None and eb.scan is not None:
        t = report.section("scan")
        t.scalar("image_range_a", list(ea.scan.image_range))
        t.scalar("image_range_b", list(eb.scan.image_range))
        t.scalar("oscillation_a", list(ea.scan.oscillation))
        t.scalar("oscillation_b", list(eb.scan.oscillation))
        t.scalar("batch_offset_a", ea.scan.batch_offset)
        t.scalar("batch_offset_b", eb.scan.batch_offset)
        if ea.scan.image_range != eb.scan.image_range:
            report.warn(
                "the scans cover different image ranges; any z comparison "
                "downstream is offset by that difference and not by an error"
            )

    if ea.beam is not None and eb.beam is not None:
        w = report.section("beam")
        w.scalar("wavelength_a", ea.beam.wavelength)
        w.scalar("wavelength_b", eb.beam.wavelength)
        w.scalar("wavelength_difference", eb.beam.wavelength - ea.beam.wavelength)

    return report
