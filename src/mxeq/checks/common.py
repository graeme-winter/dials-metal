"""Helpers shared by the boundary checks."""

from __future__ import annotations

import numpy as np

from ..expt import ExperimentList
from ..refl import ReflectionTable
from ..report import Report, Section
from ..stats import describe, resolution_bins

#: DIALS reflection flags, as far as the checks need them.
FLAGS = {
    "strong": 1 << 5,
    "used_in_refinement": 1 << 3,
    "indexed": 1 << 2,
    "integrated_sum": 1 << 11,
    "integrated_prf": 1 << 12,
    "bad_shoebox": 1 << 16,
}


def position_column(table: ReflectionTable) -> str:
    """The best available observed position column.

    Preferring the observed over the calculated is deliberate: two pipelines
    that disagree about the model still agree about where the photons landed,
    so the observed position is the one thing that can be used as a join key
    before the models have been compared.
    """
    for name in ("xyzobs.px.value", "xyzobs.mm.value", "xyzcal.px", "xyzcal.mm"):
        if name in table:
            return name
    raise KeyError(
        "no position column: looked for xyzobs.px.value, xyzobs.mm.value, "
        f"xyzcal.px, xyzcal.mm; have {', '.join(sorted(table.columns))}"
    )


def common_columns(
    a: ReflectionTable, b: ReflectionTable, section: Section
) -> set[str]:
    """Record which columns each table has that the other does not."""
    only_a = sorted(set(a.columns) - set(b.columns))
    only_b = sorted(set(b.columns) - set(a.columns))
    shared = set(a.columns) & set(b.columns)
    section.scalar("n_columns_shared", len(shared))
    if only_a:
        section.note(f"  only in A                    {', '.join(only_a)}")
        section.values["columns_only_a"] = only_a
    if only_b:
        section.note(f"  only in B                    {', '.join(only_b)}")
        section.values["columns_only_b"] = only_b
    return shared


def report_axis_offsets(
    section: Section,
    xyz_a: np.ndarray,
    xyz_b: np.ndarray,
    labels: tuple[str, str, str] = ("x", "y", "z"),
    prefix: str = "offset",
) -> None:
    """Per-axis A minus B, separately, because they fail separately.

    A systematic offset in z alone is a scan or image-numbering problem; one in
    x and y is geometry; one in all three at the sub-pixel level is arithmetic.
    A single scalar RMSD cannot tell those apart.
    """
    delta = np.asarray(xyz_a, dtype=float) - np.asarray(xyz_b, dtype=float)
    for k, label in enumerate(labels):
        section.summary(f"{prefix} {label}", describe(delta[:, k]))


def report_column_agreement(
    section: Section,
    a: ReflectionTable,
    b: ReflectionTable,
    index_a: np.ndarray,
    index_b: np.ndarray,
    columns: list[str] | None = None,
) -> None:
    """Per-column agreement over matched pairs, including bitwise equality.

    Identity is achievable at the early boundaries -- a transcribed threshold
    kernel really can produce the same integers -- so when it has been
    achieved the check should say so outright rather than print a relative
    difference of zero and leave it to be inferred.

    Row order is reported too, because it is the thing that makes two
    equivalent files fail a byte compare.  Two tables holding the same spots in
    a different order differ under ``cmp`` and agree under every metric here,
    and without the reordering count that looks like a contradiction.
    """
    shared = sorted(set(a.columns) & set(b.columns)) if columns is None else columns
    reordered = int(np.sum(np.asarray(index_a) != np.asarray(index_b)))
    section.scalar("n_rows_in_a_different_position", reordered)
    if reordered:
        section.note(
            "  The two tables hold matched rows in different positions, so they"
        )
        section.note("  will differ under a byte compare whatever the values say.")

    rows = []
    for name in shared:
        va, vb = a[name][index_a], b[name][index_b]
        bitwise = va.tobytes() == vb.tobytes()
        if va.dtype.kind == "f":
            worst = float(np.nanmax(np.abs(va - vb))) if va.size else 0.0
            rows.append([name, "yes" if bitwise else "no", f"{worst:.3e}"])
        else:
            same = bool(np.array_equal(va, vb))
            rows.append([name, "yes" if bitwise else "no", "-" if same else "differs"])
    section.table(["column", "bitwise", "max |A-B|"], rows, name="column_agreement")


def shell_table(
    section: Section,
    d: np.ndarray,
    quantities: dict[str, np.ndarray],
    n_bins: int = 10,
    name: str = "shells",
) -> None:
    """A per-resolution-shell table of medians, plus the count in each shell.

    Resolution shells are not decoration.  A pipeline difference that is
    invisible overall and confined to the outer shell is the normal way this
    goes wrong, because that is where the weak reflections are and where any
    difference in background or profile treatment shows up first.
    """
    which, edges = resolution_bins(d, n_bins)
    section.note(f"  by resolution shell, {n_bins} equal-volume bins, medians:")
    header = ["d_max", "d_min", "n", *quantities]
    rows = []
    for b in range(n_bins):
        pick = which == b
        n = int(pick.sum())
        if n == 0:
            continue
        row = [f"{edges[b]:.2f}", f"{edges[b + 1]:.2f}", str(n)]
        for values in quantities.values():
            v = np.asarray(values, dtype=float)[pick]
            v = v[np.isfinite(v)]
            row.append(f"{np.median(v): .4g}" if v.size else "-")
        rows.append(row)
    section.table(header, rows, name=name)


def crystal_setting(experiments: ExperimentList | None) -> np.ndarray | None:
    if experiments is None or len(experiments) == 0:
        return None
    crystal = experiments[0].crystal
    return None if crystal is None else crystal.setting


def note_row_counts(report: Report, a: ReflectionTable, b: ReflectionTable) -> Section:
    section = report.section("rows and columns")
    section.scalar("nrows_a", a.nrows)
    section.scalar("nrows_b", b.nrows)
    if a.nrows and b.nrows:
        section.scalar("nrows_ratio", b.nrows / a.nrows)
    if a.opaque or b.opaque:
        undecoded = sorted(set(a.opaque) | set(b.opaque))
        section.note(f"  not decoded                  {', '.join(undecoded)}")
    common_columns(a, b, section)
    return section
