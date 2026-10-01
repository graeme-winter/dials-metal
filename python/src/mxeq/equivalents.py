"""Where integration is biased, against the data itself: each observation
against its symmetry equivalents.

Equivalent observations of a reflection are measured at different places on
the detector, at different rotations, partly or fully recorded. If one kind of
observation is integrated wrongly -- partials high, reflections crossing a
module gap low, the first images of a scan high -- it disagrees, on average,
with the equivalents that are not of that kind. That needs no second program
and no reference data set, only the scaled data.

The reference for an observation is the weighted mean of its CLEAN equivalents,
leaving itself out: fully recorded (partiality at least 0.99), nothing of the
profile on masked pixels (profile.measured at least 0.999), and away from both
ends of the scan. Measuring against every equivalent instead would put the bias
into the reference -- with a fifth of the equivalents partials 10 per cent
high, the reference is 2 per cent high, and partials would read +8 and fulls
-2, neither of them true. Against clean equivalents a bias reads as itself.

Only observations both strong and well determined by their reference are
counted -- I/sigma and the reference's I/sigma both at least --min-i-sigma -- so
that the mean of the relative difference is not counting noise. Intensities are
corrected as scaling corrects them -- times lp, over qe and the partiality --
and divided by inverse_scale_factor where the table has it; rows scaling
rejected as outliers are left out.

    mxeq equivalents scaled.expt scaled.refl

Each table counts, by default, only the observations clean in every OTHER
respect, so that it measures its own cause: an observation on the first image
is usually a partial too. Read the median and the ratio of the sums, not the
mean, which noise pulls up most where observations are weakest. Read the table
of each explanatory variable down its bins: a bias that is a
property of the variable steps with it, and the clean bins sit near zero.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import gemmi
import numpy as np

from . import expt as expt_module
from . import refl

SUM = 1 << 8
PRF = 1 << 9
FOREGROUND_BAD = 1 << 14
OUTLIER_IN_SCALING = 1 << 23
EXCLUDED_FOR_SCALING = 1 << 24
SCALED = 1 << 26

CLEAN_PARTIALITY = 0.99
CLEAN_MEASURED = 0.999
CLEAN_EDGE_IMAGES = 5.0


def _encode(hkl: np.ndarray) -> np.ndarray:
    """A Miller index as one integer, 11 bits an index offset by 1024."""
    h = hkl.astype(np.int64) + 1024
    return (h[:, 0] << 22) | (h[:, 1] << 11) | h[:, 2]


def asu_keys(hkl: np.ndarray, hall: str, anomalous: bool = False) -> np.ndarray:
    """One integer per row naming its symmetry-unique reflection: the largest
    encoding of any index equivalent under the space group's rotations -- and
    their Friedel mates, unless anomalous. Vectorised over the rows, in chunks;
    which representative is chosen does not matter, only that equivalents share
    it."""
    ops = gemmi.symops_from_hall(hall)
    rotations = []
    for op in ops.sym_ops:
        r = np.array(op.rot, dtype=np.int64) // op.DEN
        rotations.append(r)
    hkl = np.asarray(hkl, dtype=np.int64).reshape(-1, 3)
    out = np.empty(len(hkl), dtype=np.int64)
    step = 200000
    for start in range(0, len(hkl), step):
        part = hkl[start : start + step]
        best = np.full(len(part), np.iinfo(np.int64).min, dtype=np.int64)
        for r in rotations:
            image = part @ r  # a reciprocal-space index is a row: h' = h R
            best = np.maximum(best, _encode(image))
            if not anomalous:
                best = np.maximum(best, _encode(-image))
        out[start : start + step] = best
    return out


@dataclass
class Explanatory:
    """One thing to bin the bias against: bin edges for a quantity, or None
    for a quantity that is already a label."""

    name: str
    values: np.ndarray
    edges: list[float] | None = None
    labels: list[str] | None = None
    #: Which observations this table counts: those clean in every OTHER
    #: respect, so that the table measures its own cause -- an observation on
    #: the first image is usually a partial too. None counts them all.
    condition: np.ndarray | None = None


@dataclass
class Bin:
    label: str
    n: int
    mean: float  # mean of I / reference - 1
    sem: float  # its standard error
    median: float
    sums: float  # sum I / sum reference - 1


@dataclass
class Table:
    explanatory: str
    bins: list[Bin] = field(default_factory=list)


@dataclass
class Result:
    intensity: str
    counted: int
    groups_with_a_reference: int
    overall: Bin
    tables: list[Table] = field(default_factory=list)


def references(
    keys: np.ndarray, intensity: np.ndarray, variance: np.ndarray, clean: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Each row's reference: the weighted mean of the clean rows sharing its
    key, itself left out, and that mean's standard deviation; NaN where no
    other clean row shares the key."""
    _, inverse = np.unique(keys, return_inverse=True)
    w = np.where(variance > 0, 1.0 / np.where(variance > 0, variance, 1.0), 0.0)
    wc = w * clean
    sw = np.bincount(inverse, weights=wc)
    swi = np.bincount(inverse, weights=wc * intensity)
    own_w = wc
    loo_w = sw[inverse] - own_w
    loo_wi = swi[inverse] - own_w * intensity
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(loo_w > 0, loo_wi / np.where(loo_w > 0, loo_w, 1.0), np.nan)
        sigma = np.where(
            loo_w > 0, 1.0 / np.sqrt(np.where(loo_w > 0, loo_w, 1.0)), np.nan
        )
    return mean, sigma


def _summary(
    label: str, rel: np.ndarray, intensity: np.ndarray, reference: np.ndarray
) -> Bin:
    n = len(rel)
    if n == 0:
        return Bin(label, 0, np.nan, np.nan, np.nan, np.nan)
    sem = float(np.std(rel, ddof=1) / np.sqrt(n)) if n > 1 else np.nan
    return Bin(
        label,
        n,
        float(np.mean(rel)),
        sem,
        float(np.median(rel)),
        float(np.sum(intensity) / np.sum(reference) - 1.0),
    )


def analyse(
    keys: np.ndarray,
    intensity: np.ndarray,
    variance: np.ndarray,
    clean: np.ndarray,
    explanatories: list[Explanatory],
    min_i_sigma: float = 5.0,
    name: str = "",
) -> Result:
    """The bias of every observation against its clean equivalents, overall
    and binned by each explanatory variable."""
    reference, ref_sigma = references(keys, intensity, variance, clean)
    sigma = np.sqrt(np.where(variance > 0, variance, np.nan))
    with np.errstate(invalid="ignore", divide="ignore"):
        counted = (
            np.isfinite(reference)
            & (reference > 0)
            & (intensity / sigma >= min_i_sigma)
            & (reference / ref_sigma >= min_i_sigma)
        )
        rel = intensity / reference - 1.0
    result = Result(
        name,
        int(counted.sum()),
        int(np.unique(keys[np.isfinite(reference)]).size),
        _summary("all", rel[counted], intensity[counted], reference[counted]),
    )
    for e in explanatories:
        table = Table(e.name)
        v = e.values
        base = counted if e.condition is None else counted & e.condition
        if e.labels is not None and e.edges is None:
            for lab in e.labels:
                sel = base & (v == lab)
                table.bins.append(
                    _summary(str(lab), rel[sel], intensity[sel], reference[sel])
                )
        else:
            edges = e.edges
            for k in range(len(edges) - 1):
                lo, hi = edges[k], edges[k + 1]
                sel = base & (v >= lo) & (v < hi)
                label = e.labels[k] if e.labels else f"[{lo:g}, {hi:g})"
                table.bins.append(
                    _summary(label, rel[sel], intensity[sel], reference[sel])
                )
        result.tables.append(table)
    return result


def from_files(
    expt_path: str,
    refl_path: str,
    kinds=("prf", "sum"),
    min_i_sigma: float = 5.0,
    anomalous: bool = False,
    others_clean: bool = True,
) -> list[Result]:
    """The analysis for a scaled (or integrated) .expt and .refl."""
    experiments = expt_module.load(expt_path)
    if not experiments.experiments or experiments[0].crystal is None:
        raise ValueError(f"{expt_path} has no crystal")
    e0 = experiments[0]
    hall = e0.crystal.hall
    if not hall:
        raise ValueError(f"{expt_path}'s crystal has no space group")
    t = refl.load(refl_path)
    c = t.columns

    def column(name, default=None):
        if name in c:
            return np.asarray(c[name])
        if default is None:
            raise ValueError(f"{refl_path} has no column {name}")
        return np.full(t.nrows, default)

    flags = column("flags").astype(np.int64)
    hkl = column("miller_index").reshape(-1, 3)
    scale = column("inverse_scale_factor", 1.0).astype(float)
    # The integrated columns are raw, as DIALS's are: the Lorentz and
    # polarisation correction, the detector's efficiency and the partiality are
    # applied as scaling applies them -- times lp, over qe, over the partiality
    # -- or equivalents at different geometry differ by the Lorentz factor alone.
    lp = column("lp", 1.0).astype(float)
    qe = column("qe", 1.0).astype(float)
    qe = np.where(qe > 0, qe, 1.0)
    partiality = column("partiality", 1.0).astype(float)
    measured = column("profile.measured", 1.0).astype(float)
    z = (
        column("xyzcal.px").reshape(-1, 3)[:, 2]
        if "xyzcal.px" in c
        else np.zeros(t.nrows)
    )
    d = column("d", np.nan).astype(float)
    images = (
        float(e0.scan.num_images) if e0.scan is not None else float(np.nanmax(z) + 1)
    )
    from_start = z
    from_end = images - z

    keep = (flags & (OUTLIER_IN_SCALING | EXCLUDED_FOR_SCALING)) == 0
    if np.any(flags & SCALED):
        keep &= (flags & SCALED) != 0
    keys = asu_keys(hkl, hall, anomalous)
    full = partiality >= CLEAN_PARTIALITY
    unmasked = (measured >= CLEAN_MEASURED) & ((flags & FOREGROUND_BAD) == 0)
    interior = (from_start >= CLEAN_EDGE_IMAGES) & (from_end >= CLEAN_EDGE_IMAGES)
    clean_base = full & unmasked & interior
    fg_bad = np.where((flags & FOREGROUND_BAD) != 0, "yes", "no")
    shells = np.nanquantile(d, np.linspace(0, 1, 7)) if np.isfinite(d).any() else [0, 1]
    shells = sorted(set(float(x) for x in shells))
    shells[-1] = shells[-1] * 1.0001 + 1e-9

    results = []
    for kind in kinds:
        bit = PRF if kind == "prf" else SUM
        value = column(f"intensity.{kind}.value").astype(float)
        var = column(f"intensity.{kind}.variance").astype(float)
        rows = keep & ((flags & bit) != 0) & (scale > 0) & (partiality > 0)
        factor = lp[rows] / qe[rows] / partiality[rows] / scale[rows]
        intensity = value[rows] * factor
        variance = var[rows] * factor**2
        isig = intensity / np.sqrt(np.where(variance > 0, variance, np.nan))
        others = others_clean
        explanatories = [
            Explanatory(
                "partiality",
                partiality[rows],
                [0.0, 0.5, 0.8, 0.9, 0.99, 1.0001],
                condition=(unmasked & interior)[rows] if others else None,
            ),
            Explanatory(
                "profile.measured: the profile on valid pixels",
                measured[rows],
                [0.0, 0.6, 0.8, 0.9, 0.95, 0.999, 1.0001],
                condition=(full & interior)[rows] if others else None,
            ),
            Explanatory(
                "the foreground on a masked pixel",
                fg_bad[rows],
                None,
                ["no", "yes"],
                condition=(full & interior)[rows] if others else None,
            ),
            Explanatory(
                "images from the scan's start",
                from_start[rows],
                [0, 1, 2, 5, 10, np.inf],
                ["first", "second", "3rd-5th", "6th-10th", "later"],
                condition=(full & unmasked)[rows] if others else None,
            ),
            Explanatory(
                "images before the scan's end",
                from_end[rows],
                [0, 1, 2, 5, 10, np.inf],
                ["last", "last but one", "3rd-5th", "6th-10th", "earlier"],
                condition=(full & unmasked)[rows] if others else None,
            ),
            Explanatory(
                "resolution d (A), a control",
                d[rows],
                shells,
                condition=clean_base[rows] if others else None,
            ),
            Explanatory(
                "I/sigma, a control",
                isig,
                [0, 10, 20, 50, 100, np.inf],
                condition=clean_base[rows] if others else None,
            ),
        ]
        results.append(
            analyse(
                keys[rows],
                intensity,
                variance,
                clean_base[rows],
                explanatories,
                min_i_sigma,
                kind,
            )
        )
    return results


def report(results: list[Result], min_i_sigma: float, others_clean: bool = True) -> str:
    lines = []
    if others_clean:
        lines.append(
            "Each table counts the observations clean in every other respect -- fully"
        )
        lines.append(
            "recorded, nothing masked, away from the scan's ends -- so that it measures its"
        )
        lines.append("own cause; --all counts every observation in every table.")
    lines.append(
        "The median and the ratio of the sums are the measures to read: the mean of a"
    )
    lines.append(
        "ratio is pulled up by noise, most where the observations are weakest."
    )
    lines.append("")
    for r in results:
        name = {"prf": "profile fitted", "sum": "summed"}.get(r.intensity, r.intensity)
        lines.append(
            f"{name} intensities: {r.counted} observations against the clean "
            f"equivalents of {r.groups_with_a_reference} reflections, I/sigma >= {min_i_sigma:g}"
        )
        o = r.overall
        lines.append(
            f"  overall: mean {100 * o.mean:+.2f} % +- {100 * o.sem:.2f}, median "
            f"{100 * o.median:+.2f} %, sums {100 * o.sums:+.2f} %"
        )
        for table in r.tables:
            lines.append("")
            lines.append(f"  {table.explanatory}")
            lines.append(
                f"    {'bin':<14} {'n':>8} {'mean %':>9} {'+-':>6} {'median %':>9} {'sums %':>8}"
            )
            for b in table.bins:
                if b.n == 0:
                    lines.append(f"    {b.label:<14} {0:>8}")
                    continue
                sem = f"{100 * b.sem:6.2f}" if np.isfinite(b.sem) else f"{'':>6}"
                lines.append(
                    f"    {b.label:<14} {b.n:>8} {100 * b.mean:+9.2f} {sem} "
                    f"{100 * b.median:+9.2f} {100 * b.sums:+8.2f}"
                )
        lines.append("")
    return "\n".join(lines)
