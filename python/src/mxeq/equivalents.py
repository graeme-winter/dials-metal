"""Where integration is biased, against the data itself: each observation
against its symmetry equivalents.

Equivalent observations of a reflection are measured at different places on
the detector, at different rotations, partly or fully recorded. If one kind of
observation is integrated wrongly -- partials high, reflections crossing a
module gap low -- it disagrees, on average, with the equivalents that are not
of that kind. That needs no second program and no reference data set, only the
scaled data.

Intensities are corrected as scaling corrects them -- times lp, over qe and the
partiality, over inverse_scale_factor -- since the integrated columns are raw,
and equivalents at different geometry would otherwise differ by the Lorentz
factor alone. The reference for an observation is the unweighted mean of its
CLEAN equivalents, itself left out: fully recorded, nothing of the profile
masked, away from the scan's ends. Clean, so that a bias does not leak into the
reference; unweighted, because weights from each observation's own variance --
computed from its own counts -- give the equivalents that happened to come out
low the most weight, and read an unbiased data set +0.7 per cent high at I/sigma
10 to 20 (simulated; +0.24 in the median and nothing in the sums unweighted).
`--reference weighted` keeps the weighted mean, to compare.

Each table counts, by default, only the observations clean in every OTHER
respect, so that it measures its own cause: an observation on the first image is
usually a partial too. Read the median and the ratio of the sums; the mean of a
ratio is pulled up by noise. Resolution shells are quantiles of what each table
counts.

From the partials, a suggested sigma_m: a partial's intensity before the
division by its partiality, against its clean reference, is its observed
partiality. Partiality is recomputed as the integrator computes it -- the
fraction of a Gaussian of width sigma_m / |zeta| in phi lying within the box's
frames -- for a trial sigma_m, and the suggestion is the sigma_m at which the
partials' median bias is zero. The sigma_m the integration used is recovered
from each partial's own partiality by inverting the same function.

    mxeq equivalents scaled.expt scaled.refl
    mxeq equivalents scaled.expt scaled.refl --worst 200    # for dials.image_viewer
"""

from __future__ import annotations

from dataclasses import dataclass, field

import gemmi
import numpy as np
from scipy.special import erf

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
    """One thing to bin the bias against: bin edges for a quantity, labels for a
    quantity already labelled, or a number of quantile bins of what is counted."""

    name: str
    values: np.ndarray
    edges: list[float] | None = None
    labels: list[str] | None = None
    #: Which observations this table counts: those clean in every OTHER
    #: respect, so that the table measures its own cause. None counts them all.
    condition: np.ndarray | None = None
    #: Bins as this many quantiles of the values this table counts, so that
    #: none is empty: resolution shells, where the strong are only some of them.
    quantiles: int | None = None


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
    #: Each analysed row's reference and relative difference, and whether it
    #: was counted, for --worst and the sigma_m suggestion.
    reference: np.ndarray | None = None
    rel: np.ndarray | None = None
    counted_rows: np.ndarray | None = None


def references(keys, intensity, variance, clean, method: str = "mean"):
    """Each row's reference -- the mean of the clean rows sharing its key,
    itself left out, unweighted or weighted by 1/variance -- and that mean's
    standard deviation; NaN where no other clean row shares the key."""
    _, inverse = np.unique(keys, return_inverse=True)
    c = clean.astype(float)
    if method == "weighted":
        w = np.where(variance > 0, 1.0 / np.where(variance > 0, variance, 1.0), 0.0) * c
        sw = np.bincount(inverse, weights=w)[inverse] - w
        swi = np.bincount(inverse, weights=w * intensity)[inverse] - w * intensity
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(sw > 0, swi / np.where(sw > 0, sw, 1.0), np.nan)
            sigma = np.where(sw > 0, 1.0 / np.sqrt(np.where(sw > 0, sw, 1.0)), np.nan)
        return mean, sigma
    n = np.bincount(inverse, weights=c)[inverse] - c
    si = np.bincount(inverse, weights=c * intensity)[inverse] - c * intensity
    v = np.where(variance > 0, variance, 0.0) * c
    sv = np.bincount(inverse, weights=v)[inverse] - v
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(n > 0, si / np.where(n > 0, n, 1.0), np.nan)
        sigma = np.where(n > 0, np.sqrt(sv) / np.where(n > 0, n, 1.0), np.nan)
    return mean, sigma


def _summary(label, rel, intensity, reference) -> Bin:
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
    keys,
    intensity,
    variance,
    clean,
    explanatories,
    min_i_sigma: float = 5.0,
    name: str = "",
    method: str = "mean",
) -> Result:
    """The bias of every observation against its clean equivalents, overall
    and binned by each explanatory variable."""
    reference, ref_sigma = references(keys, intensity, variance, clean, method)
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
        reference=reference,
        rel=rel,
        counted_rows=counted,
    )
    for e in explanatories:
        table = Table(e.name)
        v = e.values
        base = counted if e.condition is None else counted & e.condition
        if e.labels is not None and e.edges is None and e.quantiles is None:
            for lab in e.labels:
                sel = base & (v == lab)
                table.bins.append(
                    _summary(str(lab), rel[sel], intensity[sel], reference[sel])
                )
        else:
            edges = e.edges
            if e.quantiles:
                inside = v[base & np.isfinite(v)]
                if inside.size == 0:
                    result.tables.append(table)
                    continue
                edges = sorted(
                    set(
                        float(x)
                        for x in np.quantile(inside, np.linspace(0, 1, e.quantiles + 1))
                    )
                )
                edges[-1] = np.nextafter(edges[-1], np.inf)
            for k in range(len(edges) - 1):
                lo, hi = edges[k], edges[k + 1]
                sel = base & (v >= lo) & (v < hi)
                label = (
                    e.labels[k]
                    if e.labels
                    else (
                        f"[{lo:.4g}, {hi:.4g})" if e.quantiles else f"[{lo:g}, {hi:g})"
                    )
                )
                table.bins.append(
                    _summary(label, rel[sel], intensity[sel], reference[sel])
                )
        result.tables.append(table)
    return result


def partiality(sigma_m_deg, phi, zeta, lo, hi):
    """As the integrator computes it: the fraction of a Gaussian of width
    sigma_m / |zeta| in phi lying between lo and hi (radians)."""
    spread = np.radians(sigma_m_deg) / np.abs(zeta)
    scale = np.sqrt(2.0) * spread
    with np.errstate(invalid="ignore", divide="ignore"):
        f = 0.5 * (erf((hi - phi) / scale) - erf((lo - phi) / scale))
    return np.clip(np.where(spread > 0, f, 1.0), 0.0, 1.0)


def implied_sigma_m(p, phi, zeta, lo, hi):
    """The sigma_m that gives each row its partiality p: bisection, in log
    sigma_m, of the function above, which falls as sigma_m grows."""
    a = np.full(len(p), -8.0)
    b = np.full(len(p), 3.0)
    for _ in range(60):
        m = 0.5 * (a + b)
        too_wide = partiality(np.exp(m), phi, zeta, lo, hi) < p
        b = np.where(too_wide, m, b)
        a = np.where(too_wide, a, m)
    return np.exp(0.5 * (a + b))


@dataclass
class SigmaM:
    used: float  # the sigma_m the integration used, recovered from the partials
    suggested: float | None  # where the partials' median bias is zero
    n: int  # partials it rests on
    bias_used: float  # the partials' median bias as integrated
    note: str = ""


def suggest_sigma_m(counts, reference, p_model, phi, zeta, lo, hi) -> SigmaM:
    """The sigma_m at which the partials' median bias is zero. `counts` is each
    partial's intensity corrected for everything but its partiality, so that
    counts / reference is its observed partiality; `p_model` is the partiality
    the integration used."""
    n = len(counts)
    used = (
        float(np.median(implied_sigma_m(p_model, phi, zeta, lo, hi))) if n else np.nan
    )
    if n == 0:
        return SigmaM(used, None, 0, np.nan, "no partials to judge it by")

    def bias(sigma):
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = counts / (partiality(sigma, phi, zeta, lo, hi) * reference)
        return float(np.median(ratio) - 1.0)

    bias_used = float(np.median(counts / (p_model * reference)) - 1.0)
    lo_s, hi_s = used / 10.0, used * 10.0
    f_lo, f_hi = bias(lo_s), bias(hi_s)
    if not (f_lo < 0 < f_hi):
        return SigmaM(
            used,
            None,
            n,
            bias_used,
            "no sigma_m within a factor of ten makes the partials unbiased",
        )
    for _ in range(80):
        mid = np.sqrt(lo_s * hi_s)
        if bias(mid) < 0:
            lo_s = mid
        else:
            hi_s = mid
    note = "" if n >= 30 else f"only {n} partials: a rough estimate"
    return SigmaM(used, float(np.sqrt(lo_s * hi_s)), n, bias_used, note)


@dataclass
class Prepared:
    """What from_files keeps, for --worst and the sigma_m suggestion."""

    table: refl.ReflectionTable
    rows: dict[str, np.ndarray] = field(
        default_factory=dict
    )  # each kind's rows of the table
    sigma: dict[str, SigmaM] = field(default_factory=dict)


def from_files(
    expt_path,
    refl_path,
    kinds=("prf", "sum"),
    min_i_sigma: float = 5.0,
    anomalous: bool = False,
    others_clean: bool = True,
    method: str = "mean",
) -> tuple[list[Result], Prepared]:
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
    lp = column("lp", 1.0).astype(float)
    qe = column("qe", 1.0).astype(float)
    qe = np.where(qe > 0, qe, 1.0)
    partiality_col = column("partiality", 1.0).astype(float)
    measured = column("profile.measured", 1.0).astype(float)
    zeta = column("zeta", np.nan).astype(float)
    px = (
        column("xyzcal.px").reshape(-1, 3)
        if "xyzcal.px" in c
        else np.zeros((t.nrows, 3))
    )
    mm = column("xyzcal.mm").reshape(-1, 3) if "xyzcal.mm" in c else None
    bbox = column("bbox").reshape(-1, 6) if "bbox" in c else None
    z = px[:, 2]
    d = column("d", np.nan).astype(float)
    images = (
        float(e0.scan.num_images) if e0.scan is not None else float(np.nanmax(z) + 1)
    )
    from_start, from_end = z, images - z

    keep = (flags & (OUTLIER_IN_SCALING | EXCLUDED_FOR_SCALING)) == 0
    if np.any(flags & SCALED):
        keep &= (flags & SCALED) != 0
    keys = asu_keys(hkl, hall, anomalous)
    full = partiality_col >= CLEAN_PARTIALITY
    unmasked = (measured >= CLEAN_MEASURED) & ((flags & FOREGROUND_BAD) == 0)
    interior = (from_start >= CLEAN_EDGE_IMAGES) & (from_end >= CLEAN_EDGE_IMAGES)
    clean_base = full & unmasked & interior
    fg_bad = np.where((flags & FOREGROUND_BAD) != 0, "yes", "no")

    # Frames to phi, from the table itself: phi is exactly linear in z.
    phi_of_z = None
    if mm is not None and np.ptp(z) > 0:
        slope, intercept = np.polyfit(z, mm[:, 2], 1)
        phi_of_z = (slope, intercept)

    prepared = Prepared(t)
    results = []
    for kind in kinds:
        bit = PRF if kind == "prf" else SUM
        value = column(f"intensity.{kind}.value").astype(float)
        var = column(f"intensity.{kind}.variance").astype(float)
        rows = keep & ((flags & bit) != 0) & (scale > 0) & (partiality_col > 0)
        factor = lp[rows] / qe[rows] / partiality_col[rows] / scale[rows]
        intensity = value[rows] * factor
        variance = var[rows] * factor**2
        isig = intensity / np.sqrt(np.where(variance > 0, variance, np.nan))
        others = others_clean
        explanatories = [
            Explanatory(
                "partiality",
                partiality_col[rows],
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
                quantiles=6,
                condition=clean_base[rows] if others else None,
            ),
            Explanatory(
                "I/sigma, a control",
                isig,
                [0, 10, 20, 50, 100, np.inf],
                condition=clean_base[rows] if others else None,
            ),
        ]
        result = analyse(
            keys[rows],
            intensity,
            variance,
            clean_base[rows],
            explanatories,
            min_i_sigma,
            kind,
            method,
        )
        results.append(result)
        prepared.rows[kind] = np.nonzero(rows)[0]

        # The sigma_m the partials suggest: those counted, otherwise clean,
        # partial as integrated.
        if phi_of_z is not None and bbox is not None and np.isfinite(zeta).any():
            p = partiality_col[rows]
            sel = (
                result.counted_rows
                & unmasked[rows]
                & (p > 0.05)
                & (p < CLEAN_PARTIALITY)
                & np.isfinite(zeta[rows])
            )
            a = phi_of_z[1] + phi_of_z[0] * bbox[rows, 4].astype(float)
            b = phi_of_z[1] + phi_of_z[0] * bbox[rows, 5].astype(float)
            counts = intensity * p  # corrected for everything but the partiality
            prepared.sigma[kind] = suggest_sigma_m(
                counts[sel],
                result.reference[sel],
                p[sel],
                mm[rows, 2][sel],
                zeta[rows][sel],
                np.minimum(a, b)[sel],
                np.maximum(a, b)[sel],
            )
    return results, prepared


def worst(
    prepared: Prepared, result: Result, kind: str, n: int
) -> refl.ReflectionTable:
    """The n counted observations furthest from their clean equivalents, as a
    table dials.image_viewer will draw: our rows, with the reference and the
    relative difference beside them."""
    rows = prepared.rows[kind]
    counted = np.nonzero(result.counted_rows)[0]
    order = counted[np.argsort(-np.abs(result.rel[counted]))[:n]]
    order.sort()
    chosen = rows[order]
    t = prepared.table
    out = refl.ReflectionTable(nrows=len(chosen))
    for name, values in t.columns.items():
        out.columns[name] = values[chosen]
        if name in t.types:
            out.types[name] = t.types[name]
    out.identifiers = dict(t.identifiers)
    out.columns["equivalents.reference"] = result.reference[order].astype(float)
    out.types["equivalents.reference"] = "double"
    out.columns["equivalents.difference"] = result.rel[order].astype(float)
    out.types["equivalents.difference"] = "double"
    return out


def report(
    results, min_i_sigma: float, others_clean: bool = True, prepared=None
) -> str:
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
                f"    {'bin':<18} {'n':>8} {'mean %':>9} {'+-':>6} {'median %':>9} {'sums %':>8}"
            )
            for b in table.bins:
                if b.n == 0:
                    lines.append(f"    {b.label:<18} {0:>8}")
                    continue
                sem = f"{100 * b.sem:6.2f}" if np.isfinite(b.sem) else f"{'':>6}"
                lines.append(
                    f"    {b.label:<18} {b.n:>8} {100 * b.mean:+9.2f} {sem} "
                    f"{100 * b.median:+9.2f} {100 * b.sums:+8.2f}"
                )
        if prepared is not None and r.intensity in prepared.sigma:
            s = prepared.sigma[r.intensity]
            lines.append("")
            lines.append("  sigma_m, from the partials")
            lines.append(
                f"    integrated with {s.used:.4f} deg; partials, median {100 * s.bias_used:+.2f} %, "
                f"over {s.n}"
            )
            if s.suggested is not None:
                lines.append(
                    f"    suggested: {s.suggested:.4f} deg, where their median bias is zero"
                    f" -- to try, mxi_integrate --sigma-m {s.suggested:.4f}"
                )
            if s.note:
                lines.append(f"    ({s.note})")
        lines.append("")
    return "\n".join(lines)
