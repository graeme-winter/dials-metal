"""The metrics the checks report.

Version one of this package sets no thresholds.  Every function here returns a
distribution or a summary of one, and whether a number is acceptable is a
judgement to be made after looking at it on real data, not encoded in advance
by someone guessing.

One statistic needs a warning attached wherever it appears.  The pull,
(I_a - I_b) / sqrt(sigma_a^2 + sigma_b^2), is the right quantity for two
*independent* measurements of the same reflection, and its width should then be
one.  Here the two intensities come from the same photons processed twice, so
their errors are strongly correlated and the pull width will be far below one
even when the two pipelines disagree meaningfully.  A width of 0.2 is not a
pass.  Use the pull for its *shape* -- symmetry, and how heavy the tails are --
and the relative difference for its size.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

try:
    import gemmi
except ImportError:  # pragma: no cover
    gemmi = None

PERCENTILES = (1.0, 5.0, 25.0, 50.0, 75.0, 95.0, 99.0)


@dataclass
class Summary:
    """A distribution, described without assuming it is Gaussian."""

    n: int
    mean: float = float("nan")
    median: float = float("nan")
    std: float = float("nan")
    #: Median absolute deviation scaled to be comparable with a standard
    #: deviation for Gaussian data.  Reported alongside ``std`` because the
    #: gap between the two is the tail weight.
    robust_std: float = float("nan")
    minimum: float = float("nan")
    maximum: float = float("nan")
    percentiles: dict[float, float] = field(default_factory=dict)

    def line(self, name: str, fmt: str = "{: .4g}") -> str:
        if self.n == 0:
            return f"  {name:<28} (empty)"
        p = self.percentiles
        return (
            f"  {name:<28} n {self.n:>8}"
            f"  median {fmt.format(self.median)}"
            f"  mean {fmt.format(self.mean)}"
            f"  sd {fmt.format(self.std)}"
            f"  rsd {fmt.format(self.robust_std)}"
            f"  [1% {fmt.format(p.get(1.0, float('nan')))}"
            f" 99% {fmt.format(p.get(99.0, float('nan')))}]"
        )


def describe(values: np.ndarray) -> Summary:
    values = np.asarray(values, dtype=float).ravel()
    values = values[np.isfinite(values)]
    if values.size == 0:
        return Summary(n=0)
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    return Summary(
        n=int(values.size),
        mean=float(np.mean(values)),
        median=median,
        std=float(np.std(values, ddof=1)) if values.size > 1 else 0.0,
        robust_std=1.4826 * mad,
        minimum=float(values.min()),
        maximum=float(values.max()),
        percentiles={p: float(np.percentile(values, p)) for p in PERCENTILES},
    )


def correlation(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Pearson and Spearman correlation, both robust to an empty input.

    Spearman is reported next to Pearson because intensity distributions are
    long-tailed: a Pearson coefficient over unbinned intensities is dominated
    by the few strongest reflections and can look excellent while the bulk
    disagrees.
    """
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if a.size < 3:
        return float("nan"), float("nan")

    def pearson(x, y):
        x = x - x.mean()
        y = y - y.mean()
        denominator = float(np.sqrt((x * x).sum() * (y * y).sum()))
        return float((x * y).sum() / denominator) if denominator > 0 else float("nan")

    rank_a = np.argsort(np.argsort(a)).astype(float)
    rank_b = np.argsort(np.argsort(b)).astype(float)
    return pearson(a, b), pearson(rank_a, rank_b)


def pull(
    i_a: np.ndarray, sigma_a: np.ndarray, i_b: np.ndarray, sigma_b: np.ndarray
) -> np.ndarray:
    """(I_a - I_b) / sqrt(sigma_a^2 + sigma_b^2).  See the module docstring."""
    i_a, i_b = np.asarray(i_a, dtype=float), np.asarray(i_b, dtype=float)
    va = np.asarray(sigma_a, dtype=float) ** 2
    vb = np.asarray(sigma_b, dtype=float) ** 2
    denominator = np.sqrt(va + vb)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (i_a - i_b) / denominator
    return np.where(denominator > 0, out, np.nan)


def relative_difference(i_a: np.ndarray, i_b: np.ndarray) -> np.ndarray:
    """2(I_a - I_b) / (|I_a| + |I_b|), which stays finite through zero."""
    i_a, i_b = np.asarray(i_a, dtype=float), np.asarray(i_b, dtype=float)
    denominator = np.abs(i_a) + np.abs(i_b)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = 2.0 * (i_a - i_b) / denominator
    return np.where(denominator > 0, out, np.nan)


def d_spacing(hkl: np.ndarray, setting: np.ndarray) -> np.ndarray:
    """Resolution of each reflection, from the A matrix with a*, b*, c* as columns."""
    hkl = np.atleast_2d(np.asarray(hkl, dtype=float))
    q = hkl @ np.asarray(setting, dtype=float).T
    length = np.linalg.norm(q, axis=1)
    with np.errstate(divide="ignore"):
        return np.where(length > 0, 1.0 / length, np.inf)


def resolution_bins(d: np.ndarray, n_bins: int = 10) -> tuple[np.ndarray, np.ndarray]:
    """Bin assignments and bin edges, equal in reciprocal volume.

    Equal-volume shells rather than equal counts, so that a shell means the
    same thing between two runs that selected different numbers of
    reflections.  Edges are returned in d, high resolution last.
    """
    d = np.asarray(d, dtype=float)
    finite = d[np.isfinite(d) & (d > 0)]
    if finite.size == 0:
        return np.zeros(len(d), dtype=int), np.array([np.inf, 0.0])
    lo, hi = finite.max(), finite.min()  # lo is low resolution, large d
    edges_cubed = np.linspace(1.0 / lo**3, 1.0 / hi**3, n_bins + 1)
    edges = edges_cubed ** (-1.0 / 3.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        inverse = np.where(np.isfinite(d) & (d > 0), 1.0 / d**3, np.nan)
    which = np.digitize(inverse, edges_cubed[1:-1], right=False)
    which = np.where(np.isfinite(inverse), which, -1)
    return which.astype(int), edges


def to_asu(hkl: np.ndarray, hall: str | None) -> np.ndarray | None:
    """Map Miller indices into the reciprocal asymmetric unit.

    Returns ``None`` if the space group is unavailable, which the callers treat
    as "cannot merge" rather than as "no symmetry".
    """
    if hall is None or gemmi is None:
        return None
    try:
        sg = gemmi.find_spacegroup_by_ops(gemmi.symops_from_hall(hall))
    except (ValueError, RuntimeError):
        return None
    if sg is None:
        return None
    asu = gemmi.ReciprocalAsu(sg)
    ops = sg.operations()
    hkl = np.atleast_2d(np.asarray(hkl, dtype=int))
    out = np.empty_like(hkl)
    for i, h in enumerate(hkl):
        out[i] = asu.to_asu((int(h[0]), int(h[1]), int(h[2])), ops)[0]
    return out


@dataclass
class Merged:
    """Merged intensities keyed by ASU index, with the multiplicity kept."""

    hkl: np.ndarray
    intensity: np.ndarray
    sigma: np.ndarray
    multiplicity: np.ndarray
    #: Half-dataset means, for CC-half.  NaN where multiplicity is below two.
    half_a: np.ndarray
    half_b: np.ndarray
    r_meas: float


def merge(
    asu_hkl: np.ndarray,
    intensity: np.ndarray,
    sigma: np.ndarray,
    seed: int = 0,
) -> Merged:
    """Weighted merge with half-dataset splits and Rmeas.

    The split is random with a fixed seed rather than by parity of the
    observation index, because observation order is a pipeline-dependent thing
    and splitting on it would make CC-half compare two pipelines' orderings
    instead of their data.
    """
    asu_hkl = np.atleast_2d(np.asarray(asu_hkl, dtype=np.int64))
    intensity = np.asarray(intensity, dtype=float)
    sigma = np.asarray(sigma, dtype=float)

    unique, inverse = np.unique(asu_hkl, axis=0, return_inverse=True)
    inverse = inverse.ravel()
    n_groups = len(unique)

    weight = np.where(sigma > 0, 1.0 / np.maximum(sigma, 1e-12) ** 2, 0.0)
    sum_w = np.bincount(inverse, weights=weight, minlength=n_groups)
    sum_wi = np.bincount(inverse, weights=weight * intensity, minlength=n_groups)
    multiplicity = np.bincount(inverse, minlength=n_groups)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(sum_w > 0, sum_wi / sum_w, np.nan)
        merged_sigma = np.where(
            sum_w > 0, 1.0 / np.sqrt(np.maximum(sum_w, 1e-300)), np.nan
        )

    # Rmeas: the multiplicity-corrected merging residual.  Groups of one
    # contribute nothing to either sum, which is what the n/(n-1) factor is
    # for and why they are excluded rather than given a residual of zero.
    deviation = np.abs(intensity - mean[inverse])
    many = multiplicity[inverse] > 1
    with np.errstate(invalid="ignore", divide="ignore"):
        factor = np.sqrt(
            multiplicity[inverse] / np.maximum(multiplicity[inverse] - 1.0, 1.0)
        )
    numerator = float(np.nansum((factor * deviation)[many]))
    denominator = float(np.nansum(np.abs(intensity)[many]))
    r_meas = numerator / denominator if denominator > 0 else float("nan")

    rng = np.random.default_rng(seed)
    side = rng.integers(0, 2, size=len(intensity))
    halves = []
    for s in (0, 1):
        pick = side == s
        w = np.where(pick, weight, 0.0)
        sw = np.bincount(inverse, weights=w, minlength=n_groups)
        swi = np.bincount(inverse, weights=w * intensity, minlength=n_groups)
        with np.errstate(invalid="ignore", divide="ignore"):
            halves.append(np.where(sw > 0, swi / sw, np.nan))

    return Merged(
        hkl=unique,
        intensity=mean,
        sigma=merged_sigma,
        multiplicity=multiplicity,
        half_a=halves[0],
        half_b=halves[1],
        r_meas=r_meas,
    )


def cc_half(merged: Merged) -> float:
    ok = np.isfinite(merged.half_a) & np.isfinite(merged.half_b)
    if ok.sum() < 3:
        return float("nan")
    return correlation(merged.half_a[ok], merged.half_b[ok])[0]
