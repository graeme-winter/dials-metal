"""Binning a disagreement against what might explain it."""

import numpy as np
import pytest

from mxeq import refl, trends


def table(n, seed=0, ratio=1.0, ratio_with_d=None):
    """A small integrated table, with an optional systematic in the ratio."""
    rng = np.random.default_rng(seed)
    hkl = np.stack(
        [np.arange(n) % 17 - 8, (np.arange(n) // 17) % 17 - 8, np.arange(n) // 289],
        axis=1,
    ).astype(np.int64)
    d = np.linspace(1.5, 8.0, n)
    intensity = rng.gamma(2.0, 200.0, n)
    if ratio_with_d is not None:
        intensity = intensity * (1.0 + ratio_with_d * (d - d.mean()) / np.ptp(d))
    t = refl.ReflectionTable(nrows=n)
    t.columns["miller_index"] = hkl
    t.columns["entering"] = np.zeros(n, dtype=np.int64)
    t.columns["xyzcal.px"] = np.stack(
        [rng.uniform(0, 2000, n), rng.uniform(0, 2000, n), np.arange(n) % 300.0],
        axis=1,
    )
    t.columns["d"] = d
    t.columns["intensity.sum.value"] = intensity * ratio
    t.columns["intensity.sum.variance"] = np.maximum(intensity, 1.0)
    t.columns["zeta"] = rng.uniform(0.1, 1.0, n)
    t.columns["partiality"] = np.full(n, 0.999)
    return t


def test_a_flat_ratio_shows_as_flat():
    # The control: two files differing by a constant scale should show that
    # scale in every bin of every variable and no trend in any of them.
    a = table(600, seed=1, ratio=0.9)
    b = table(600, seed=1, ratio=1.0)
    report = trends.compare(a, b, ["intensity.sum.value"], n_bins=5)
    assert "matched 600" in report
    ratios = [
        float(line.split()[-2])
        for line in report.splitlines()
        if line.startswith("  ") and len(line.split()) == 7 and "from" not in line
    ]
    assert ratios, report
    assert all(abs(r - 0.9) < 1e-9 for r in ratios), report


def test_a_ratio_that_trends_with_resolution_shows_as_a_trend():
    # The thing the tool is for: a disagreement that depends on resolution must
    # be visible as a monotonic column, not averaged away into one number.
    a = table(600, seed=2, ratio_with_d=0.6)
    b = table(600, seed=2)
    report = trends.compare(a, b, ["intensity.sum.value"], n_bins=5)
    # Only the resolution block: everything after it is another variable's
    # table, and reading into those was this test's own first bug.
    after = report.split("against resolution")[1].splitlines()
    block = []
    for line in after[1:]:
        if not line.startswith("  "):
            break
        block.append(line)
    ratios = [
        float(line.split()[-2])
        for line in block
        if len(line.split()) == 7 and "from" not in line
    ]
    assert len(ratios) >= 4, report
    # Rows are highest resolution last, so the ratio should fall down the table.
    assert ratios[0] > ratios[-1] + 0.2, ratios


def test_nothing_matching_says_what_to_check():
    # The failure that has happened repeatedly in this work: comparing a file
    # written before a reindexing step with one written after.
    a = table(100, seed=3)
    b = table(100, seed=3)
    b.columns["miller_index"] = -b.columns["miller_index"] - 7
    report = trends.compare(a, b, ["intensity.sum.value"])
    assert "Nothing matched" in report
    assert "reindexing" in report


def test_a_missing_column_is_said_rather_than_skipped():
    a = table(200, seed=4)
    b = table(200, seed=4)
    report = trends.compare(a, b, ["intensity.prf.value"], n_bins=4)
    assert "not in both files" in report


def test_bins_hold_equal_populations():
    # Equal-width bins on a quantity like I/sigma put almost everything in the
    # first one and say nothing.
    values = np.concatenate([np.random.default_rng(5).gamma(1.0, 5.0, 1000)])
    edges = trends.bin_edges(values, 10)
    counts = np.histogram(values, edges)[0]
    assert counts.min() > 0.5 * counts.max()
