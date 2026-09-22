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


def ratios_under(report, heading):
    """The ratio column of one block, found by its header rather than by
    counting from the end: adding a column to the table should not break every
    test that reads it, and once did."""
    after = report.split(heading)[1].splitlines()
    header = after[1].split()
    where = header.index("ratio")
    out = []
    for line in after[2:]:
        if not line.startswith("  "):
            break
        parts = line.split()
        if len(parts) != len(header):
            break
        out.append(float(parts[where]))
    return out


def test_a_flat_ratio_shows_as_flat():
    # The control: two files differing by a constant scale should show that
    # scale in every bin of every variable and no trend in any of them.
    a = table(600, seed=1, ratio=0.9)
    b = table(600, seed=1, ratio=1.0)
    report = trends.compare(a, b, ["intensity.sum.value"], n_bins=5)
    assert "matched 600" in report
    ratios = ratios_under(report, "against resolution")
    assert ratios, report
    assert all(abs(r - 0.9) < 1e-9 for r in ratios), report


def test_a_ratio_that_trends_with_resolution_shows_as_a_trend():
    # The thing the tool is for: a disagreement that depends on resolution must
    # be visible as a monotonic column, not averaged away into one number.
    a = table(600, seed=2, ratio_with_d=0.6)
    b = table(600, seed=2)
    report = trends.compare(a, b, ["intensity.sum.value"], n_bins=5)
    ratios = ratios_under(report, "against resolution")
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


def test_spearman_survives_what_pearson_does_not():
    # The reading this was added for: a bin where almost everything agrees and
    # a handful of values are wildly wrong. Pearson collapses, Spearman does
    # not, and the difference between them is the diagnosis.
    rng = np.random.default_rng(11)
    truth = rng.uniform(100.0, 200.0, 500)
    same = truth * rng.normal(1.0, 0.01, 500)
    same[:5] = truth[:5] * 400.0
    pearson, _ = __import__("mxeq.stats", fromlist=["stats"]).correlation(truth, same)
    rho = trends.spearman(truth, same)
    assert pearson < 0.5, pearson
    assert rho > 0.9, rho


def test_the_outlier_count_finds_planted_outliers():
    a = table(1000, seed=7)
    b = table(1000, seed=7)
    spoiled = a.columns["intensity.sum.value"].copy()
    spoiled[:20] *= 50.0
    a.columns["intensity.sum.value"] = spoiled
    report = trends.compare(a, b, ["intensity.sum.value"], n_bins=4)
    after = report.split("against resolution")[1].splitlines()
    header = after[1].split()
    where = header.index("out")
    counts = []
    for line in after[2:]:
        parts = line.split()
        if len(parts) != len(header):
            break
        counts.append(int(parts[where]))
    assert sum(counts) >= 15, report


def test_disagreeing_reflections_come_out_as_a_table():
    from mxeq import disagree

    a = table(400, seed=21)
    b = table(400, seed=21)
    spoiled = a.columns["intensity.sum.value"].copy()
    spoiled[:30] *= 40.0
    a.columns["intensity.sum.value"] = spoiled
    out, report = disagree.select(a, b, factor=2.0, absolute=1.0)
    assert out.nrows == 30, report
    # The reference's value and the ratio travel with it, so the viewer's table
    # shows both without a second file being opened.
    assert "reference.intensity" in out.columns
    assert "disagreement.ratio" in out.columns
    assert np.allclose(out.columns["disagreement.ratio"], 40.0)
    # And everything the image viewer needs to draw a box.
    for needed in ("miller_index", "xyzcal.px", "bbox" if "bbox" in a.columns else "d"):
        assert needed in out.columns


def test_pairs_that_are_both_tiny_are_not_called_a_disagreement():
    from mxeq import disagree

    a = table(200, seed=22)
    b = table(200, seed=22)
    a.columns["intensity.sum.value"] = np.full(200, 0.6)
    b.columns["intensity.sum.value"] = np.full(200, 0.2)
    out, _ = disagree.select(a, b, factor=2.0, absolute=5.0)
    assert out.nrows == 0
    # But the same ratio between numbers that matter is a disagreement.
    a.columns["intensity.sum.value"] = np.full(200, 600.0)
    b.columns["intensity.sum.value"] = np.full(200, 200.0)
    out, _ = disagree.select(a, b, factor=2.0, absolute=5.0)
    assert out.nrows == 200


def test_a_sign_flip_disagrees_however_small():
    from mxeq import disagree

    a = table(100, seed=23)
    b = table(100, seed=23)
    a.columns["intensity.sum.value"] = np.full(100, 20.0)
    b.columns["intensity.sum.value"] = np.full(100, -20.0)
    out, _ = disagree.select(a, b, factor=100.0, absolute=1.0)
    assert out.nrows == 100
