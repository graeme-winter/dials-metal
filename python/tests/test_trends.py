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
    # It refuses rather than reporting on nothing, and says where to look.
    with pytest.raises(ValueError, match="reindex"):
        trends.compare(a, b, ["intensity.sum.value"])


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
    out, report = disagree.select(a, b, factor=2.0, floor=1.0)
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
    out, _ = disagree.select(a, b, factor=2.0, floor=5.0)
    assert out.nrows == 0
    # But the same ratio between numbers that matter is a disagreement.
    a.columns["intensity.sum.value"] = np.full(200, 600.0)
    b.columns["intensity.sum.value"] = np.full(200, 200.0)
    out, _ = disagree.select(a, b, factor=2.0, floor=5.0)
    assert out.nrows == 200


def test_a_sign_flip_disagrees_however_small():
    from mxeq import disagree

    a = table(100, seed=23)
    b = table(100, seed=23)
    a.columns["intensity.sum.value"] = np.full(100, 20.0)
    b.columns["intensity.sum.value"] = np.full(100, -20.0)
    out, _ = disagree.select(a, b, factor=100.0, floor=1.0)
    assert out.nrows == 100


def test_absolute_and_relative_find_different_reflections():
    from mxeq import disagree

    # A constant offset of fifty counts: a huge relative error on a weak
    # reflection and a trivial one on a strong reflection, but the same
    # absolute error on both. The two criteria should disagree about which
    # reflections are the problem, which is the whole reason for having both.
    a = table(400, seed=51)
    b = table(400, seed=51)
    a.columns["intensity.sum.value"] = b.columns["intensity.sum.value"] + 50.0

    relative, _ = disagree.select(a, b, factor=2.0, floor=1.0)
    absolute, _ = disagree.select(a, b, difference=25.0)
    # Every reflection is off by fifty counts, so absolute catches all of them.
    assert absolute.nrows == 400
    # And relative catches only the ones where fifty counts is a lot.
    assert 0 < relative.nrows < 400
    weak = b.columns["intensity.sum.value"] < 50.0
    assert relative.nrows <= int(weak.sum()) + 20


def test_the_criteria_narrow_rather_than_widen():
    from mxeq import disagree

    a = table(400, seed=52)
    b = table(400, seed=52)
    a.columns["intensity.sum.value"] = b.columns["intensity.sum.value"] * 3.0

    only_factor, _ = disagree.select(a, b, factor=2.0, floor=1.0)
    both, _ = disagree.select(a, b, factor=2.0, difference=500.0)
    assert both.nrows < only_factor.nrows
    assert both.nrows > 0


def test_sigma_needs_the_variance_and_says_so():
    from mxeq import disagree

    a = table(100, seed=53)
    b = table(100, seed=53)
    del a.columns["intensity.sum.variance"]
    with pytest.raises(ValueError, match="variance"):
        disagree.select(a, b, sigmas=3.0)


def test_a_reflection_seen_many_times_is_matched_every_time():
    from mxeq import match

    # A Miller index is not a key on a scan that goes round more than once. Ten
    # rotations record every reflection about twenty times, all sharing their
    # index and their entering flag, and matching on those alone pairs one of
    # each and throws the rest away: on a real ten rotation sweep that left
    # 647334 of 6977214 rows matched and called 722427 keys duplicates.
    turns = 10
    hkl = np.tile(np.array([[1, 2, 3], [4, 5, 6]], dtype=np.int64), (turns, 1))
    entering = np.zeros(len(hkl), dtype=np.int64)
    # One observation per turn, 3600 frames apart, as a real sweep gives.
    z = np.repeat(np.arange(turns) * 3600.0, 2) + np.tile([10.0, 20.0], turns)

    ia, ib, unpartnered = match.match_observations(
        hkl, entering, z, hkl, entering, z + 0.1
    )
    assert len(ia) == len(hkl), (len(ia), len(hkl))
    assert unpartnered == 0
    # And each is paired with its own turn, not with another turn's copy.
    assert np.array_equal(z[ia], z[ib])


def test_an_observation_with_no_partner_is_left_out_rather_than_mispaired():
    from mxeq import match

    hkl = np.array([[1, 2, 3]] * 3, dtype=np.int64)
    entering = np.zeros(3, dtype=np.int64)
    a_z = np.array([10.0, 3610.0, 7210.0])
    b_z = np.array([10.0, 7210.0])          # the middle turn is missing
    ia, ib, unpartnered = match.match_observations(
        hkl, entering, a_z, hkl[:2], entering[:2], b_z
    )
    assert len(ia) == 2
    assert np.allclose(a_z[ia], b_z[ib])
    assert unpartnered == 1


def test_the_radius_is_a_sanity_check_and_not_a_discriminator():
    from mxeq import match

    # Within a group the observations are a whole turn apart, so there is
    # nothing for a tight radius to protect against and plenty for it to lose.
    # On a real pair, half a frame matched 81 per cent of the smaller table and
    # two frames matched 96, the 99th percentile of the frame difference being
    # 1.6 against a median of 0.001.
    turns = 6
    hkl = np.tile(np.array([[1, 2, 3]], dtype=np.int64), (turns, 1))
    entering = np.zeros(turns, dtype=np.int64)
    z = np.arange(turns) * 3600.0
    # One side predicts each observation 1.4 frames later: well inside a turn
    # and well outside half a frame.
    offset = z + 1.4

    tight = match.match_observations(hkl, entering, z, hkl, entering, offset,
                                     radius=0.5)
    assert len(tight[0]) == 0, "a tight radius loses every one of them"

    loose = match.match_observations(hkl, entering, z, hkl, entering, offset)
    assert len(loose[0]) == turns
    # And each is still paired with its own turn rather than a neighbouring one.
    assert np.allclose(offset[loose[1]] - z[loose[0]], 1.4)


def test_a_radius_wider_than_a_turn_still_pairs_correctly():
    from mxeq import match

    # Because the walk is in frame order, not nearest-first: even a radius that
    # spans several turns cannot pair an observation with the wrong turn's copy
    # while its own is available.
    turns = 5
    hkl = np.tile(np.array([[2, 0, 1]], dtype=np.int64), (turns, 1))
    entering = np.zeros(turns, dtype=np.int64)
    z = np.arange(turns) * 100.0
    ia, ib, _ = match.match_observations(hkl, entering, z, hkl, entering,
                                         z + 0.2, radius=10000.0)
    assert len(ia) == turns
    assert np.allclose(z[ia] + 0.2, (z + 0.2)[ib])


def test_unpartnered_observations_are_explained_by_cause():
    from mxeq import match

    # Each cause planted once, so the explanation has to tell them apart.
    turn = 3600.0
    a_hkl = np.array([[1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0], [5, 0, 0]],
                     dtype=np.int64)
    a_ent = np.zeros(5, dtype=np.int64)
    a_z = np.array([100.0, 100.0, 100.0, 100.0, 100.0])

    # B: 1 absent entirely; 2 on the other side of the sphere; 3 a turn later;
    # 4 thirty frames away; 5 a true match so the matcher has something.
    b_hkl = np.array([[2, 0, 0], [3, 0, 0], [4, 0, 0], [5, 0, 0]], dtype=np.int64)
    b_ent = np.array([1, 0, 0, 0], dtype=np.int64)
    b_z = np.array([100.5, 100.0 + turn, 130.0, 100.2])

    ia, ib, _ = match.match_observations(a_hkl, a_ent, a_z, b_hkl, b_ent, b_z,
                                         radius=5.0)
    assert len(ia) == 1
    why = match.explain_unpartnered(a_hkl, a_ent, a_z, b_hkl, b_ent, b_z,
                                    ia, 5.0, turn)
    assert why["not predicted by the other"] == 1
    assert why["entering flag disagrees"] == 1
    assert why["a whole number of turns apart"] == 1
    assert why["further apart than the radius"] == 1
    assert why["within the radius yet unpaired"] == 0


def test_the_turn_is_measured_from_the_repeats():
    from mxeq import match

    turns = 8
    hkl = np.tile(np.array([[1, 2, 3], [3, 2, 1], [0, 1, 4]], dtype=np.int64),
                  (turns, 1))
    entering = np.zeros(len(hkl), dtype=np.int64)
    z = np.repeat(np.arange(turns) * 3600.0, 3) + np.tile([5.0, 50.0, 500.0], turns)
    assert match.estimate_turn(hkl, entering, z) == pytest.approx(3600.0)
    # And less than a turn has no repeats to measure it from.
    assert match.estimate_turn(hkl[:3], entering[:3], z[:3]) is None
