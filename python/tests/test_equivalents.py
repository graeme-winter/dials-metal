"""mxeq equivalents measures a planted bias as itself.

Each observation is compared with the clean equivalents of its reflection --
fully recorded, nothing masked -- leaving itself out, so that a bias in some
observations does not leak into the reference the others are measured
against. Planted: partials 10 per cent high and observations crossing a mask
3 per cent low, a fifth and a tenth of all observations; the tool must read
+10, -3 and 0, where measuring against every equivalent would read the partials
low and the clean high. And the symmetry grouping must partition indices as
gemmi's own reciprocal asymmetric unit does.
"""

import gemmi
import numpy as np
import pytest

from mxeq import equivalents as eq


@pytest.mark.parametrize(
    "hall", ["P 2ac 2ab", "I 2 2 3", "P 4nw 2abw", "P 6c 2c", '-R 3 2"', "P 1"]
)
def test_equivalent_indices_share_a_key_as_gemmi_groups_them(hall):
    rng = np.random.default_rng(5)
    hkl = rng.integers(-12, 13, size=(4000, 3))
    hkl = hkl[np.abs(hkl).sum(axis=1) > 0]
    ops = gemmi.symops_from_hall(hall)
    sg = gemmi.find_spacegroup_by_ops(ops)
    asu = gemmi.ReciprocalAsu(sg)
    theirs = eq._encode(
        np.array([asu.to_asu(list(map(int, h)), sg.operations())[0] for h in hkl])
    )
    mine = eq.asu_keys(hkl, hall)
    _, a = np.unique(theirs, return_inverse=True)
    _, b = np.unique(mine, return_inverse=True)
    assert len(set(zip(a.tolist(), b.tolist()))) == a.max() + 1 == b.max() + 1


def planted():
    rng = np.random.default_rng(11)
    groups, per = 3000, 8
    keys = np.repeat(np.arange(groups), per)
    truth = np.repeat(1000.0 * 10 ** rng.uniform(0, 2, groups), per)
    n = len(keys)
    kind = rng.uniform(size=n)
    partial = kind < 0.2
    masked = (kind >= 0.2) & (kind < 0.3)
    partiality = np.where(partial, rng.uniform(0.5, 0.89, n), 1.0)
    measured = np.where(masked, rng.uniform(0.8, 0.94, n), 1.0)
    bias = np.where(partial, 1.10, np.where(masked, 0.97, 1.0))
    sigma = 0.01 * truth
    intensity = truth * bias + sigma * rng.normal(size=n)
    variance = sigma**2
    clean = ~partial & ~masked
    explanatories = [
        eq.Explanatory(
            "partiality", partiality, [0.0, 0.9, 0.99, 1.0001], condition=~masked
        ),
        eq.Explanatory(
            "measured", measured, [0.0, 0.95, 0.999, 1.0001], condition=~partial
        ),
    ]
    return keys, intensity, variance, clean, explanatories


def bins(result, name):
    return {b.label: b for t in result.tables if t.explanatory == name for b in t.bins}


def test_a_planted_bias_reads_as_itself_against_clean_equivalents():
    keys, intensity, variance, clean, explanatories = planted()
    r = eq.analyse(keys, intensity, variance, clean, explanatories, min_i_sigma=5.0)
    p = bins(r, "partiality")
    m = bins(r, "measured")
    assert abs(p["[0, 0.9)"].median - 0.10) < 0.005, p["[0, 0.9)"]
    assert abs(m["[0, 0.95)"].median + 0.03) < 0.005, m["[0, 0.95)"]
    assert abs(p["[0.99, 1.0001)"].median) < 0.003, p["[0.99, 1.0001)"]


def test_against_every_equivalent_the_bias_would_leak_into_the_reference():
    # The reason for clean references: with all equivalents as the reference,
    # the partials read low and the clean observations high.
    keys, intensity, variance, clean, explanatories = planted()
    r = eq.analyse(
        keys, intensity, variance, np.ones_like(clean), explanatories, min_i_sigma=5.0
    )
    p = bins(r, "partiality")
    assert p["[0, 0.9)"].median < 0.09
    assert p["[0.99, 1.0001)"].median < -0.005


def test_an_unweighted_reference_reads_unbiased_data_as_unbiased():
    # Variances from each observation's own counts, as integration's are:
    # weighted by them, the equivalents that came out low weigh most and the
    # reference reads low -- an unbiased data set reads high. Unweighted, the
    # sums read zero.
    rng = np.random.default_rng(3)
    groups, per = 20000, 6
    truth = np.repeat(50.0 * 10 ** rng.uniform(0, 3, groups), per)
    keys = np.repeat(np.arange(groups), per)
    sd = np.sqrt(truth + 30.0 + (0.03 * truth) ** 2)
    intensity = truth + sd * rng.normal(size=len(truth))
    variance = np.maximum(intensity, 0) + 30.0 + (0.03 * intensity) ** 2
    clean = np.ones(len(truth), dtype=bool)
    mean = eq.analyse(
        keys, intensity, variance, clean, [], min_i_sigma=10.0, method="mean"
    )
    weighted = eq.analyse(
        keys, intensity, variance, clean, [], min_i_sigma=10.0, method="weighted"
    )
    assert abs(mean.overall.sums) < 0.001, mean.overall
    assert weighted.overall.sums > mean.overall.sums + 0.001, (
        weighted.overall,
        mean.overall,
    )


def test_the_partials_suggest_the_sigma_m_that_unbiases_them():
    # Partials integrated with sigma_m 0.08 degrees where the crystal's is 0.05:
    # their partiality computed too small, so dividing by it overshoots. The
    # tool recovers the 0.08 used from each row's own partiality and suggests
    # the 0.05 at which their bias is zero.
    rng = np.random.default_rng(7)
    n = 400
    zeta = rng.uniform(0.3, 1.0, n) * rng.choice([-1, 1], n)
    phi = rng.uniform(0.0, 0.5, n)
    used, truth = 0.08, 0.05
    width = np.radians(used) / np.abs(zeta)
    lo = phi - 3.0 * width
    hi = phi + rng.uniform(-1.0, 1.0, n) * width  # the box cut short: a scan's end
    p_model = eq.partiality(used, phi, zeta, lo, hi)
    keep = (p_model > 0.05) & (p_model < 0.99)
    reference = np.full(n, 1000.0)
    counts = reference * eq.partiality(truth, phi, zeta, lo, hi)
    s = eq.suggest_sigma_m(
        counts[keep],
        reference[keep],
        p_model[keep],
        phi[keep],
        zeta[keep],
        lo[keep],
        hi[keep],
    )
    assert abs(s.used - used) < 0.001, s
    assert s.suggested is not None and abs(s.suggested - truth) < 0.001, s
    assert s.bias_used > 0.0, s


def test_resolution_shells_are_quantiles_of_what_is_counted():
    keys, intensity, variance, clean, _ = planted()
    d = np.linspace(1.0, 5.0, len(keys))
    weak = np.arange(len(keys)) % 3 == 0
    variance = np.where(
        weak, (intensity / 2.0) ** 2, variance
    )  # I/sigma 2: never counted
    e = eq.Explanatory("d", d, quantiles=4, condition=clean)
    r = eq.analyse(keys, intensity, variance, clean, [e], min_i_sigma=5.0)
    assert len(r.tables[0].bins) == 4 and all(b.n > 0 for b in r.tables[0].bins)


def test_the_worst_are_written_with_their_reference_beside_them():
    from mxeq import refl

    n = 50
    table = refl.ReflectionTable(nrows=n)
    table.columns["miller_index"] = np.arange(3 * n, dtype=np.int64).reshape(n, 3)
    table.types["miller_index"] = "cctbx::miller::index<>"
    rel = np.linspace(-0.5, 0.45, n)  # no two the same size
    counted = np.ones(n, dtype=bool)
    counted[0] = False  # the furthest of all, but not counted
    result = eq.Result(
        "prf",
        n,
        1,
        eq.Bin("all", n, 0, 0, 0, 0),
        reference=np.full(n, 10.0),
        rel=rel,
        counted_rows=counted,
    )
    prepared = eq.Prepared(table, rows={"prf": np.arange(n)})
    out = eq.worst(prepared, result, "prf", 4)
    # The four counted furthest, in table order: rows 1, 2 and 3 (-0.48, -0.46,
    # -0.44) and 49 (+0.45).
    assert out.nrows == 4
    assert np.allclose(out.columns["equivalents.difference"], rel[[1, 2, 3, 49]])
    assert np.array_equal(
        out.columns["miller_index"], table.columns["miller_index"][[1, 2, 3, 49]]
    )
    assert np.allclose(out.columns["equivalents.reference"], 10.0)
