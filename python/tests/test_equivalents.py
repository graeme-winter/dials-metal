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
