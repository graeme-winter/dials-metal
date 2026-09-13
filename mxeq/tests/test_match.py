"""The joins, with hand-worked cases."""

from __future__ import annotations

import numpy as np
import pytest

import fixtures
from mxeq import match


def test_identical_positions_all_match():
    xyz = fixtures.strong_table(n=50)["xyzobs.px.value"]
    m = match.match_positions(xyz, xyz, radius=1.0)
    assert m.n_matched == 50
    assert m.fraction_matched == 1.0
    assert np.allclose(m.distance, 0.0)


def test_matching_is_mutual_so_one_row_cannot_claim_two():
    # Two rows of A sit on top of one row of B. Only one pairing is possible,
    # and the other row of A must be reported unmatched rather than sharing.
    a = np.array([[10.0, 10.0, 0.0], [10.2, 10.0, 0.0]])
    b = np.array([[10.0, 10.0, 0.0]])
    m = match.match_positions(a, b, radius=1.0)
    assert m.n_matched == 1
    assert m.index_a.tolist() == [0]
    assert m.n_only_a == 1


def test_nothing_beyond_the_radius_matches():
    a = np.array([[0.0, 0.0, 0.0]])
    b = np.array([[5.0, 0.0, 0.0]])
    assert match.match_positions(a, b, radius=2.0).n_matched == 0
    assert match.match_positions(a, b, radius=6.0).n_matched == 1


def test_fraction_matched_uses_the_larger_table():
    a = np.zeros((1, 3))
    b = np.vstack([np.zeros((1, 3)), np.arange(3, 30).reshape(9, 3) * 100.0])
    m = match.match_positions(a, b, radius=1.0)
    assert m.n_matched == 1
    assert m.fraction_matched == pytest.approx(0.1)


def test_z_scale_separates_the_scan_axis():
    # One image apart in z and nothing else. With z counted as a pixel it is
    # inside a radius of 2; scaled so that one image is ten pixels, it is not.
    a = np.array([[0.0, 0.0, 0.0]])
    b = np.array([[0.0, 0.0, 1.0]])
    assert match.match_positions(a, b, radius=2.0, scale=(1, 1, 1)).n_matched == 1
    assert match.match_positions(a, b, radius=2.0, scale=(1, 1, 0.1)).n_matched == 0


def test_empty_input_is_not_an_error():
    m = match.match_positions(np.empty((0, 3)), np.zeros((3, 3)), radius=1.0)
    assert m.n_matched == 0 and m.n_only_b == 3


def test_keyed_join_pairs_on_the_key_not_the_order():
    hkl = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    order = [2, 0, 1]
    m, duplicates = match.match_keys([hkl], [hkl[order]])
    assert duplicates == 0
    assert m.n_matched == 3
    for i, j in zip(m.index_a, m.index_b):
        assert hkl[i].tolist() == hkl[order][j].tolist()


def test_duplicate_keys_are_paired_in_tie_break_order():
    # The same reflection recorded twice, as on a scan that turns more than
    # once. Frame order must decide which pairs with which.
    hkl = np.array([[1, 0, 0], [1, 0, 0]])
    frames_a = np.array([10.0, 400.0])
    frames_b = np.array([402.0, 12.0])  # deliberately the other way round
    m, duplicates = match.match_keys([hkl], [hkl], frames_a, frames_b)
    # One *key* is duplicated, carrying two rows on each side.
    assert duplicates == 1
    assert m.n_matched == 2
    pairs = {(int(i), int(j)) for i, j in zip(m.index_a, m.index_b)}
    assert pairs == {(0, 1), (1, 0)}


def test_keys_present_on_one_side_only_are_counted_not_matched():
    a = np.array([[1, 0, 0], [2, 0, 0]])
    b = np.array([[2, 0, 0], [3, 0, 0]])
    m, _ = match.match_keys([a], [b])
    assert m.n_matched == 1 and m.n_only_a == 1 and m.n_only_b == 1


def test_multiple_key_columns_all_take_part():
    hkl = np.array([[1, 0, 0], [1, 0, 0]])
    entering = np.array([0, 1])
    m, _ = match.match_keys([hkl, entering], [hkl, entering])
    assert m.n_matched == 2
    m, _ = match.match_keys([hkl, entering], [hkl, 1 - entering])
    # Still two pairs, but crossed over: entering is part of the identity.
    assert sorted(zip(m.index_a.tolist(), m.index_b.tolist())) == [(0, 1), (1, 0)]
