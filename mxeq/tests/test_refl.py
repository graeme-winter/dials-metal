"""The reflection table reader, and what counts as a prediction."""

from __future__ import annotations

import numpy as np

from mxeq.refl import ReflectionTable, has_prediction


def test_a_row_with_no_prediction_is_recognised_however_it_was_left():
    """Two ways a row carries no prediction, and both occur in real files.

    Left as written by whoever allocated the column: exactly zero, which is
    what this package writes and what a failed prediction leaves behind. Or
    never written at all, in which case it holds whatever was in the memory --
    on a real DIALS indexed.refl those rows read [1.5e-320, 5.2e-310, 0.0],
    denormals, which compare unequal to zero.

    Testing only against zero missed the second kind and reported prediction
    offsets of fifteen hundred pixels between two pipelines that agree.
    """
    table = ReflectionTable(nrows=4)
    table.columns["miller_index"] = np.array(
        [[1, 2, 3], [4, 5, 6], [0, 0, 0], [7, 8, 9]], dtype=np.int64
    )
    table.columns["xyzcal.px"] = np.array(
        [
            [100.0, 200.0, 3.0],  # a real prediction
            [1.5e-320, 5.2e-310, 0.0],  # never written
            [0.0, 0.0, 0.0],  # not indexed
            [0.0, 0.0, 0.0],  # indexed, prediction failed
        ]
    )
    got = has_prediction(table)
    assert list(got) == [True, False, False, False]


def test_has_prediction_can_be_taken_at_a_set_of_rows():
    table = ReflectionTable(nrows=3)
    table.columns["miller_index"] = np.array(
        [[1, 0, 0], [0, 0, 0], [0, 1, 0]], dtype=np.int64
    )
    table.columns["xyzcal.px"] = np.array(
        [[1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [2.0, 2.0, 2.0]]
    )
    assert list(has_prediction(table, [2, 0])) == [True, True]
