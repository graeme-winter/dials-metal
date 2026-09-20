"""The readers, which are the part with the most assumptions in it."""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

import fixtures
from mxeq import expt, refl


def test_refl_round_trip_preserves_every_column():
    table = fixtures.integrated_table(n=64)
    back = refl.loads(refl.dumps(table))
    assert back.nrows == table.nrows
    assert set(back.columns) == set(table.columns)
    assert back.identifiers == table.identifiers
    for name, values in table.columns.items():
        if values.dtype == bool:
            assert np.array_equal(back[name], values), name
        else:
            assert np.allclose(
                np.asarray(back[name], dtype=float), np.asarray(values, dtype=float)
            ), name


def test_refl_compound_columns_keep_their_width():
    table = fixtures.integrated_table(n=16)
    back = refl.loads(refl.dumps(table))
    assert back["miller_index"].shape == (16, 3)
    assert back["xyzcal.px"].shape == (16, 3)
    assert back.types["miller_index"] == "cctbx::miller::index<>"


def test_refl_file_starts_with_the_dials_tag():
    # If this ever fails against a real file, the tag is the first thing to
    # check and the fixture writer is what is wrong.
    raw = refl.dumps(fixtures.strong_table(n=4))
    assert refl.TAG.encode() in raw[:64]


def test_refl_rejects_a_truncated_compound_payload():
    table = fixtures.integrated_table(n=8)
    doc = refl.dumps(table)
    broken = refl.loads(doc)
    # Chop one scalar off a vec3 payload and re-encode by hand.
    broken.columns["xyzcal.px"] = broken["xyzcal.px"].ravel()[:-1]
    with pytest.raises(refl.ReflFormatError):
        refl._decode_column(
            "xyzcal.px",
            "vec3<double>",
            broken.columns["xyzcal.px"].astype("<f8").tobytes(),
        )


def test_refl_refuses_a_document_that_is_not_a_table():
    import msgpack

    with pytest.raises(refl.ReflFormatError):
        refl.loads(msgpack.packb({"nrows": 3}))


def test_structure_describes_an_unreadable_document():
    import msgpack

    lines = refl.structure(msgpack.packb({"surprise": [1, 2, 3]}))
    assert any("surprise" in line for line in lines)
    assert any("array" in line for line in lines)


def test_structure_never_raises_on_rubbish():
    assert refl.structure(b"\xc1\xc1\xc1not msgpack")


def test_expt_reads_cell_and_scan(tmp_path):
    path = fixtures.write_experiments(tmp_path / "a.expt")
    experiments = expt.load(path)
    assert len(experiments) == 1
    crystal = experiments[0].crystal
    assert crystal is not None
    a, b, c, alpha, beta, gamma = crystal.cell
    assert a == pytest.approx(fixtures.CELL)
    assert alpha == pytest.approx(90.0)
    assert experiments[0].scan.num_images == 600
    assert experiments[0].beam.wavelength == pytest.approx(0.9537)


def test_setting_matrix_is_the_inverse_not_the_inverse_transpose():
    """The A matrix convention, on a cell where the two differ.

    A triclinic cell is used deliberately: for an orthogonal cell the inverse
    and the inverse transpose are both diagonal and the wrong one passes.
    """
    real = np.array([[10.0, 0.0, 0.0], [3.0, 20.0, 0.0], [1.0, 2.0, 30.0]])
    crystal = expt.Crystal(a=real)
    # By definition a . a* = 1 and a . b* = 0, i.e. real @ setting == I.
    assert np.allclose(real @ crystal.setting, np.eye(3), atol=1e-12)


def test_cell_from_real_space_matches_hand_calculation():
    real = np.array([[10.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 30.0]])
    a, b, c, alpha, beta, gamma = expt.cell_from_real_space(real)
    assert (a, b, c) == pytest.approx((10.0, 20.0, 30.0))
    assert (alpha, beta, gamma) == pytest.approx((90.0, 90.0, 90.0))

    # 60 degrees between a and b, chosen so the answer is exact.
    real = np.array([[1.0, 0.0, 0.0], [0.5, math.sqrt(3) / 2, 0.0], [0.0, 0.0, 1.0]])
    *_, gamma = expt.cell_from_real_space(real)
    assert gamma == pytest.approx(60.0)


def test_expt_handles_a_missing_model(tmp_path):
    doc = fixtures.experiments_dict()
    doc["experiment"][0]["crystal"] = -1
    path = tmp_path / "nocrystal.expt"
    path.write_text(json.dumps(doc))
    assert expt.load(str(path))[0].crystal is None


def test_expt_reads_scan_varying_cells(tmp_path):
    path = fixtures.write_experiments(tmp_path / "sv.expt", n_scan_points=20)
    crystal = expt.load(path)[0].crystal
    assert crystal.scan_varying
    cells = crystal.cells_at_scan_points()
    assert cells.shape == (20, 6)
    # Planted as a monotonic expansion, so the last must exceed the first.
    assert cells[-1, 0] > cells[0, 0]
