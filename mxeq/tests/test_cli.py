"""The command line, end to end over files on disk."""

from __future__ import annotations

import json

import numpy as np
import pytest

import fixtures
from mxeq import refl
from mxeq.cli import main


@pytest.fixture
def strong_files(tmp_path):
    a = fixtures.strong_table(n=300)
    b = fixtures.drop(a, fraction=0.05, seed=12)
    refl.write(tmp_path / "a.refl", a)
    refl.write(tmp_path / "b.refl", b)
    return str(tmp_path / "a.refl"), str(tmp_path / "b.refl")


@pytest.fixture
def integrated_files(tmp_path):
    a = fixtures.integrated_table(n=500)
    refl.write(tmp_path / "ia.refl", a)
    refl.write(tmp_path / "ib.refl", a)
    path = fixtures.write_experiments(tmp_path / "i.expt")
    return str(tmp_path / "ia.refl"), str(tmp_path / "ib.refl"), path


def test_check_strong(strong_files, capsys):
    a, b = strong_files
    assert main(["check", "strong", a, b]) == 0
    out = capsys.readouterr().out
    assert "mxeq strong" in out
    assert "spatial match" in out


def test_check_emits_json(strong_files, capsys):
    a, b = strong_files
    assert main(["check", "strong", a, b, "--json"]) == 0
    decoded = json.loads(capsys.readouterr().out)
    assert decoded["sections"]["spatial match"]["n_matched"] > 0


def test_check_integrated_with_experiments(integrated_files, capsys):
    a, b, e = integrated_files
    assert main(["check", "integrated", a, b, "-e", e, "--bins", "5"]) == 0
    assert "intensity.sum.value" in capsys.readouterr().out


def test_check_refined(tmp_path, capsys):
    a = fixtures.write_experiments(tmp_path / "a.expt")
    b = fixtures.write_experiments(tmp_path / "b.expt", cell=78.2)
    assert main(["check", "refined", a, b]) == 0
    assert "unit cell" in capsys.readouterr().out


def test_auto_guesses_the_boundary(integrated_files, capsys):
    a, b, _ = integrated_files
    assert main(["check", "auto", a, b]) == 0
    assert "mxeq scaled" in capsys.readouterr().out


def test_operator_is_parsed_and_applied(tmp_path, capsys):
    operator = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.int64)
    a = fixtures.integrated_table(n=300)
    b = fixtures.reindex(a, operator)
    refl.write(tmp_path / "a.refl", a)
    refl.write(tmp_path / "b.refl", b)
    argv = ["check", "integrated", str(tmp_path / "a.refl"), str(tmp_path / "b.refl")]
    main(argv + ["--json"])
    without = json.loads(capsys.readouterr().out)
    main(argv + ["--json", "--operator", "0,1,0,0,0,1,1,0,0"])
    with_operator = json.loads(capsys.readouterr().out)
    assert without["sections"]["keyed match"]["fraction_matched"] < 0.05
    assert with_operator["sections"]["keyed match"]["fraction_matched"] == 1.0


def test_operator_of_the_wrong_length_is_rejected(strong_files):
    a, b = strong_files
    with pytest.raises(SystemExit):
        main(["check", "integrated", a, b, "--operator", "1,0,0"])


def test_inspect_a_reflection_table(strong_files, capsys):
    a, _ = strong_files
    assert main(["inspect", a]) == 0
    out = capsys.readouterr().out
    assert "reflection table" in out
    assert "xyzobs.px.value" in out
    assert "vec3<double>" in out


def test_inspect_an_experiment_list(tmp_path, capsys):
    path = fixtures.write_experiments(tmp_path / "a.expt")
    assert main(["inspect", path]) == 0
    out = capsys.readouterr().out
    assert "experiment list" in out
    assert "cell" in out


def test_inspect_falls_back_to_a_structure_dump(tmp_path, capsys):
    """The escape hatch: a file that is msgpack but not what was expected."""
    import msgpack

    path = tmp_path / "odd.refl"
    path.write_bytes(msgpack.packb({"something": "else"}))
    assert main(["inspect", str(path)]) == 2
    assert "something" in capsys.readouterr().out
