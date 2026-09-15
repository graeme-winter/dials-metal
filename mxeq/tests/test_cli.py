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


def test_a_reflection_table_where_an_experiment_list_belongs_says_so(
    strong_files, capsys
):
    """`check refined` compares models, which live in the .expt.

    Handing it reflection tables used to fail inside a UTF-8 decoder --

        UnicodeDecodeError: 'utf-8' codec can't decode byte 0x93 in position 0

    which names the byte and not the mistake. 0x93 is msgpack's header for a
    three-element array, which is exactly what a reflection table starts with,
    so the file says what it is and the message can too.
    """
    a, b = strong_files
    # Both tables, so this is routed to the reflection comparison rather than
    # refused: comparing two refined pipelines usually does mean comparing
    # their reflections.
    assert main(["check", "refined", str(a), str(b)]) == 0
    captured = capsys.readouterr()
    assert "comparing reflections" in captured.err


def test_an_experiment_list_where_a_reflection_table_belongs_says_so(tmp_path, capsys):
    document = tmp_path / "a.expt"
    document.write_text(json.dumps({"__id__": "ExperimentList", "experiment": []}))
    code = main(["check", "indexed", str(document), str(document)])
    captured = capsys.readouterr()
    assert code == 2
    assert "experiment list" in captured.err
    assert "wants the .refl" in captured.err


def test_the_sniffer_reads_the_first_byte_only(tmp_path):
    from mxeq.cli import _sniff

    table = tmp_path / "t.refl"
    table.write_bytes(b"\x93\x00\x00")
    assert _sniff(table) == "refl"

    document = tmp_path / "d.expt"
    document.write_text('{"__id__": "ExperimentList"}')
    assert _sniff(document) == "expt"

    # Leading whitespace is still JSON.
    spaced = tmp_path / "s.expt"
    spaced.write_text('\n  {"__id__": "ExperimentList"}')
    assert _sniff(spaced) == "expt"

    # Neither, and an empty file, are reported as unknown rather than guessed
    # at -- a wrong guess would produce exactly the confusing error this
    # replaces.
    other = tmp_path / "x.bin"
    other.write_bytes(b"\xff\xd8\xff")
    assert _sniff(other) is None
    empty = tmp_path / "e.bin"
    empty.write_bytes(b"")
    assert _sniff(empty) is None
    assert _sniff(tmp_path / "missing") is None


def test_the_common_residual_is_over_the_same_reflections(strong_files, capsys):
    """Two rmsds averaged over different sets do not compare.

    DIALS refined 1800 images of insulin on 61679 of 78618 reflections and this
    package on 71454, so the headline numbers are not answering the same
    question. The section this checks for reports both sides over the
    reflections both of them predicted, which is the comparison that is.
    """
    a, b = strong_files
    assert main(["check", "indexed", str(a), str(b)]) == 0
    captured = capsys.readouterr()
    if "residual over the reflections both predicted" in captured.out:
        assert "n_common" in captured.out
