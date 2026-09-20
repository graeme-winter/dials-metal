"""The on-disk format, pinned against DIALS-written files.

These tests exist because the first version of the reader was wrong about the
layout in two ways -- a two-element top level instead of three, and a bare
payload blob instead of a (count, blob) pair -- and the full suite passed
anyway.  It passed because the only thing checking the format was a round-trip
against this package's own writer, and a writer built from the same wrong
assumption agrees with it perfectly.

A round-trip proves self-consistency. Only a file someone else wrote proves
the format. So the assertions below are split in two: the ones about
``tests/data`` are evidence, and the ones constructed by hand are pins that
fail loudly if the reader is loosened.
"""

from __future__ import annotations

import json
import pathlib

import msgpack
import numpy as np
import pytest

from mxeq import refl
from mxeq.checks import strong

DATA = pathlib.Path(__file__).parent / "data"

#: Cut from `dials.find_spots` and `dials-metal-find-spots` output on the same
#: 300-image sweep by `tests/make_real_fixture.py`; see its docstring for why
#: the structure is sliced rather than re-written.
DIALS_FILE = DATA / "dials_strong.refl"
OTHER_FILE = DATA / "other_strong.refl"


def test_real_file_header_is_a_three_element_array():
    raw = DIALS_FILE.read_bytes()
    # 0x93 is a msgpack fixarray of three. This is the byte that was wrong.
    assert raw[0] == 0x93
    assert raw[2:29] == refl.TAG.encode()


def test_real_file_reads():
    table = refl.load(str(DIALS_FILE))
    assert table.nrows == 48
    assert table.version == 2
    assert len(table.identifiers) == 1
    assert "xyzobs.px.value" in table
    assert table["xyzobs.px.value"].shape == (48, 3)
    assert table["bbox"].shape == (48, 6)


def test_real_file_element_widths():
    """The dtype widths, checked against a file rather than assumed."""
    raw = DIALS_FILE.read_bytes()
    data = msgpack.unpackb(raw, raw=True, strict_map_key=False)[-1][b"data"]
    expected = {
        b"int": 4,
        b"std::size_t": 8,
        b"double": 8,
        b"vec3<double>": 24,
        b"int6": 24,
    }
    for name, (type_name, payload) in data.items():
        count, blob = payload
        assert count == 48, name
        assert len(blob) == 48 * expected[type_name], name


def test_real_file_values_are_physically_sensible():
    """A reader can decode a file and still have the byte order wrong."""
    table = refl.load(str(DIALS_FILE))
    xyz = table["xyzobs.px.value"]
    bbox = table["bbox"]
    # Centroids inside the detector, and inside their own bounding boxes.
    assert (xyz[:, 0] > 0).all() and (xyz[:, 0] < 5000).all()
    assert (xyz[:, 1] > 0).all() and (xyz[:, 1] < 5000).all()
    assert (xyz[:, 0] >= bbox[:, 0]).all() and (xyz[:, 0] <= bbox[:, 1]).all()
    assert (xyz[:, 1] >= bbox[:, 2]).all() and (xyz[:, 1] <= bbox[:, 3]).all()
    # Strong spots carry at least min_count pixels and a positive intensity.
    assert (table["n_signal"] >= 2).all()
    assert (table["intensity.sum.value"] > 0).all()
    # The strong flag, bit five.
    assert (table["flags"] & (1 << 5)).all()


def test_real_pair_is_equivalent_and_reordered():
    """The finding from the full files, reduced to something that fits in git.

    DIALS and the Metal spot finder produce the same spots with bitwise
    identical intensities and bounding boxes, in a different row order, with
    centroids differing only by floating-point summation order.
    """
    a = refl.load(str(DIALS_FILE))
    b = refl.load(str(OTHER_FILE))
    report = strong.check(a, b)
    spatial = next(s for s in report.sections if s.title == "spatial match")
    assert spatial.values["n_matched"] == 48
    assert spatial.values["fraction_matched"] == 1.0
    assert spatial.values["separation"]["max"] < 1e-9

    shared = next(
        s for s in report.sections if s.title == "shared columns over matched pairs"
    )
    agreement = dict(
        (row[0], row[1]) for row in shared.values["column_agreement"]["rows"]
    )
    for column in ("bbox", "flags", "n_signal", "intensity.sum.value"):
        assert agreement[column] == "yes", column
    # The centroid is the one thing that is not bitwise equal.
    assert agreement["xyzobs.px.value"] == "no"


def test_shoebox_is_carried_as_opaque_not_decoded():
    """A shoebox must not be decoded, and must not stop the file reading."""
    table = refl.load(str(DIALS_FILE))
    # The fixture has it stripped, so construct the case rather than skip it.
    raw = DIALS_FILE.read_bytes()
    doc = msgpack.unpackb(raw, raw=True, strict_map_key=False)
    doc[-1][b"data"][b"shoebox"] = [b"Shoebox<>", [48, b"\x00" * 16]]
    rebuilt = refl.loads(msgpack.packb(doc, use_bin_type=True))
    assert rebuilt.opaque == {"shoebox": "Shoebox<>"}
    assert set(rebuilt.columns) == set(table.columns)


# --------------------------------------------------------------------------
# pins: these fail if the reader is loosened
# --------------------------------------------------------------------------


def _document(body: dict, version: int | None = 2) -> bytes:
    parts = [refl.TAG] + ([version] if version is not None else []) + [body]
    return msgpack.packb(parts, use_bin_type=True)


def test_writer_emits_the_three_element_form_with_wrapped_payloads():
    """Pin the writer's output shape, not just its round-trip behaviour."""
    table = refl.ReflectionTable(nrows=2)
    table.columns["id"] = np.array([0, 0], dtype=np.int32)
    table.types["id"] = "int"
    doc = msgpack.unpackb(refl.dumps(table), raw=True, strict_map_key=False)
    assert len(doc) == 3
    assert doc[1] == refl.FORMAT_VERSION
    type_name, payload = doc[-1][b"data"][b"id"]
    assert type_name == b"int"
    assert payload[0] == 2
    assert isinstance(payload[1], bytes) and len(payload[1]) == 8


def test_bare_payload_without_a_count_is_still_accepted():
    body = {
        "nrows": 2,
        "identifiers": {},
        "data": {"id": ["int", np.array([1, 2], dtype="<i4").tobytes()]},
    }
    assert refl.loads(_document(body))["id"].tolist() == [1, 2]


def test_two_element_form_without_a_version_is_still_accepted():
    body = {
        "nrows": 2,
        "identifiers": {},
        "data": {"id": ["int", [2, np.array([3, 4], dtype="<i4").tobytes()]]},
    }
    table = refl.loads(_document(body, version=None))
    assert table["id"].tolist() == [3, 4]
    assert table.version is None


def test_a_payload_count_disagreeing_with_nrows_is_rejected():
    body = {
        "nrows": 4,
        "identifiers": {},
        "data": {"id": ["int", [2, np.array([1, 2], dtype="<i4").tobytes()]]},
    }
    with pytest.raises(refl.ReflFormatError, match="declares 2 elements"):
        refl.loads(_document(body))


def test_structure_survives_a_binary_payload():
    """The escape hatch's own bug: it decoded every blob as UTF-8.

    A column of packed doubles is not valid UTF-8, so the one tool meant to
    explain an unreadable file raised on the first real file it was given.
    """
    body = {
        "nrows": 1,
        "data": {"x": ["double", [1, np.array([1.5], dtype="<f8").tobytes()]]},
    }
    lines = refl.structure(_document(body))
    assert any("binary, 8 bytes" in line for line in lines)
    assert any("'x'" in line for line in lines)


def test_structure_survives_a_non_utf8_map_key():
    lines = refl.structure(msgpack.packb({b"\x9d\x00bad": 1}, use_bin_type=False))
    assert lines
    assert any("bytes" in line for line in lines)


def test_structure_describes_the_real_file():
    lines = refl.structure(DIALS_FILE.read_bytes(), max_items=6)
    assert any("array, 3 items" in line for line in lines)
    assert any("dials::af::reflection_table" in line for line in lines)


# --------------------------------------------------------------------------
# .expt shapes that changed under us
# --------------------------------------------------------------------------

REAL_SCAN = DATA / "real_scan.json"


def _expt_with(scan: dict) -> str:
    import fixtures

    doc = fixtures.experiments_dict()
    doc["scan"] = [scan]
    return json.dumps(doc)


def test_scan_reads_the_per_image_oscillation_array():
    """Current DIALS writes properties.oscillation, not [start, width].

    Reading only the old key against a current file gives (0.0, 0.0) for both
    sides of a comparison, which then agree perfectly and say nothing. A silent
    pass is worse than a failure, so this is pinned against a real scan block.
    """
    from mxeq import expt

    real = json.loads(REAL_SCAN.read_text())["scan"]
    scan = expt.loads(_expt_with(real))[0].scan
    assert scan.image_range == (1, 6)
    assert scan.oscillation[0] == pytest.approx(0.0)
    assert scan.oscillation[1] == pytest.approx(0.1, abs=1e-12)
    # Accumulation round-off only.
    assert scan.max_width_deviation < 1e-10


def test_scan_still_reads_the_old_oscillation_pair():
    from mxeq import expt

    scan = expt.loads(_expt_with({"image_range": [1, 100], "oscillation": [5.0, 0.2]}))[
        0
    ].scan
    assert scan.oscillation == (5.0, 0.2)


def test_scan_width_comes_from_the_endpoints():
    from mxeq import expt

    values = [0.1 * i for i in range(300)]
    values[1] += 0.05  # one bad element near the start
    scan = expt.loads(
        _expt_with({"image_range": [1, 300], "properties": {"oscillation": values}})
    )[0].scan
    assert scan.oscillation[1] == pytest.approx(0.1, abs=1e-12)
    # ...and the perturbation is reported, not swallowed.
    assert scan.max_width_deviation > 0.1


def test_millimetre_columns_are_not_used_as_a_position_key():
    """xyzobs.mm is import-time and goes stale; it must never be a join key.

    It carries the inverse parallax correction under whatever geometry was
    current at import and is never recomputed, so after refinement it describes
    a detector that no longer exists. Two pipelines that refined to slightly
    different detectors would be compared in two different millimetre frames,
    and the difference would read as a centroid disagreement.
    """
    import numpy as np

    from mxeq.checks.common import position_column

    table = refl.ReflectionTable(nrows=2)
    table.columns["xyzobs.mm.value"] = np.zeros((2, 3))
    table.types["xyzobs.mm.value"] = "vec3<double>"
    with pytest.raises(KeyError, match="stale"):
        position_column(table)

    table.columns["xyzobs.px.value"] = np.zeros((2, 3))
    table.types["xyzobs.px.value"] = "vec3<double>"
    assert position_column(table) == "xyzobs.px.value"
