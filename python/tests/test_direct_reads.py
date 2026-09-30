"""Frames read by chunk offset, outside HDF5's lock, are the frames HDF5 gives.

The NXmx reader finds every chunk's place in its file once, and then reads a
frame with pread rather than H5Dread_chunk, so that threads do not queue for
HDF5's one lock. The bytes must be the same: here the spot finder runs on the
planted fixtures -- a virtual dataset across two data files, and one with a
frame never written -- both ways -- by offset with MXI_HDF5_DIRECT=1, and through HDF5 by default -- and
the tables must be identical, the report saying which path the frames took.

Needs MXI_FIND, and h5py, hdf5plugin and numpy for the fixtures.
"""

import os
import re
import subprocess
import sys

import pytest

FIND = os.environ.get("MXI_FIND")


def fixture(tmp_path, name, *extra):
    out = tmp_path / name
    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(out), "12", *extra],
        check=True,
        capture_output=True,
    )
    return out / "series.nxs"


def find(tmp_path, master, name, direct):
    env = dict(os.environ)
    env.pop("MXI_HDF5_DIRECT", None)
    if direct:
        env["MXI_HDF5_DIRECT"] = "1"
    run = subprocess.run(
        [
            FIND,
            "-j",
            "3",
            "--timing",
            "-o",
            str(tmp_path / (name + ".refl")),
            str(master),
        ],
        capture_output=True,
        text=True,
        env=env,
    )
    assert run.returncode == 0, run.stderr
    counts = re.search(r"outside HDF5's lock: (\d+); through HDF5: (\d+)", run.stdout)
    assert counts, run.stdout
    return (tmp_path / (name + ".refl")).read_bytes(), int(counts[1]), int(counts[2])


@pytest.mark.skipif(not FIND, reason="set MXI_FIND")
@pytest.mark.parametrize(
    "extra", [(), ("--holey",)], ids=["virtual", "a frame never written"]
)
def test_direct_reads_are_the_frames_hdf5_gives(tmp_path, extra):
    pytest.importorskip("h5py")
    pytest.importorskip("hdf5plugin")
    master = fixture(tmp_path, "series", *extra)
    direct, by_offset, _ = find(tmp_path, master, "direct", True)
    through, by_offset_forced, by_hdf5 = find(tmp_path, master, "hdf5", False)
    assert direct == through, "the tables differ"
    assert by_offset > 0, "no frame was read by its chunk offset"
    assert (
        by_offset_forced == 0 and by_hdf5 > 0
    ), "MXI_HDF5_DIRECT=0 still read by offset"
