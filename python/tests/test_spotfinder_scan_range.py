"""The spot finder reads the scan's frames, and gives them the scan's z.

Two layouts, because they need opposite offsets and the old code handled one:

* a master file for the whole run, the scan covering images 6 to 10 -- which
  arrive as frames 5 to 9, already the right z;
* a file holding only the slice, the same images arriving as frames 0 to 4,
  which need five adding.

`single_file_indices` in the .expt is what tells them apart.

Found by spot finding 1800 images of a 36000 image run, which warned that z
would be wrong outside the scan and then read all 36000.  Restricting it to the
scan then showed the z offset had been added twice.
"""

import json
import os
import shutil
import subprocess

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
hdf5plugin = pytest.importorskip("hdf5plugin")

from mxeq import refl

BINARY = os.environ.get("MXI_FIND") or shutil.which("mxi_find")
pytestmark = pytest.mark.skipif(
    not BINARY, reason="mxi_find is not built or not on PATH"
)
TEMPLATE_EXPT = os.environ.get("MXI_TEMPLATE_EXPT")


def series(directory, frames, planted):
    """A series of 64x64 frames with one spot on each, where `planted` says
    which rows they sit at so no two merge."""
    data = directory / "data_000001.h5"
    master = directory / "master.nxs"
    rng = np.random.default_rng(0)
    with h5py.File(data, "w") as f:
        d = f.create_dataset(
            "data",
            shape=(frames, 64, 64),
            dtype="u2",
            chunks=(1, 64, 64),
            **hdf5plugin.Bitshuffle(nelems=0, cname="lz4"),
        )
        for i in range(frames):
            frame = rng.poisson(0.1, (64, 64)).astype("u2")
            y = planted(i)
            frame[y : y + 3, 20:23] += 400
            d[i] = frame
    with h5py.File(master, "w") as f:
        layout = h5py.VirtualLayout(shape=(frames, 64, 64), dtype="u2")
        layout[...] = h5py.VirtualSource(data.name, "data", shape=(frames, 64, 64))
        g = f.create_group("entry")
        g.attrs["NX_class"] = "NXentry"
        dg = g.create_group("data")
        dg.attrs["NX_class"] = "NXdata"
        dg.create_virtual_dataset("data", layout)
    return master


def experiment(directory, template, master, indices):
    e = json.load(open(template))
    panel = e["detector"][0]["panels"][0]
    panel["image_size"] = [64, 64]
    scan = e["scan"][0]
    scan["image_range"] = [6, 10]
    scan["properties"]["oscillation"] = [0.1 * i for i in range(5)]
    for k in list(scan["properties"]):
        if k != "oscillation" and isinstance(scan["properties"][k], list):
            scan["properties"][k] = scan["properties"][k][:5]
    e["imageset"][0]["template"] = str(master)
    e["imageset"][0]["single_file_indices"] = list(indices)
    e["crystal"] = []
    e["experiment"][0].pop("crystal", None)
    path = directory / "scan.expt"
    json.dump(e, open(path, "w"))
    return path


def find(expt, directory):
    out = directory / "strong.refl"
    run = subprocess.run(
        [BINARY, "-e", str(expt), "-o", str(out)], capture_output=True, text=True
    )
    assert run.returncode == 0, run.stderr
    # The report is on standard output; standard error is for errors, and is
    # what the assertion above shows if the run fails.
    return refl.load(str(out)), run.stdout


@pytest.mark.skipif(not TEMPLATE_EXPT, reason="needs MXI_TEMPLATE_EXPT")
def test_a_master_for_the_whole_run_reads_only_the_scan(tmp_path):
    # Twenty frames, the scan covering five of them. Eight rows apart on
    # successive frames: a spot three rows tall on consecutive frames three
    # rows apart touches its neighbour in three dimensions and the two merge,
    # which is what this test's first version did.
    master = series(tmp_path, 20, lambda i: 4 + 8 * (i % 7))
    expt = experiment(tmp_path, TEMPLATE_EXPT, master, range(5, 10))
    table, said = find(expt, tmp_path)
    assert "Reading 5 of the 20 images" in said
    z = table.columns["xyzobs.px.value"][:, 2]
    assert table.nrows == 5
    assert ((z >= 5) & (z < 10)).all(), z


@pytest.mark.skipif(not TEMPLATE_EXPT, reason="needs MXI_TEMPLATE_EXPT")
def test_a_file_holding_only_the_slice_still_gets_the_scan_z(tmp_path):
    # The same five images, as frames 0 to 4 of a file of their own.
    master = series(tmp_path, 5, lambda i: 6 + 8 * i)
    expt = experiment(tmp_path, TEMPLATE_EXPT, master, range(0, 5))
    table, _ = find(expt, tmp_path)
    z = table.columns["xyzobs.px.value"][:, 2]
    assert table.nrows == 5
    assert ((z >= 5) & (z < 10)).all(), z
