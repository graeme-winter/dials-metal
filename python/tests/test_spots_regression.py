"""The spot finder's regression script, run so that it cannot rot.

tests/spots/regression.sh checks the whole path from HDF5 to a reflection table
on planted data -- every planted reflection one spot where it was planted, the
summary accounting for every spot, 3D grouping, a missing frame, one thread and
eight byte for byte. Nothing ran it, and when the Python tools moved into mxeq
and mxi_find's report moved to standard output it failed four of its ten
checks unnoticed. Run here beside the suites, which have run ctest already.

Needs MXI_FIND, the build's mxi_find, and h5py, hdf5plugin and numpy for the
planted data.
"""

import os
import pathlib
import subprocess

import pytest

FIND = os.environ.get("MXI_FIND")
SCRIPT = (
    pathlib.Path(__file__).resolve().parents[2] / "tests" / "spots" / "regression.sh"
)


@pytest.mark.skipif(not FIND, reason="set MXI_FIND")
def test_the_spot_finder_regression_script_passes():
    pytest.importorskip("h5py")
    pytest.importorskip("hdf5plugin")
    env = dict(os.environ, REGRESSION_NO_CTEST="1")
    run = subprocess.run(
        ["bash", str(SCRIPT), str(pathlib.Path(FIND).resolve().parent)],
        capture_output=True,
        text=True,
        env=env,
        timeout=600,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    assert "0 failed" in run.stdout, run.stdout
