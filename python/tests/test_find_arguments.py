"""mxi_find's argument: the images, or an experiment list naming them.

As dials.find_spots takes it, `mxi_find imported.expt` finds spots on the images
the experiment list's imageset names, as `-e imported.expt` does; given it in
the images' place, mxi_find once handed the JSON to HDF5, which printed its
diagnostic stack before the one useful line. A file that is not HDF5 is now said
to be so plainly, and two experiment lists are refused. Needs MXI_FIND.
"""

import os
import subprocess

import pytest

FIND = os.environ.get("MXI_FIND")
needs_find = pytest.mark.skipif(not FIND, reason="set MXI_FIND")


def run(tmp_path, *args):
    return subprocess.run([FIND, *args], capture_output=True, text=True, cwd=tmp_path)


@needs_find
def test_a_file_that_is_not_hdf5_is_said_to_be_so_plainly(tmp_path):
    (tmp_path / "text.nxs").write_text("hello, not images\n")
    result = run(tmp_path, "text.nxs", "-o", "x.refl")
    assert result.returncode != 0
    both = result.stdout + result.stderr
    assert "is not an HDF5 file" in both
    assert "HDF5-DIAG" not in both


@needs_find
def test_json_in_the_images_place_is_taken_for_an_experiment_list(tmp_path):
    (tmp_path / "fake.nxs").write_text('{"not": "an experiment list"}\n')
    result = run(tmp_path, "fake.nxs", "-o", "x.refl")
    assert result.returncode != 0
    both = result.stdout + result.stderr
    assert "is not an experiment list" in both
    assert "HDF5-DIAG" not in both


@needs_find
@pytest.mark.parametrize(
    "order", [("a.expt", "-e", "b.expt"), ("-e", "a.expt", "b.expt")]
)
def test_two_experiment_lists_are_refused_either_way_round(tmp_path, order):
    for name in ("a.expt", "b.expt"):
        (tmp_path / name).write_text("{}\n")
    result = run(tmp_path, *order)
    assert result.returncode == 2
    assert "two experiment lists" in result.stderr
