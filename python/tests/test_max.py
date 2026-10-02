"""mxi_max finds the largest valid count, and mxi_find --gpu-force narrows
32-bit frames to 16 bits losing nothing -- or refuses.

Metal's threshold takes 16 bits only. A 32-bit series whose counts never reach
0xFFFD -- the 16-bit markers above -- can be narrowed: the bad-pixel marker to
the 16-bit one, every other count as it is. Planted series: 16 bits; 32 bits
with a tile join's column (0xffffffff) and a dead pixel (0xfffffffe); and 32
bits with one real count of 70000. mxi_max must find the largest count apart from the marker and say
whether it fits; --gpu-force must give the 32-bit threshold's spots byte for
byte, and stop at the count that does not fit. Needs MXI_FIND and MXI_MAX.
"""

import os
import subprocess
import sys

import pytest

FIND = os.environ.get("MXI_FIND")
MAX = os.environ.get("MXI_MAX")
needs = pytest.mark.skipif(not (FIND and MAX), reason="set MXI_FIND and MXI_MAX")


def series(tmp_path, name, *extra):
    pytest.importorskip("h5py")
    pytest.importorskip("hdf5plugin")
    out = tmp_path / name
    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(out), "12", *extra],
        check=True,
        capture_output=True,
    )
    return str(out / "series.nxs")


def run(*args):
    return subprocess.run(list(args), capture_output=True, text=True)


@needs
def test_the_largest_valid_count_and_whether_it_fits(tmp_path):
    r16 = run(MAX, series(tmp_path, "s16"), "-j", "2")
    r32 = run(MAX, series(tmp_path, "s32", "--bits", "32"), "-j", "2")
    hot = run(MAX, series(tmp_path, "hot", "--bits", "32", "--hot"), "-j", "2")
    assert r16.returncode == 0 and "16-bit already" in r16.stdout
    assert r32.returncode == 0 and "fits in 16 bits" in r32.stdout
    largest = (
        [line for line in r16.stdout.splitlines() if "largest count" in line][0]
        .split(":")[1]
        .split()[0]
    )
    assert (
        f"largest count on a valid pixel: {largest} " in r32.stdout
    )  # the marker is not a count
    # A dead pixel on 12 frames, and 12 frames of a 256 pixel tile join: each
    # marker counted apart, neither a count.
    assert "pixels marked bad: 12; tile joins: 3072" in r32.stdout
    assert hot.returncode == 1 and "70000" in hot.stdout and "on image 7" in hot.stdout


@needs
def test_gpu_force_gives_the_32_bit_thresholds_spots_or_refuses(tmp_path):
    s32 = series(tmp_path, "s32", "--bits", "32")
    wide = run(FIND, s32, "-j", "2", "-o", str(tmp_path / "wide.refl"))
    narrow = run(
        FIND, s32, "-j", "2", "--gpu-force", "-o", str(tmp_path / "narrow.refl")
    )
    assert wide.returncode == 0 and narrow.returncode == 0, narrow.stderr
    assert (tmp_path / "wide.refl").read_bytes() == (
        tmp_path / "narrow.refl"
    ).read_bytes()
    hot = series(tmp_path, "hot", "--bits", "32", "--hot")
    refused = run(FIND, hot, "-j", "2", "--gpu-force", "-o", str(tmp_path / "hot.refl"))
    assert refused.returncode != 0
    assert "count of 70000" in refused.stderr + refused.stdout
    assert run(FIND, hot, "-j", "2", "-o", str(tmp_path / "hot2.refl")).returncode == 0
