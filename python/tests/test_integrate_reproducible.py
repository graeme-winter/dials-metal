"""mxi_integrate gives the same bytes whatever the thread count.

Profile learning summed per-thread partials, and that was a data race: threads
are new on every parallel call but the calling thread takes part in each, so a
lane it chose once outlived the call, while each later call's new threads were
numbered from zero again -- two threads adding into the same partial profiles.
Even without the race the additions happened in an order the scheduler chose,
and two four-thread runs differed by 3.6e-11 in intensity.prf.value. Learning
now sums fixed blocks by reflection index, in order, and adds them in order.

This runs the real program on real data, so it needs:

    MXI_INTEGRATE      the mxi_integrate binary
    MXI_TEST_EXPT      a refined .expt
    MXI_TEST_REFL      its strong or refined .refl
    MXI_TEST_IMAGES    the images, if not found from the .expt

and is skipped, saying so, without them.
"""

import os
import subprocess

import pytest

BINARY = os.environ.get("MXI_INTEGRATE")
EXPT = os.environ.get("MXI_TEST_EXPT")
REFL = os.environ.get("MXI_TEST_REFL")
IMAGES = os.environ.get("MXI_TEST_IMAGES")
pytestmark = pytest.mark.skipif(
    not (BINARY and EXPT and REFL),
    reason="needs MXI_INTEGRATE, MXI_TEST_EXPT and MXI_TEST_REFL",
)


def integrate(tmp_path, threads, tag):
    out = tmp_path / f"integrated_{tag}.refl"
    command = [
        BINARY,
        EXPT,
        REFL,
        "-o",
        str(out),
        "--threads",
        str(threads),
        "--d-min",
        "2.5",
        "--scan-blocks",
        "2",
        "--first-image",
        "0",
        "--last-image",
        "40",
        # short chunks: many parallel reductions, which is where the
        # race lived and where an order chosen by the scheduler shows
        "--window",
        "2",
        "--max-boxes",
        "400",
    ]
    if IMAGES:
        command += ["--images", IMAGES]
    run = subprocess.run(command, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr[-2000:]
    return out.read_bytes()


def test_the_thread_count_does_not_change_a_byte(tmp_path):
    one = integrate(tmp_path, 1, "1")
    four = integrate(tmp_path, 4, "4a")
    again = integrate(tmp_path, 4, "4b")
    assert four == again, "two four-thread runs differ"
    assert four == one, "four threads differ from one"
