"""One pass over the images gives the two passes' answer, byte for byte.

mxi_integrate fits each reflection as soon as the scan blocks its profile is
interpolated from are complete, holding its shoebox until then, and reads every
frame once; --two-pass reads them twice, to learn and then to fit. The profiles
are learned in the same order either way and a held shoebox is the one the
second pass would rebuild, so the tables must be the same bytes. Sixty images
in six scan blocks, so that blocks are finalised and shoeboxes released all
through the stream, and with the shoeboxes saved, which is where an ordering
mistake once overwrote them.

Needs the matched data the post-refinement test uses: MXI_INTEGRATE,
MXI_POSTREFINE_EXPT and MXI_POSTREFINE_REFL.
"""

import os
import subprocess

import pytest

BINARY = os.environ.get("MXI_INTEGRATE")
EXPT = os.environ.get("MXI_POSTREFINE_EXPT")
REFL = os.environ.get("MXI_POSTREFINE_REFL")

needs_data = pytest.mark.skipif(
    not (BINARY and EXPT and REFL),
    reason="set MXI_INTEGRATE, MXI_POSTREFINE_EXPT and MXI_POSTREFINE_REFL",
)


def integrate(tmp_path, name, *extra):
    run = subprocess.run(
        [
            BINARY,
            EXPT,
            REFL,
            "--last-image",
            "60",
            "--scan-blocks",
            "6",
            "--save-shoeboxes",
            "--threads",
            "2",
            "-o",
            name + ".refl",
            "--output-expt",
            name + ".expt",
            "--timing",
            *extra,
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    return run.stdout


@needs_data
def test_one_pass_is_the_two_passes_byte_for_byte(tmp_path):
    one = integrate(tmp_path, "one")
    two = integrate(tmp_path, "two", "--two-pass")
    a = (tmp_path / "one.refl").read_bytes()
    b = (tmp_path / "two.refl").read_bytes()
    assert a == b, "one pass and two gave different tables"
    # And one pass read each frame once, where two read it twice.
    assert "each read 1.00 times" in one, one
    assert "each read 2.00 times" in two, two
    assert "one pass: at most" in one
