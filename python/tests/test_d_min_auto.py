"""mxi_scale --d-min-auto is the --d-min it reports, byte for byte.

--d-min-auto scales everything, estimates the resolution from CC half, and
scales again to the limit, rounded to a hundredth of an angstrom so that what
it chose can be said and repeated. It says "as --d-min X would"; this holds it
to that: the table from --d-min-auto and the one from --d-min X are the same
bytes. At a CC half limit of 0.99 so that the 300 image sweep crosses it; at the
default 0.3 its CC half never falls that far and nothing is cut.

The equivalence needs data to scale, MXI_SCALE_EXPT and MXI_SCALE_REFL (a
symmetrized.expt and .refl); the refusal of both options needs only
MXI_SCALE, the program.
"""

import os
import re
import subprocess

import pytest

BINARY = os.environ.get("MXI_SCALE")
EXPT = os.environ.get("MXI_SCALE_EXPT")
REFL = os.environ.get("MXI_SCALE_REFL")

needs_binary = pytest.mark.skipif(not BINARY, reason="set MXI_SCALE")
needs_data = pytest.mark.skipif(
    not (BINARY and EXPT and REFL),
    reason="set MXI_SCALE, MXI_SCALE_EXPT and MXI_SCALE_REFL",
)


def scale(tmp_path, name, *extra):
    run = subprocess.run(
        [
            BINARY,
            EXPT,
            REFL,
            "-o",
            name + ".refl",
            "--output-expt",
            name + ".expt",
            "--cc-half-limit",
            "0.99",
            *extra,
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    return run.stdout


@needs_binary
def test_both_limits_are_refused(tmp_path):
    run = subprocess.run(
        [BINARY, "any.expt", "any.refl", "--d-min", "2", "--d-min-auto"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 2, run.stderr
    assert "two answers to one question" in run.stderr


@needs_data
def test_the_automatic_limit_is_the_limit_it_reports(tmp_path):
    out = scale(tmp_path, "auto", "--d-min-auto")
    found = re.search(r"as --d-min ([0-9.]+) would", out)
    assert found, out
    scale(tmp_path, "given", "--d-min", found.group(1))
    auto = (tmp_path / "auto.refl").read_bytes()
    given = (tmp_path / "given.refl").read_bytes()
    assert auto == given, "--d-min-auto and the --d-min it reported differ"
    assert "Suggested" not in out
