"""Profile fitting in single precision -- mxi_integrate --gpu -- agrees with the
double-precision fit.

--gpu-emulate runs on the CPU exactly the single-precision steps a device runs
(src/fit_device.hh), without the CPU's exact recomputation of a subdivision near
a cell boundary, which float could not make exact anyway. On the 300 image sweep
it differs from the double fit by a median of 9e-6 sigma and at most 0.036, the
same reflections fitted. This holds it to thirty times looser than that, so
that it catches the two drifting apart and not rounding.

Needs the matched data the post-refinement test uses: MXI_INTEGRATE,
MXI_POSTREFINE_EXPT and MXI_POSTREFINE_REFL.
"""

import os
import subprocess

import numpy as np
import pytest

from mxeq import refl

BINARY = os.environ.get("MXI_INTEGRATE")
EXPT = os.environ.get("MXI_POSTREFINE_EXPT")
REFL = os.environ.get("MXI_POSTREFINE_REFL")
PRF = 1 << 9

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
            "--threads",
            "2",
            "-o",
            name + ".refl",
            "--output-expt",
            name + ".expt",
            *extra,
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    return refl.load(str(tmp_path / (name + ".refl")))


def column(table, name):
    try:
        return np.asarray(table[name])
    except Exception:
        return np.asarray(table.columns[name])


@needs_data
def test_single_precision_fitting_agrees_with_double(tmp_path):
    double = integrate(tmp_path, "double")
    single = integrate(tmp_path, "single", "--gpu-emulate")
    fa = (column(double, "flags").astype(np.int64) & PRF) != 0
    fb = (column(single, "flags").astype(np.int64) & PRF) != 0
    assert fa.sum() > 100, "a silent test: nothing profile fitted"
    assert (
        fa == fb
    ).all(), f"{int((fa != fb).sum())} reflections fitted by one and not the other"
    ia = column(double, "intensity.prf.value").astype(float)[fa]
    ib = column(single, "intensity.prf.value").astype(float)[fa]
    sigma = np.sqrt(column(double, "intensity.prf.variance").astype(float)[fa])
    d = np.abs(ia - ib) / sigma
    assert np.median(d) < 1e-3, f"median {np.median(d)} sigma"
    assert d.max() < 0.1, f"largest {d.max()} sigma"
