"""mxi_integrate --postrefine: integrate, refine against the centres
integration measured, integrate again.

Refinement fitted to the spot finder's centres leaves strong reflections about
0.1 images early, because those centres depend on strength. Post-refinement
removes it: on a 300 image insulin sweep, the median z offset over I/sigma of
ten or more fell from -0.115 images to -0.006.

The end-to-end test needs matched data: MXI_POSTREFINE_EXPT and
MXI_POSTREFINE_REFL, a refined model and its indexed reflections, with
MXI_POSTREFINE_IMAGES if the model's template does not find the images.
"""

import json
import os
import subprocess

import numpy as np
import pytest

from mxeq import refl

BINARY = os.environ.get("MXI_INTEGRATE")
EXPT = os.environ.get("MXI_POSTREFINE_EXPT")
REFL = os.environ.get("MXI_POSTREFINE_REFL")
IMAGES = os.environ.get("MXI_POSTREFINE_IMAGES")

needs_binary = pytest.mark.skipif(not BINARY, reason="set MXI_INTEGRATE")
needs_data = pytest.mark.skipif(
    not (BINARY and EXPT and REFL),
    reason="set MXI_INTEGRATE, MXI_POSTREFINE_EXPT and MXI_POSTREFINE_REFL",
)


@needs_binary
def test_its_options_are_refused_without_it(tmp_path):
    # Refused rather than ignored: a run that silently drops an option has used
    # settings nobody chose.
    for flag, value in (("--output-expt", "x.expt"), ("--postrefine-points", "5")):
        run = subprocess.run(
            [BINARY, "any.expt", flag, value],
            capture_output=True,
            text=True,
            cwd=tmp_path,
        )
        assert run.returncode == 2, run.stderr
        assert "is for --postrefine" in run.stderr


def z_offset(path):
    """Median observed minus predicted z over reflections summed with I/sigma
    of at least ten."""
    t = refl.load(str(path))
    flags = t["flags"].ravel().astype(np.int64)
    i = t["intensity.sum.value"].ravel()
    v = t["intensity.sum.variance"].ravel()
    strong = ((flags & 256) != 0) & (i / np.sqrt(np.maximum(v, 1e-9)) >= 10)
    dz = t["xyzobs.px.value"][:, 2] - t["xyzcal.px"][:, 2]
    return float(np.median(dz[strong]))


def integrate(tmp_path, *extra):
    command = [BINARY, EXPT, REFL, *extra]
    if IMAGES:
        command += ["--images", IMAGES]
    run = subprocess.run(command, capture_output=True, text=True, cwd=tmp_path)
    assert run.returncode == 0, run.stderr
    return run


@needs_data
def test_post_refinement_removes_the_z_offset(tmp_path):
    integrate(tmp_path, "-o", "plain.refl")
    integrate(tmp_path, "-o", "post.refl", "--postrefine", "--output-expt", "post.expt")
    before = z_offset(tmp_path / "plain.refl")
    after = z_offset(tmp_path / "post.refl")
    # A third of what it was, or within 0.02 images of zero where there was
    # little to remove in the first place.
    assert abs(after) <= max(abs(before) / 3, 0.02), (before, after)

    # The models it integrated with are written, as samples DIALS can read.
    written = json.load(open(tmp_path / "post.expt"))
    first, last = written["scan"][0]["image_range"]
    points = written["crystal"][0]["A_at_scan_points"]
    assert len(points) == last - first + 2

    # And the first integration, which was only a means, is not left behind.
    assert not list(tmp_path.glob("*.before-postrefinement.refl"))

    # The cell a post-refinement settles on is the cell it began from, near
    # enough: what it corrects is where spots are in rotation, not their size.
    given = json.load(open(EXPT))
    a0 = np.array(given["crystal"][0]["real_space_a"])
    a1 = np.array(written["crystal"][0]["real_space_a"])
    assert abs(np.linalg.norm(a1) / np.linalg.norm(a0) - 1) < 0.002
