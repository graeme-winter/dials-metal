"""The position residual report."""

import json
import re

import numpy as np
import pytest

from mxeq import refl, residuals


def table(
    n=4000, seed=0, offset=(0.0, 0.0, 0.0), error=(0.2, 0.2, 0.1), drift=0.0, sigma=0.05
):
    """An integrated table with planted residuals of known size.

    `offset` is a systematic, `error` the prediction scatter, `drift` a
    residual in fast that grows along the scan, and `sigma` the counting noise.
    """
    rng = np.random.default_rng(seed)
    t = refl.ReflectionTable(nrows=n)
    image = rng.uniform(0, 1800, n)
    fast = rng.uniform(0, 2000, n)
    slow = rng.uniform(0, 2000, n)
    t.columns["xyzcal.px"] = np.stack([fast, slow, image], axis=1)
    counting = rng.normal(0.0, sigma, (n, 3))
    prediction = rng.normal(0.0, 1.0, (n, 3)) * np.asarray(error)
    r = np.asarray(offset) + prediction + counting
    r[:, 0] += drift * (image - 900.0) / 900.0
    t.columns["xyzres.px.value"] = r
    t.columns["xyzres.px.variance"] = np.full((n, 3), sigma * sigma)
    t.columns["intensity.sum.value"] = np.full(n, 1000.0)
    t.columns["intensity.sum.variance"] = np.full(n, 1000.0)
    t.columns["d"] = rng.uniform(1.5, 10.0, n)
    return t


def payload(path):
    text = path.read_text()
    return json.loads(re.search(r"const R = (.*?);\n", text, re.S).group(1))


def test_the_summary_separates_the_offset_the_prediction_and_the_counting(tmp_path):
    # A planted systematic, a planted prediction scatter and planted counting
    # noise, each of known size: the report must pull them apart, since that
    # separation is the whole reason for having a variance in the table.
    t = table(offset=(0.3, -0.1, 0.05), error=(0.25, 0.15, 0.08), sigma=0.05)
    out = tmp_path / "r.html"
    residuals.write(t, str(out), least_signal=0.0)
    summary = {s["axis"]: s for s in payload(out)["summary"]}
    for axis, off, err in (
        ("fast", 0.3, 0.25),
        ("slow", -0.1, 0.15),
        ("image", 0.05, 0.08),
    ):
        s = summary[axis]
        assert s["median"] == pytest.approx(off, abs=0.02), axis
        assert s["counting"] == pytest.approx(0.05, rel=0.01), axis
        # What the spread has beyond counting is the prediction error.
        assert s["prediction"] == pytest.approx(err, rel=0.1), axis


def test_a_drift_along_the_scan_shows_as_one(tmp_path):
    # The plot against image number exists to show a scan-varying model that
    # has drifted, so a planted drift has to come out as a monotonic series.
    t = table(drift=0.5, error=(0.05, 0.05, 0.05), sigma=0.01)
    out = tmp_path / "r.html"
    residuals.write(t, str(out), least_signal=0.0, n_bins=10)
    fast = payload(out)["against_image"][0]
    assert fast["median"][0] < -0.3
    assert fast["median"][-1] > 0.3
    assert fast["median"] == sorted(fast["median"])


def test_a_perfect_prediction_gives_a_unit_pull(tmp_path):
    # With no prediction error the only scatter is counting, and residual over
    # sigma is a unit Gaussian. The pull panel is how a reader tells whether
    # prediction error is there at all.
    t = table(error=(0.0, 0.0, 0.0), sigma=0.05)
    out = tmp_path / "r.html"
    residuals.write(t, str(out), least_signal=0.0)
    for pull in payload(out)["pulls"]:
        assert pull["rms"] == pytest.approx(1.0, abs=0.05)


def test_weak_reflections_are_left_out_unless_asked_for(tmp_path):
    t = table(n=1000)
    t.columns["intensity.sum.value"][:600] = 10.0  # I/sigma about 0.3
    out = tmp_path / "r.html"
    residuals.write(t, str(out), least_signal=10.0)
    assert payload(out)["n"] == 400
    residuals.write(t, str(out), least_signal=0.0)
    assert payload(out)["n"] == 1000


def test_a_table_without_residuals_says_what_it_needs():
    t = refl.ReflectionTable(nrows=3)
    t.columns["xyzcal.px"] = np.zeros((3, 3))
    with pytest.raises(ValueError, match="xyzres"):
        residuals.write(t, "/dev/null")


def test_rows_with_no_centre_are_skipped_not_plotted_as_zero(tmp_path):
    # mxi_integrate writes NaN where it found no centre of mass. Those must be
    # left out; drawn as zero they would pull every median toward a perfect
    # prediction.
    t = table(n=2000, offset=(0.5, 0.5, 0.5), error=(0.01, 0.01, 0.01), sigma=0.01)
    t.columns["xyzres.px.value"][:1000] = np.nan
    t.columns["xyzres.px.variance"][:1000] = np.nan
    out = tmp_path / "r.html"
    residuals.write(t, str(out), least_signal=0.0)
    data = payload(out)
    assert data["n"] == 1000
    for s in data["summary"]:
        assert s["median"] == pytest.approx(0.5, abs=0.02)
