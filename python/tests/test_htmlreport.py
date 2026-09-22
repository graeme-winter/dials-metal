"""The HTML comparison report."""

import json
import re

import numpy as np
import pytest

from mxeq import htmlreport
from test_trends import table


def payload(path):
    """The data the page carries, pulled back out of it."""
    text = path.read_text()
    found = re.search(r"const REPORT = (.*?);\n", text, re.S)
    assert found, "the page has no data in it"
    return json.loads(found.group(1)), text


def test_the_page_carries_the_numbers_it_draws(tmp_path):
    a = table(600, seed=31, ratio=0.9)
    b = table(600, seed=31)
    out = tmp_path / "r.html"
    htmlreport.write(a, b, str(out), ["intensity.sum.value"], n_bins=5)
    data, text = payload(out)
    assert len(data["values"]) == 1
    value = data["values"][0]
    assert value["value"] == "intensity.sum.value"
    assert value["n"] == 600
    assert value["median_ratio"] == pytest.approx(0.9)
    # Every block has one point per bin and a spread for each.
    for block in value["blocks"]:
        n = len(block["x"])
        assert n > 1
        for field in ("ratio", "spread", "rho", "corr", "outliers", "count"):
            assert len(block[field]) == n, field
    assert "plot.ly" in text or "plotly" in text


def test_a_trend_reaches_the_page_as_a_trend(tmp_path):
    # The whole point of drawing it: a disagreement that depends on resolution
    # must arrive as a monotonic series, not as one averaged number.
    a = table(600, seed=32, ratio_with_d=0.6)
    b = table(600, seed=32)
    out = tmp_path / "r.html"
    htmlreport.write(a, b, str(out), ["intensity.sum.value"], n_bins=5)
    data, _ = payload(out)
    blocks = {b["against"]: b for b in data["values"][0]["blocks"]}
    assert "resolution" in blocks
    ratio = blocks["resolution"]["ratio"]
    assert ratio[0] > ratio[-1] + 0.2, ratio


def test_nothing_matching_is_refused_with_the_reason(tmp_path):
    a = table(100, seed=33)
    b = table(100, seed=33)
    b.columns["miller_index"] = -b.columns["miller_index"] - 7
    with pytest.raises(ValueError, match="reindex"):
        htmlreport.write(a, b, str(tmp_path / "r.html"), ["intensity.sum.value"])


def test_the_scatter_is_sampled_rather_than_complete(tmp_path):
    # A page holding a million points does not open. Four thousand is enough to
    # see whether a disagreement is a few reflections or all of them.
    a = table(9000, seed=34)
    b = table(9000, seed=34)
    out = tmp_path / "r.html"
    htmlreport.write(a, b, str(out), ["intensity.sum.value"], n_bins=4)
    data, _ = payload(out)
    scatter = data["values"][0]["scatter"]
    assert len(scatter["ours"]) == 4000
    assert len(scatter["theirs"]) == 4000


def test_a_column_missing_from_one_side_is_left_out_rather_than_crashing(tmp_path):
    a = table(200, seed=35)
    b = table(200, seed=35)
    out = tmp_path / "r.html"
    htmlreport.write(
        a, b, str(out), ["intensity.sum.value", "intensity.prf.value"], n_bins=4
    )
    data, _ = payload(out)
    assert [v["value"] for v in data["values"]] == ["intensity.sum.value"]


def test_resolution_is_drawn_against_inverse_d_squared(tmp_path):
    # Shells of equal 1/d^2 hold equal volumes of reciprocal space, so that is
    # the axis on which a resolution trend looks like what it is. A linear axis
    # in d compresses every high-resolution shell into the left of the plot,
    # which is where most of the reflections are.
    #
    # The ticks still read in Angstroms, because that is what anyone means by
    # resolution.
    a = table(600, seed=41)
    b = table(600, seed=41)
    out = tmp_path / "r.html"
    htmlreport.write(a, b, str(out), ["intensity.sum.value"], n_bins=5)
    data, _ = payload(out)
    blocks = {block["against"]: block for block in data["values"][0]["blocks"]}
    resolution = blocks["resolution"]
    assert resolution["ticks"] is not None
    # x is 1/d^2 and the labels are d, so each label squared times its x is one.
    for x, label in zip(resolution["x"], resolution["ticks"]["text"]):
        d = float(label)
        assert x == pytest.approx(1.0 / (d * d), rel=1e-6)
    # Ascending in 1/d^2, which is descending in d: low resolution on the left.
    assert resolution["x"] == sorted(resolution["x"])
    assert float(resolution["ticks"]["text"][0]) > float(
        resolution["ticks"]["text"][-1]
    )


def test_counts_that_span_decades_are_drawn_on_a_log_axis(tmp_path):
    a = table(600, seed=42)
    b = table(600, seed=42)
    # A pixel count spanning three decades, as a real one does between a
    # reflection near the rotation axis and one far from it.
    counts = np.geomspace(20, 20000, 600)
    a.columns["num_pixels.foreground"] = counts.astype(np.int64)
    b.columns["num_pixels.foreground"] = counts.astype(np.int64)
    out = tmp_path / "r.html"
    htmlreport.write(a, b, str(out), ["intensity.sum.value"], n_bins=5)
    data, _ = payload(out)
    blocks = {block["against"]: block for block in data["values"][0]["blocks"]}
    assert blocks["foreground pixels (ours)"]["axis"] == "log"
    # And a variable that does not span decades is left alone.
    assert blocks["I/sigma (reference)"]["axis"] == "linear"
