"""A report you can look at, rather than a table you have to read.

The trends are numbers about numbers: a ratio, a spread, two correlations and
an outlier count, for nine explanatory variables and five columns.  That is
several hundred figures, and the shape of a disagreement -- flat, trending,
kinked at one end -- is the thing that identifies it.  A shape is not something
a table conveys.

Plotly is loaded from its CDN rather than bundled, and nothing here imports it:
this module writes JSON into a page.  So `mxeq` gains no dependency, the report
is one file that can be mailed to someone, and it needs a network connection
once to draw.  `--offline` inlines the library instead, for a machine that has
none, if `plotly` happens to be installed to take it from.
"""

from __future__ import annotations

import html
import json
from dataclasses import asdict

import numpy as np

from . import match, refl, stats, trends

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.27.0.min.js"


def _figures(a, ia, b, ib, values, n_bins):
    """One figure per value per explanatory variable, as plotly data."""
    variables = trends.explanatory_variables(a, ia, b, ib)
    out = []
    for value in values:
        if value not in a.columns or value not in b.columns:
            continue
        va = a.columns[value].ravel()[ia]
        vb = b.columns[value].ravel()[ib]
        finite = np.isfinite(va) & np.isfinite(vb) & (np.abs(vb) > 0)
        overall, _ = stats.correlation(va[finite], vb[finite])
        blocks = []
        for against in variables:
            rows = trends.trend(va, vb, against, n_bins)
            # A variable that does not vary gives one bin, and a panel with one
            # point in it says nothing and takes up half a screen saying it.
            # The text report still lists them, where a short table costs
            # nothing.
            if len(rows) < 2:
                continue
            centres = [0.5 * (r.low + r.high) for r in rows]
            # Resolution is drawn against 1/d^2, which is what makes a
            # resolution axis even: shells of equal 1/d^2 hold equal volumes
            # of reciprocal space. The ticks still read in Angstroms.
            axis = "linear"
            ticks = None
            if against.scale == "inverse_square":
                x = [1.0 / (c * c) if c > 0 else 0.0 for c in centres]
                ticks = {
                    "vals": x,
                    "text": [f"{c:.2f}" for c in centres],
                }
            elif against.scale == "log":
                x = centres
                axis = "log"
            else:
                x = centres
            blocks.append(
                {
                    "against": against.name,
                    "unit": against.unit if against.scale != "inverse_square" else "A",
                    "axis": axis,
                    "ticks": ticks,
                    "x": x,
                    "ratio": [r.median_ratio for r in rows],
                    "spread": [r.spread for r in rows],
                    "rho": [r.rank_correlation for r in rows],
                    "corr": [r.correlation for r in rows],
                    "outliers": [r.outliers for r in rows],
                    "count": [r.count for r in rows],
                }
            )
        # A sample of the pairs themselves, because a scatter shows what a
        # binned median cannot: whether the disagreement is a few points or
        # all of them.
        take = np.where(finite)[0]
        if take.size > 4000:
            take = np.random.default_rng(0).choice(take, 4000, replace=False)
        out.append(
            {
                "value": value,
                "correlation": None if not np.isfinite(overall) else float(overall),
                "median_ratio": float(np.median(va[finite] / vb[finite]))
                if finite.any()
                else float("nan"),
                "n": int(finite.sum()),
                "blocks": blocks,
                "scatter": {
                    "ours": va[take].tolist(),
                    "theirs": vb[take].tolist(),
                },
            }
        )
    return out


PAGE = """<!doctype html>
<meta charset="utf-8">
<title>__TITLE__</title>
<script src="__PLOTLY__"></script>
<style>
 body {{ font: 14px/1.5 system-ui, sans-serif; margin: 2rem auto; max-width: 1100px;
        color: #222; }}
 h1 {{ font-size: 1.4rem; }}
 h2 {{ font-size: 1.1rem; margin-top: 2.5rem; border-top: 1px solid #ddd;
       padding-top: 1rem; }}
 .summary {{ background: #f6f6f6; padding: 0.6rem 1rem; border-radius: 4px; }}
 .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 0.5rem; }}
 .note {{ color: #666; font-size: 0.9rem; }}
 code {{ background: #f0f0f0; padding: 0 0.2rem; }}
</style>
<h1>__TITLE__</h1>
<p class="note">__SUBTITLE__</p>
<div id="body"></div>
<script>
const REPORT = __DATA__;

// Resolution is drawn against 1/d^2 with the ticks still reading in Angstroms,
// and counts that span decades on a log axis. A linear axis in d compresses
// every high-resolution shell into the left of the plot, which is where most
// of the reflections are.
function axis(block) {
  const unit = block.unit ? " (" + block.unit + ")" : "";
  const out = {title: block.against + unit};
  if (block.axis === "log") out.type = "log";
  if (block.ticks) {
    out.tickmode = "array";
    out.tickvals = block.ticks.vals;
    out.ticktext = block.ticks.text;
  }
  return out;
}

function panel(where, spec) {
  const d = document.createElement("div");
  where.appendChild(d);
  Plotly.newPlot(d, spec.traces, spec.layout,
                 {displayModeBar: false, responsive: true});
}

const body = document.getElementById("body");
for (const value of REPORT.values) {
  const h = document.createElement("h2");
  h.textContent = value.value;
  body.appendChild(h);

  const s = document.createElement("p");
  s.className = "summary";
  s.textContent =
    value.n + " reflections, correlation " +
    (value.correlation === null ? "n/a" : value.correlation.toFixed(4)) +
    ", median ratio " + value.median_ratio.toFixed(4);
  body.appendChild(s);

  // Ours against theirs, on log axes so four decades fit.
  panel(body, {
    traces: [
      {x: value.scatter.theirs, y: value.scatter.ours, mode: "markers",
       type: "scattergl", name: "reflections",
       marker: {size: 3, opacity: 0.35, color: "#1f77b4"}},
      {x: [0.01, 1e6], y: [0.01, 1e6], mode: "lines", name: "equality",
       line: {dash: "dot", width: 1, color: "#999"}}
    ],
    layout: {
      height: 380, title: "ours against theirs",
      xaxis: {title: "theirs", type: "log"},
      yaxis: {title: "ours", type: "log"},
      margin: {t: 40, l: 60, r: 20, b: 50}
    }
  });

  const grid = document.createElement("div");
  grid.className = "grid";
  body.appendChild(grid);

  for (const block of value.blocks) {
    // The ratio with its robust spread as the band: what the bulk does, and
    // how widely, as against what the worst of it does.
    panel(grid, {
      traces: [
        {x: block.x, y: block.ratio, mode: "lines+markers", name: "median ratio",
         error_y: {type: "data", array: block.spread, visible: true,
                   thickness: 1, width: 3, color: "#aaa"},
         line: {color: "#1f77b4"}},
        {x: [block.x[0], block.x[block.x.length - 1]], y: [1, 1],
         mode: "lines", name: "agreement",
         line: {dash: "dot", width: 1, color: "#999"}}
      ],
      layout: {
        height: 300, title: "ratio against " + block.against,
        xaxis: axis(block),
        yaxis: {title: "ours / theirs"},
        showlegend: false, margin: {t: 40, l: 60, r: 20, b: 50}
      }
    });

    // Rank correlation and outlier count together: a bin where rho is high and
    // the outlier count is too is a few disasters, not general scatter.
    panel(grid, {
      traces: [
        {x: block.x, y: block.rho, mode: "lines+markers", name: "rank corr",
         line: {color: "#2ca02c"}},
        {x: block.x, y: block.corr, mode: "lines+markers", name: "pearson",
         line: {color: "#d62728", dash: "dash"}},
        {x: block.x, y: block.outliers, mode: "lines+markers", name: "outliers",
         yaxis: "y2", line: {color: "#ff7f0e", width: 1}}
      ],
      layout: {
        height: 300, title: "agreement against " + block.against,
        xaxis: axis(block),
        yaxis: {title: "correlation", range: [-0.05, 1.05]},
        yaxis2: {title: "outliers", overlaying: "y", side: "right",
                 showgrid: false},
        legend: {orientation: "h", y: -0.25},
        margin: {t: 40, l: 60, r: 60, b: 70}
      }
    });
  }
}
</script>
"""


def write(
    ours: refl.ReflectionTable,
    theirs: refl.ReflectionTable,
    path: str,
    values: list[str],
    n_bins: int = 12,
    radius: float = 5.0,
    title: str = "integration comparison",
    plotly_src: str = PLOTLY_CDN,
) -> str:
    ia, ib, unpartnered = trends.paired(ours, theirs, radius)

    data = {"values": _figures(ours, ia, theirs, ib, values, n_bins)}
    subtitle = (
        f"{ours.nrows} rows against {theirs.nrows}, matched {len(ia)}"
        + (f", {unpartnered} unpartnered" if unpartnered else "")
        + ". Bins hold equal populations. The band on the ratio is the robust "
        "spread, not the standard error."
    )
    page = (
        PAGE.replace("__TITLE__", html.escape(title))
        .replace("__SUBTITLE__", html.escape(subtitle))
        .replace("__PLOTLY__", html.escape(plotly_src, quote=True))
        .replace("__DATA__", json.dumps(data))
    )
    with open(path, "w") as f:
        f.write(page)
    return f"{subtitle}\nwrote {path}"
