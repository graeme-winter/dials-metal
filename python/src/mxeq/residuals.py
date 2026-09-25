"""Pictures of how well positions were predicted, from xyzres.px.

`mxi_integrate` writes, for every reflection it could find a centre of mass
for, the observed centre less the predicted one and its variance.  That is a
residual in the images' own frame -- pixels on the detector and images along
the scan, parallax included -- and it answers a different question depending on
what it is plotted against:

* against IMAGE NUMBER, whether the scan-varying model has drifted: a median
  that wanders along the scan is the crystal or the beam moving in a way the
  refinement did not follow;
* against RESOLUTION, whether the cell is right: a residual that grows outward
  is a cell or a distance error, since both scale positions with 1/d;
* across the DETECTOR, whether the detector model is right: a pattern that
  follows the panel rather than the reflections is a tilt, a distance, or the
  parallax;
* and as a PULL, residual over sigma, whether the prediction error is larger
  than the counting noise at all.

Strong reflections by default, since for them counting noise is a hundredth of
a pixel and what is left is the prediction.  The page loads plotly from its
CDN, as `mxeq html` does, and nothing here imports it.
"""

from __future__ import annotations

import html
import json

import numpy as np

from . import refl
from .htmlreport import PLOTLY_CDN


def _binned(x: np.ndarray, y: np.ndarray, n_bins: int):
    """Median of y and its robust spread in equal-population bins of x."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 2 * n_bins:
        n_bins = max(1, x.size // 20)
    edges = np.unique(np.percentile(x, np.linspace(0, 100, n_bins + 1)))
    centres, medians, spreads, counts = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        last = hi == edges[-1]
        inside = (x >= lo) & ((x <= hi) if last else (x < hi))
        if inside.sum() < 5:
            continue
        values = y[inside]
        middle = float(np.median(values))
        centres.append(float(np.median(x[inside])))
        medians.append(middle)
        spreads.append(float(1.4826 * np.median(np.abs(values - middle))))
        counts.append(int(inside.sum()))
    return centres, medians, spreads, counts


def _map(fast, slow, value, image_size, cells):
    """Median of value over a grid of detector cells, as a heatmap."""
    fx = np.clip((fast / image_size[0] * cells).astype(int), 0, cells - 1)
    sy = np.clip((slow / image_size[1] * cells).astype(int), 0, cells - 1)
    grid = np.full((cells, cells), np.nan)
    count = np.zeros((cells, cells), dtype=int)
    for i in range(cells):
        for j in range(cells):
            inside = (sy == i) & (fx == j)
            count[i, j] = int(inside.sum())
            # A cell with too few reflections shows as empty rather than as a
            # confident number drawn from three of them.
            if count[i, j] >= 10:
                grid[i, j] = float(np.median(value[inside]))
    return grid, count


def gather(table: refl.ReflectionTable, least_signal: float):
    """The rows worth plotting, and what to plot them against."""
    if "xyzres.px.value" not in table.columns:
        raise ValueError(
            "no xyzres.px.value: this table was not written by an mxi_integrate "
            "new enough to record position residuals"
        )
    r = table.columns["xyzres.px.value"]
    v = table.columns["xyzres.px.variance"]
    ok = np.all(np.isfinite(r), axis=1) & np.all(np.isfinite(v), axis=1)
    ok &= np.all(v > 0, axis=1)

    isig = None
    if "intensity.sum.value" in table.columns:
        I = table.columns["intensity.sum.value"].ravel()
        V = table.columns["intensity.sum.variance"].ravel()
        with np.errstate(invalid="ignore", divide="ignore"):
            isig = I / np.sqrt(np.maximum(V, 1e-12))
        ok &= isig >= least_signal

    rows = np.where(ok)[0]
    return {
        "rows": rows,
        "residual": r[rows],
        "variance": v[rows],
        "fast": table.columns["xyzcal.px"][rows, 0],
        "slow": table.columns["xyzcal.px"][rows, 1],
        "image": table.columns["xyzcal.px"][rows, 2],
        "d": table.columns["d"].ravel()[rows] if "d" in table.columns else None,
        "total": table.nrows,
        "with_centre": int(
            (np.all(np.isfinite(r), axis=1)).sum()
        ),
    }


AXES = ("fast", "slow", "image")
UNITS = ("px", "px", "images")


def build(table: refl.ReflectionTable, least_signal: float = 10.0, n_bins: int = 30,
          cells: int = 12, image_size=None) -> dict:
    g = gather(table, least_signal)
    n = len(g["rows"])
    if n == 0:
        raise ValueError(
            f"no reflections with a residual and I/sigma of at least "
            f"{least_signal:g}; lower --least-signal"
        )
    r, v = g["residual"], g["variance"]
    out = {"n": n, "total": g["total"], "with_centre": g["with_centre"],
           "least_signal": least_signal, "axes": list(AXES), "units": list(UNITS)}

    # The summary a reader wants first: the systematic offset and how much
    # scatter there is beyond counting noise.
    median = np.median(r, axis=0)
    spread = 1.4826 * np.median(np.abs(r - median), axis=0)
    counting = np.sqrt(np.median(v, axis=0))
    out["summary"] = [
        {"axis": a, "unit": u, "median": float(median[k]), "spread": float(spread[k]),
         "counting": float(counting[k]),
         "prediction": float(np.sqrt(max(spread[k] ** 2 - counting[k] ** 2, 0.0)))}
        for k, (a, u) in enumerate(zip(AXES, UNITS))
    ]

    out["against_image"] = [
        dict(zip(("x", "median", "spread", "count"), _binned(g["image"], r[:, k], n_bins)))
        for k in range(3)
    ]
    if g["d"] is not None:
        inverse = 1.0 / np.maximum(g["d"], 1e-9) ** 2
        blocks = []
        for k in range(3):
            x, med, spr, cnt = _binned(inverse, r[:, k], n_bins)
            blocks.append({"x": x, "median": med, "spread": spr, "count": cnt,
                           "ticks": [f"{1.0 / np.sqrt(t):.2f}" if t > 0 else "" for t in x]})
        out["against_resolution"] = blocks

    if image_size is None:
        image_size = (float(np.max(g["fast"])) + 1.0, float(np.max(g["slow"])) + 1.0)
    maps = []
    for k in range(2):
        grid, count = _map(g["fast"], g["slow"], r[:, k], image_size, cells)
        maps.append({"z": [[None if not np.isfinite(c) else c for c in row] for row in grid],
                     "count": count.tolist()})
    out["maps"] = maps
    out["image_size"] = list(image_size)

    pulls = []
    for k in range(3):
        p = r[:, k] / np.sqrt(v[:, k])
        p = p[np.isfinite(p)]
        clipped = p[np.abs(p) < 20]
        pulls.append({"values": clipped.tolist() if clipped.size <= 20000
                      else np.random.default_rng(0).choice(clipped, 20000, replace=False).tolist(),
                      "rms": float(np.sqrt(np.mean(p ** 2))) if p.size else float("nan")})
    out["pulls"] = pulls
    return out


PAGE = """<!doctype html>
<meta charset="utf-8">
<title>__TITLE__</title>
<script src="__PLOTLY__"></script>
<style>
 body { font: 14px/1.5 system-ui, sans-serif; margin: 2rem auto; max-width: 1100px; color: #222; }
 h1 { font-size: 1.4rem; } h2 { font-size: 1.1rem; margin-top: 2.2rem; border-top: 1px solid #ddd; padding-top: 1rem; }
 table { border-collapse: collapse; margin: 0.5rem 0 1rem; }
 td, th { padding: 0.2rem 0.8rem; text-align: right; border-bottom: 1px solid #eee; }
 th { font-weight: 600; } .note { color: #666; font-size: 0.9rem; }
 .grid { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 0.4rem; }
 .grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 0.4rem; }
</style>
<h1>__TITLE__</h1>
<p class="note">__SUBTITLE__</p>
<div id="body"></div>
<script>
const R = __DATA__;
const body = document.getElementById("body");
const cfg = {displayModeBar: false, responsive: true};
function el(tag, cls, text) { const e = document.createElement(tag); if (cls) e.className = cls; if (text) e.textContent = text; return e; }
function plot(where, traces, layout) { const d = el("div"); where.appendChild(d); Plotly.newPlot(d, traces, layout, cfg); }

body.appendChild(el("h2", null, "Summary"));
const t = el("table");
t.innerHTML = "<tr><th></th><th>median</th><th>spread</th><th>counting</th><th>prediction</th></tr>" +
  R.summary.map(s => `<tr><td style="text-align:left">${s.axis} (${s.unit})</td>` +
    `<td>${s.median.toFixed(3)}</td><td>${s.spread.toFixed(3)}</td>` +
    `<td>${s.counting.toFixed(3)}</td><td>${s.prediction.toFixed(3)}</td></tr>`).join("");
body.appendChild(t);
body.appendChild(el("p", "note",
  "Median is the systematic offset. Spread is the robust scatter; counting is the median sigma " +
  "from the counts; prediction is what the spread has beyond counting, sqrt(spread^2 - counting^2)."));

function against(title, blocks, xtitle, useTicks) {
  body.appendChild(el("h2", null, title));
  const g = el("div", "grid"); body.appendChild(g);
  blocks.forEach((b, k) => {
    const upper = b.median.map((m, i) => m + b.spread[i]);
    const lower = b.median.map((m, i) => m - b.spread[i]);
    const xa = {title: xtitle};
    if (useTicks) { xa.tickmode = "array"; xa.tickvals = b.x.filter((_, i) => i % 5 === 0); xa.ticktext = b.ticks.filter((_, i) => i % 5 === 0); }
    plot(g, [
      {x: b.x.concat(b.x.slice().reverse()), y: upper.concat(lower.reverse()), fill: "toself",
       fillcolor: "rgba(31,119,180,0.15)", line: {width: 0}, hoverinfo: "skip", showlegend: false},
      {x: b.x, y: b.median, mode: "lines+markers", line: {color: "#1f77b4"}, marker: {size: 4}, showlegend: false},
      {x: [b.x[0], b.x[b.x.length - 1]], y: [0, 0], mode: "lines", line: {dash: "dot", color: "#999", width: 1}, showlegend: false}
    ], {height: 280, title: R.axes[k] + " residual (" + R.units[k] + ")", xaxis: xa,
        margin: {t: 36, l: 50, r: 10, b: 45}});
  });
}
against("Against image number", R.against_image, "image", false);
if (R.against_resolution) against("Against resolution", R.against_resolution, "d (A)", true);

body.appendChild(el("h2", null, "Across the detector"));
const m = el("div", "grid2"); body.appendChild(m);
R.maps.forEach((map, k) => {
  const flat = map.z.flat().filter(v => v !== null);
  const lim = Math.max(...flat.map(Math.abs), 1e-6);
  plot(m, [{z: map.z, type: "heatmap", colorscale: "RdBu", zmin: -lim, zmax: lim, reversescale: true,
            colorbar: {title: "px"}}],
       {height: 380, title: "median " + R.axes[k] + " residual",
        xaxis: {title: "fast cell", constrain: "domain"}, yaxis: {title: "slow cell", autorange: "reversed", scaleanchor: "x"},
        margin: {t: 36, l: 50, r: 10, b: 45}});
});

body.appendChild(el("h2", null, "Pull: residual over sigma"));
body.appendChild(el("p", "note",
  "For a perfect prediction and an honest sigma this is a unit Gaussian. Wider means prediction error " +
  "larger than the counting noise; the rms says by how much."));
const p = el("div", "grid"); body.appendChild(p);
const gx = [], gy = [];
for (let x = -6; x <= 6; x += 0.1) { gx.push(x); gy.push(Math.exp(-x * x / 2) / Math.sqrt(2 * Math.PI)); }
R.pulls.forEach((pull, k) => {
  plot(p, [
    {x: pull.values, type: "histogram", histnorm: "probability density", xbins: {start: -10, end: 10, size: 0.25},
     marker: {color: "rgba(31,119,180,0.6)"}, name: "pull"},
    {x: gx, y: gy, mode: "lines", line: {color: "#d62728", width: 1.5}, name: "unit Gaussian"}
  ], {height: 280, title: R.axes[k] + ": rms " + pull.rms.toFixed(2), xaxis: {title: "pull", range: [-10, 10]},
      showlegend: false, margin: {t: 36, l: 50, r: 10, b: 45}});
});
</script>
"""


def write(table: refl.ReflectionTable, path: str, least_signal: float = 10.0,
          n_bins: int = 30, cells: int = 12, title: str = "position residuals",
          plotly_src: str = PLOTLY_CDN) -> str:
    data = build(table, least_signal, n_bins, cells)
    subtitle = (
        f"{data['n']} reflections with I/sigma of at least {least_signal:g}, of "
        f"{data['with_centre']} with a centre of mass and {data['total']} in all. "
        "Observed centre less predicted, in the images' own frame. Bands are the "
        "robust spread."
    )
    page = (PAGE.replace("__TITLE__", html.escape(title))
            .replace("__SUBTITLE__", html.escape(subtitle))
            .replace("__PLOTLY__", html.escape(plotly_src, quote=True))
            .replace("__DATA__", json.dumps(data)))
    with open(path, "w") as f:
        f.write(page)
    lines = [subtitle, ""]
    lines.append("  %-8s %9s %9s %9s %11s" % ("axis", "median", "spread", "counting", "prediction"))
    for s in data["summary"]:
        lines.append("  %-8s %9.3f %9.3f %9.3f %11.3f" % (
            s["axis"], s["median"], s["spread"], s["counting"], s["prediction"]))
    lines.append("")
    lines.append(f"wrote {path}")
    return "\n".join(lines)
