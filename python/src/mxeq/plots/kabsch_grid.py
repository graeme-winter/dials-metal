"""Plot the spot density in Kabsch space, from mxi_grid.

Three sections through each grid: the detector plane (eps1, eps2) summed over
the rotation direction, and the two planes that contain it. A single strong
spot is only a handful of pixels across, so its grid is coarse and noisy; the
reference beside it is the average of the spots nearest it on the detector,
which is what says whether a feature belongs to the spot or to the sampling.

Contours are drawn at fractions of each panel's own maximum, so shapes can be
compared between panels whose totals differ.
"""

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read(path):
    """Read the grids mxi_grid wrote: a header, then labelled cubes."""
    header = {}
    grids = {}
    with open(path) as handle:
        lines = [line.rstrip("\n") for line in handle]
    words = lines[0].lstrip("# ").split()
    for i in range(0, len(words) - 1, 2):
        header[words[i]] = float(words[i + 1])
    side = int(header["side"])
    n = side**3
    at = 1
    while at < len(lines):
        parts = lines[at].split()
        label, spots = parts[0], int(parts[1])
        values = np.array([float(v) for v in lines[at + 1 : at + 1 + n]])
        grids[label] = dict(
            spots=spots,
            # Written with eps1 fastest, then eps2, then eps3.
            cube=values.reshape(side, side, side),
        )
        at += 1 + n
    return header, grids


def draw(ax, image, extent, title, xlabel, ylabel):
    if image.max() <= 0:
        ax.set_axis_off()
        return
    ax.imshow(
        image.T,
        origin="lower",
        extent=extent,
        aspect="auto",
        cmap="magma",
        interpolation="nearest",
    )
    # Contours at fractions of this panel's own maximum, so the shape can be
    # compared with a panel holding a different number of counts.
    levels = [0.1, 0.3, 0.5, 0.7, 0.9]
    ax.contour(
        np.linspace(extent[0], extent[1], image.shape[0]),
        np.linspace(extent[2], extent[3], image.shape[1]),
        image.T / image.max(),
        levels=levels,
        colors="w",
        linewidths=0.6,
        alpha=0.7,
    )
    ax.axhline(0, color="w", lw=0.4, alpha=0.4)
    ax.axvline(0, color="w", lw=0.4, alpha=0.4)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(labelsize=7)


def main(path, out):
    header, grids = read(path)
    sigma_b = header["sigma_b"]
    sigma_m = header["sigma_m"]
    half = header["half_width"]

    # In sigmas, which is the only scale on which the two directions are
    # comparable: eps1 and eps2 are measured against sigma_b and eps3 against
    # sigma_m, and they differ by a factor of four here.
    d = (-half, half, -half, half)

    names = [k for k in grids if k.startswith("spot_")]
    names.sort()
    rows = [("all", "average of all %d spots" % grids["all"]["spots"])]
    for i, name in enumerate(names):
        rows.append((name, name.replace("_", " ")))
        rows.append(("reference_%d" % i, "average of its %d nearest neighbours"
                     % grids["reference_%d" % i]["spots"]))

    fig, axes = plt.subplots(len(rows), 3, figsize=(11, 3.1 * len(rows)),
                             squeeze=False)
    for r, (key, label) in enumerate(rows):
        cube = grids[key]["cube"]
        # cube is indexed [eps3, eps2, eps1].
        draw(axes[r][0], cube.sum(axis=0).T, d,
             "%s\ndetector plane" % label, "eps1 / sigma_b", "eps2 / sigma_b")
        draw(axes[r][1], cube.sum(axis=1).T, d,
             "rotation against eps1", "eps1 / sigma_b", "eps3 / sigma_m")
        draw(axes[r][2], cube.sum(axis=2).T, d,
             "rotation against eps2", "eps2 / sigma_b", "eps3 / sigma_m")

    fig.suptitle(
        "spot density in Kabsch space, sigma_b = %.4f deg, sigma_m = %.4f deg\n"
        "axes in sigmas; a Gaussian would fall to 0.1 of its peak at 2.1 sigma"
        % (sigma_b, sigma_m),
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=110)
    print("wrote", out)

    # The picture as a number: how far out, in sigmas, the average profile
    # falls to a tenth of its peak. A Gaussian does that at 2.146 sigma.
    cube = grids["all"]["cube"]
    side = cube.shape[0]
    axis = np.linspace(-half, half, side)
    print()
    print("average profile, radius at which it falls to 0.1 of its peak")
    print("  a Gaussian would be at 2.146 sigma")
    for name, marginal in (
        ("eps1 (sigma_b)", cube.sum(axis=(0, 1))),
        ("eps2 (sigma_b)", cube.sum(axis=(0, 2))),
        ("eps3 (sigma_m)", cube.sum(axis=(1, 2))),
    ):
        profile = marginal / marginal.max()
        above = np.abs(axis[profile >= 0.1])
        print("  %-16s %.2f sigma" % (name, above.max() if above.size else 0.0))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
