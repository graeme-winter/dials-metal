"""Spot width against obliquity, with the sensor smear that might explain it.

eps2 lies in the scattering plane -- radially on the detector -- and eps1
across it. A photon absorbed at a random depth in the sensor is recorded a
little further out than where its ray entered, so that smear is radial and
belongs to eps2 alone. If it is the cause of the anisotropy then
width2^2 - width1^2 should equal the predicted smear squared, and should grow
with obliquity the way absorption predicts.
"""

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main(path, out, pixel_degrees):
    d = np.loadtxt(path)
    x, y, r, ob, w1, w2, w3, counts, pred = d.T
    keep = (counts > 50) & (w1 > 0) & (w2 > 0)
    x, y, r, ob, w1, w2, pred = (a[keep] for a in (x, y, r, ob, w1, w2, pred))

    edges = np.percentile(ob, np.linspace(0, 100, 13))
    mid, m1, m2, mp, excess = [], [], [], [], []
    for i in range(len(edges) - 1):
        s = (ob >= edges[i]) & (ob <= edges[i + 1])
        if s.sum() < 50:
            continue
        mid.append(np.median(ob[s]))
        m1.append(np.median(w1[s]))
        m2.append(np.median(w2[s]))
        mp.append(np.median(pred[s]))
        excess.append(np.median(w2[s] ** 2 - w1[s] ** 2))
    mid, m1, m2, mp, excess = map(np.array, (mid, m1, m2, mp, excess))

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

    ax[0].plot(mid, m1 / pixel_degrees, "o-", label="width1, tangential")
    ax[0].plot(mid, m2 / pixel_degrees, "s-", label="width2, radial")
    ax[0].axhline(1 / np.sqrt(12), color="k", ls=":",
                  label="one lit pixel, 1/sqrt(12)")
    ax[0].set_xlabel("obliquity (degrees)")
    ax[0].set_ylabel("width (pixels)")
    ax[0].set_title("the spots are about one pixel wide")
    ax[0].legend(fontsize=8)
    ax[0].set_ylim(0, None)

    ax[1].plot(mid, excess, "o-", label=r"measured $w_2^2 - w_1^2$")
    ax[1].plot(mid, mp ** 2, "s--", label="predicted sensor smear$^2$")
    ax[1].set_xlabel("obliquity (degrees)")
    ax[1].set_ylabel("degrees squared")
    ax[1].set_title("the excess does not grow as absorption predicts")
    ax[1].legend(fontsize=8)
    ax[1].set_ylim(0, None)

    # Where on the detector, to see whether the anisotropy is smooth.
    ratio = w2 / w1
    grid = 24
    xi = np.linspace(x.min(), x.max(), grid)
    yi = np.linspace(y.min(), y.max(), grid)
    image = np.full((grid, grid), np.nan)
    for i in range(grid - 1):
        for j in range(grid - 1):
            s = (x >= xi[i]) & (x < xi[i + 1]) & (y >= yi[j]) & (y < yi[j + 1])
            if s.sum() >= 20:
                image[i, j] = np.median(ratio[s])
    im = ax[2].imshow(image.T, origin="lower", cmap="coolwarm",
                      extent=(x.min(), x.max(), y.min(), y.max()),
                      vmin=1.0, vmax=1.5)
    fig.colorbar(im, ax=ax[2], label="width2 / width1")
    ax[2].set_xlabel("fast (mm)")
    ax[2].set_ylabel("slow (mm)")
    ax[2].set_title("and it is not radial on the face")

    fig.suptitle(
        "spot width against obliquity: is the anisotropy the sensor?", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]))
