"""Pictures of the reference profiles, and of where two integrations differ.

Numbers in a table say a disagreement trends with something.  They do not show
its shape, and a profile is a shape: whether it is off-centre, lopsided,
broader in one direction than the model allows, or has a hole in it where a
module gap fell.  Those look alike in a correlation and not at all alike on a
page.

matplotlib is optional everywhere else in this package and is optional here:
nothing in the checker imports this module.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class ReferenceProfiles:
    """What `mxi_integrate --save-profiles` wrote."""

    side: int
    sigma_b: float
    sigma_m: float
    half_width: float
    divisions: int
    blocks: int
    panels: int
    profiles: list[np.ndarray] = field(default_factory=list)
    spots: list[int] = field(default_factory=list)

    def grid(self, which: int) -> np.ndarray:
        """One profile as (e3, e2, e1), which is the order it was written in."""
        n = self.side
        return self.profiles[which].reshape(n, n, n)


def read_profiles(path: str) -> ReferenceProfiles:
    header: dict[str, str] = {}
    profiles: list[list[float]] = []
    spots: list[int] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0] == "profile":
                profiles.append([])
                spots.append(int(parts[3]))
            elif len(parts) == 2 and not profiles:
                header[parts[0]] = parts[1]
            elif profiles:
                try:
                    profiles[-1].append(float(parts[0]))
                except ValueError as exc:
                    raise ValueError(
                        f"{path}: expected a number in profile "
                        f"{len(profiles) - 1}, found {line!r}"
                    ) from exc
            else:
                raise ValueError(
                    f"{path}: found {line!r} before any 'profile' line, so the "
                    "header and the body disagree about where one starts"
                )
    out = ReferenceProfiles(
        side=int(header["side"]),
        sigma_b=float(header["sigma_b"]),
        sigma_m=float(header["sigma_m"]),
        half_width=float(header["half_width"]),
        divisions=int(header["divisions"]),
        # Older files have no scan blocks; one block is what they meant.
        blocks=int(header.get("blocks", 1)),
        panels=int(header["panels"]),
        profiles=[np.asarray(p, dtype=float) for p in profiles],
        spots=spots,
    )
    expected = out.side**3
    for i, p in enumerate(out.profiles):
        if p.size != expected:
            raise ValueError(
                f"profile {i} has {p.size} points, not {expected}: the file and "
                "the side length disagree"
            )
    return out


def draw_profiles(
    reference: ReferenceProfiles, path: str, block: int | None = None
) -> str:
    """Every region's profile, as three central sections each.

    Sections rather than projections: a projection hides a profile that is
    hollow or double-peaked, which is exactly the sort of thing worth seeing.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # With a divided scan there can be forty-five profiles, which is a page
    # nobody reads. One block at a time is the useful view.
    per_panel = reference.divisions * reference.divisions
    if block is None:
        which = list(range(len(reference.profiles)))
    else:
        which = [
            r
            for r in range(len(reference.profiles))
            if (r // per_panel) % max(reference.blocks, 1) == block
        ]
        if not which:
            raise ValueError(f"no block {block}: the file has {reference.blocks}")
    n = len(which)
    middle = reference.side // 2
    fig, axes = plt.subplots(n, 3, figsize=(7.5, 2.5 * n), squeeze=False)
    for row, r in enumerate(which):
        grid = reference.grid(r)
        sections = [
            ("e1 e2", grid[middle, :, :]),
            ("e1 e3", grid[:, middle, :]),
            ("e2 e3", grid[:, :, middle]),
        ]
        for c, (name, plane) in enumerate(sections):
            ax = axes[row][c]
            ax.imshow(plane, origin="lower", interpolation="nearest")
            ax.set_xticks([])
            ax.set_yticks([])
            if row == 0:
                ax.set_title(name, fontsize=9)
            if c == 0:
                per_panel = reference.divisions * reference.divisions
                block = (r // per_panel) % max(reference.blocks, 1)
                cell = r % per_panel
                label = (
                    f"block {block}\ncell {cell}"
                    if reference.blocks > 1
                    else f"region {r}"
                )
                ax.set_ylabel(f"{label}\n{reference.spots[r]} spots", fontsize=7)
    fig.suptitle("reference profiles, central sections", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def draw_profile_widths(reference: ReferenceProfiles, path: str) -> str:
    """The second moment of each region's profile, along each axis.

    A profile that is wider in one region than another is the detector
    changing across its face; a profile wider in `e3` than the model allows is
    the rocking curve, not the detector.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(reference.profiles)
    axis_names = ("e1", "e2", "e3")
    widths = np.zeros((n, 3))
    coord = np.arange(reference.side) - reference.side // 2
    for r in range(n):
        grid = reference.grid(r)
        total = grid.sum()
        if total <= 0:
            continue
        for axis in range(3):
            marginal = grid.sum(axis=tuple(i for i in range(3) if i != axis))
            mean = (marginal * coord).sum() / total
            widths[r, 2 - axis] = np.sqrt(
                max((marginal * (coord - mean) ** 2).sum() / total, 0.0)
            )
    fig, ax = plt.subplots(figsize=(6.5, 3.5))
    for axis in range(3):
        ax.plot(np.arange(n), widths[:, axis], "o-", label=axis_names[axis])
    ax.set_xlabel("region")
    ax.set_ylabel("width, grid points")
    ax.set_title("profile width by region and axis")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path
