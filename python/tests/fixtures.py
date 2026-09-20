"""Synthetic reflection tables and experiment lists with known answers.

These are not imitations of real data and are not trying to be.  Each fixture
is built so that the answer to the check being tested is known exactly by
construction: a planted reindexing operator, a planted intensity ratio, a
planted misorientation.  A test that asserts against a number nobody chose is
a test that passes for reasons nobody knows.
"""

from __future__ import annotations

import json

import numpy as np

from mxeq.refl import ReflectionTable

# Insulin-like: cubic I 2 3, a = 78 A.  Used because the numbers are familiar
# enough that a wrong one looks wrong.
CELL = 78.0
HALL = " I 2 2 3"


def experiments_dict(
    cell: float = CELL,
    hall: str = HALL,
    image_range: tuple[int, int] = (1, 600),
    oscillation: tuple[float, float] = (0.0, 0.1),
    rotation: np.ndarray | None = None,
    n_scan_points: int = 0,
) -> dict:
    """A minimal but complete experiment list document."""
    a = np.eye(3) * cell
    if rotation is not None:
        a = a @ np.asarray(rotation, dtype=float).T

    crystal = {
        "__id__": "crystal",
        "real_space_a": list(a[0]),
        "real_space_b": list(a[1]),
        "real_space_c": list(a[2]),
        "space_group_hall_symbol": hall,
    }
    if n_scan_points:
        # A small monotonic cell expansion over the scan, which is what
        # radiation damage looks like and gives the scan-varying branch
        # something with a known sign to find.
        setting = []
        for i in range(n_scan_points):
            scale = 1.0 + 1e-4 * i / max(n_scan_points - 1, 1)
            setting.append(list(np.linalg.inv(a * scale).flatten()))
        crystal["A_at_scan_points"] = setting

    return {
        "__id__": "ExperimentList",
        "experiment": [
            {
                "__id__": "Experiment",
                "identifier": "fixture",
                "beam": 0,
                "detector": 0,
                "crystal": 0,
                "scan": 0,
                "goniometer": 0,
                "imageset": -1,
            }
        ],
        "crystal": [crystal],
        "beam": [{"direction": [0.0, 0.0, 1.0], "wavelength": 0.9537}],
        "detector": [
            {
                "panels": [
                    {
                        "name": "panel0",
                        "fast_axis": [1.0, 0.0, 0.0],
                        "slow_axis": [0.0, -1.0, 0.0],
                        "origin": [-155.0, 155.0, -200.0],
                        "image_size": [4148, 4362],
                        "pixel_size": [0.075, 0.075],
                        "trusted_range": [0.0, 65535.0],
                    }
                ]
            }
        ],
        "goniometer": [
            {
                "rotation_axis": [1.0, 0.0, 0.0],
                "fixed_rotation": [1, 0, 0, 0, 1, 0, 0, 0, 1],
                "setting_rotation": [1, 0, 0, 0, 1, 0, 0, 0, 1],
            }
        ],
        "scan": [
            {
                "image_range": list(image_range),
                "oscillation": list(oscillation),
                "batch_offset": 0,
            }
        ],
    }


def write_experiments(path, **kwargs) -> str:
    with open(path, "w") as f:
        json.dump(experiments_dict(**kwargs), f)
    return str(path)


def strong_table(n: int = 500, seed: int = 1) -> ReflectionTable:
    rng = np.random.default_rng(seed)
    xyz = np.column_stack(
        [
            rng.uniform(50, 4000, n),
            rng.uniform(50, 4200, n),
            rng.uniform(0, 600, n),
        ]
    )
    # Long-tailed intensities, because the weak end is where two spot finders
    # differ and a uniform distribution would hide that.
    intensity = np.exp(rng.normal(5.0, 1.5, n))
    n_signal = np.maximum(3, (intensity / 40).astype(np.int32))
    depth = rng.integers(1, 4, n).astype(np.int32)
    bbox = np.column_stack(
        [
            (xyz[:, 0] - 3).astype(np.int32),
            (xyz[:, 0] + 3).astype(np.int32),
            (xyz[:, 1] - 3).astype(np.int32),
            (xyz[:, 1] + 3).astype(np.int32),
            xyz[:, 2].astype(np.int32),
            xyz[:, 2].astype(np.int32) + depth,
        ]
    )
    table = ReflectionTable(nrows=n)
    table.columns["xyzobs.px.value"] = xyz
    table.types["xyzobs.px.value"] = "vec3<double>"
    table.columns["intensity.sum.value"] = intensity
    table.types["intensity.sum.value"] = "double"
    table.columns["intensity.sum.variance"] = intensity + 25.0
    table.types["intensity.sum.variance"] = "double"
    table.columns["n_signal"] = n_signal
    table.types["n_signal"] = "int"
    table.columns["bbox"] = bbox
    table.types["bbox"] = "int6"
    table.columns["panel"] = np.zeros(n, dtype=np.uint64)
    table.types["panel"] = "std::size_t"
    table.columns["flags"] = np.full(n, 32, dtype=np.uint64)
    table.types["flags"] = "std::size_t"
    table.identifiers = {0: "fixture"}
    return table


def _unique_hkl(n: int, rng: np.random.Generator, limit: int = 30) -> np.ndarray:
    """``n`` distinct non-zero Miller indices."""
    pool = set()
    while len(pool) < n:
        candidate = rng.integers(-limit, limit + 1, (n, 3))
        for h in candidate:
            if h.any():
                pool.add(tuple(int(v) for v in h))
            if len(pool) >= n:
                break
    return np.array(sorted(pool)[:n], dtype=np.int32)


def integrated_table(
    n: int = 800, seed: int = 2, cell: float = CELL
) -> ReflectionTable:
    rng = np.random.default_rng(seed)
    hkl = _unique_hkl(n, rng)
    setting = np.linalg.inv(np.eye(3) * cell)
    d = 1.0 / np.linalg.norm(hkl @ setting.T, axis=1)

    # Intensity falling off with resolution, so the shell tables have a
    # gradient in them and a shell-dependent bug would be visible.
    intensity = np.exp(rng.normal(6.0, 1.2, n)) * np.exp(-2.0 / np.maximum(d, 0.5) ** 2)
    variance = np.maximum(intensity, 1.0) + 30.0

    table = ReflectionTable(nrows=n)

    def add(name, values, type_name):
        table.columns[name] = values
        table.types[name] = type_name

    add("miller_index", hkl.astype(np.int32), "cctbx::miller::index<>")
    add("entering", rng.integers(0, 2, n).astype(bool), "bool")
    add("id", np.zeros(n, dtype=np.int32), "int")
    add("intensity.sum.value", intensity, "double")
    add("intensity.sum.variance", variance, "double")
    add("intensity.prf.value", intensity * rng.normal(1.0, 0.02, n), "double")
    add("intensity.prf.variance", variance * 0.9, "double")
    add("partiality", np.clip(rng.beta(8, 2, n), 0.05, 1.0), "double")
    add("lp", rng.uniform(0.5, 1.5, n), "double")
    add("d", d, "double")
    add("background.mean", rng.uniform(0.2, 3.0, n), "double")
    add(
        "xyzcal.px",
        np.column_stack(
            [
                rng.uniform(50, 4000, n),
                rng.uniform(50, 4200, n),
                rng.uniform(0, 600, n),
            ]
        ),
        "vec3<double>",
    )
    add("inverse_scale_factor", rng.uniform(0.8, 1.2, n), "double")
    add("intensity.scale.value", intensity, "double")
    add("intensity.scale.variance", variance, "double")
    table.identifiers = {0: "fixture"}
    return table


def perturb_positions(
    table: ReflectionTable, sigma: float = 0.05, seed: int = 3
) -> ReflectionTable:
    rng = np.random.default_rng(seed)
    out = table.select(np.arange(table.nrows))
    out.columns["xyzobs.px.value"] = table["xyzobs.px.value"] + rng.normal(
        0.0, sigma, table["xyzobs.px.value"].shape
    )
    return out


def drop(table: ReflectionTable, fraction: float, seed: int = 4) -> ReflectionTable:
    rng = np.random.default_rng(seed)
    keep = rng.random(table.nrows) >= fraction
    return table.select(keep)


def reindex(table: ReflectionTable, operator: np.ndarray) -> ReflectionTable:
    out = table.select(np.arange(table.nrows))
    out.columns["miller_index"] = (
        table["miller_index"].astype(np.int64) @ np.asarray(operator, dtype=np.int64)
    ).astype(np.int32)
    return out
