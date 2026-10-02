"""mxi_import reads an NXmx master's geometry right, and its overrides do
what they say.

A master is planted with h5py -- metadata only, since mxi_import reads no
pixels -- awkwardly: lengths in metres, a detector chain with a two-theta
rotation in it, a fixed goniometer axis at an angle, a relative depends_on.
What the experiment list must hold is worked out here, independently, with
numpy's own rotations, and the overrides are each held to what they ask for.
And, where the data is to hand, the insulin master's mu, goniometer and scan
against dials.import's own output. Needs MXI_IMPORT.
"""

import json
import os
import subprocess

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

IMPORT = os.environ.get("MXI_IMPORT")
needs_import = pytest.mark.skipif(not IMPORT, reason="set MXI_IMPORT")


def rotation(axis, degrees):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    t = np.radians(degrees)
    k = np.array(
        [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
    )
    return np.eye(3) + np.sin(t) * k + (1 - np.cos(t)) * (k @ k)


def imgcif(v):
    return np.array([-v[0], v[1], -v[2]])


TWO_THETA = 15.0
DISTANCE_M = 0.25
OFFSET_M = 0.12
OFFSET_VEC = np.array([0.6, 0.8, 0.0])
CHI = 20.0
PIVOT_M = np.array([0.0, 0.005, 0.002])  # the two-theta arm's offset: not at the sample


def plant(path):
    with h5py.File(path, "w") as f:
        e = f.create_group("entry")
        e.attrs["NX_class"] = "NXentry"
        beam = e.create_group("instrument/beam")
        beam.create_dataset("incident_wavelength", data=1.0).attrs["units"] = "angstrom"
        det = e.create_group("instrument/detector")
        det.attrs["NX_class"] = "NXdetector"
        det.create_dataset("sensor_material", data=b"Silicon")
        det.create_dataset("sensor_thickness", data=0.00045).attrs["units"] = "m"
        det.create_dataset("saturation_value", data=20000)
        det.create_dataset("count_time", data=0.01)
        t = e.create_group("instrument/transformations")
        tt = t.create_dataset("two_theta", data=TWO_THETA)
        tt.attrs.update(
            {
                "transformation_type": "rotation",
                "vector": [-1.0, 0, 0],
                "units": "deg",
                "depends_on": ".",
                "offset": PIVOT_M,
                "offset_units": "m",
            }
        )
        dz = t.create_dataset("det_z", data=DISTANCE_M)
        dz.attrs.update(
            {
                "transformation_type": "translation",
                "vector": [0, 0, 1.0],
                "units": "m",
                "depends_on": "two_theta",
            }
        )
        mod = det.create_group("module")
        mod.attrs["NX_class"] = "NXdetector_module"
        mod.create_dataset("data_origin", data=[0, 0])
        mod.create_dataset("data_size", data=[300, 200])
        mo = mod.create_dataset("module_offset", data=OFFSET_M)
        mo.attrs.update(
            {
                "transformation_type": "translation",
                "vector": OFFSET_VEC,
                "units": "m",
                "depends_on": "/entry/instrument/transformations/det_z",
            }
        )
        for name, vec in (
            ("fast_pixel_direction", [-1.0, 0, 0]),
            ("slow_pixel_direction", [0, -1.0, 0]),
        ):
            d = mod.create_dataset(name, data=7.5e-5)
            d.attrs.update(
                {
                    "transformation_type": "translation",
                    "vector": vec,
                    "units": "m",
                    "depends_on": "module_offset",
                }
            )
        s = e.create_group("sample")
        s.create_dataset("depends_on", data=b"/entry/sample/transformations/phi")
        g = s.create_group("transformations")
        om = g.create_dataset("omega", data=np.arange(10) * 0.5 + 10.0)
        om.attrs.update(
            {
                "transformation_type": "rotation",
                "vector": [-1.0, 0, 0],
                "units": "deg",
                "depends_on": ".",
            }
        )
        chi = g.create_dataset("chi", data=[CHI])
        chi.attrs.update(
            {
                "transformation_type": "rotation",
                "vector": [0, 0, 1.0],
                "units": "deg",
                "depends_on": "omega",
            }
        )
        phi = g.create_dataset("phi", data=[0.0])
        phi.attrs.update(
            {
                "transformation_type": "rotation",
                "vector": [-1.0, -0.002, 0.001],
                "units": "deg",
                "depends_on": "chi",
            }
        )


def expected_panel():
    # p -> two_theta(det_z(module_offset(p))), McStas, then imgCIF; a
    # transformation with an offset is R p + offset, the offset after.
    r = rotation([-1, 0, 0], TWO_THETA)
    inner = (
        np.array([0, 0, DISTANCE_M * 1000])
        + OFFSET_VEC / np.linalg.norm(OFFSET_VEC) * OFFSET_M * 1000
    )
    origin = r @ inner + PIVOT_M * 1000
    return imgcif(origin), imgcif(r @ [-1.0, 0, 0]), imgcif(r @ [0, -1.0, 0])


def run(tmp_path, *extra):
    plant(tmp_path / "master.nxs")
    result = subprocess.run(
        [
            IMPORT,
            str(tmp_path / "master.nxs"),
            "-o",
            str(tmp_path / "imported.expt"),
            *extra,
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    return json.load(open(tmp_path / "imported.expt")), result.stdout


@needs_import
def test_a_planted_masters_geometry_is_read_as_worked_out_here(tmp_path):
    e, _ = run(tmp_path)
    p = e["detector"][0]["panels"][0]
    origin, fast, slow = expected_panel()
    assert np.allclose(p["origin"], origin, atol=1e-9)
    assert np.allclose(p["fast_axis"], fast, atol=1e-12)
    assert np.allclose(p["slow_axis"], slow, atol=1e-12)
    assert p["image_size"] == [200, 300] and np.allclose(
        p["pixel_size"], [0.075, 0.075]
    )
    assert p["trusted_range"] == [0.0, 20000.0] and p["thickness"] == pytest.approx(
        0.45
    )
    assert p["material"] == "Si" and p["mu"] > 0
    g = e["goniometer"][0]
    assert g["names"] == ["phi", "chi", "omega"] and g["scan_axis"] == 2
    assert np.allclose(g["axes"][0], imgcif([-1.0, -0.002, 0.001]))
    assert np.allclose(g["axes"][1], imgcif([0, 0, 1.0]))
    assert g["angles"] == [0.0, CHI, 0.0]
    s = e["scan"][0]
    assert s["image_range"] == [1, 10]
    assert np.allclose(s["properties"]["oscillation"], np.arange(10) * 0.5 + 10.0)
    assert np.allclose(s["properties"]["exposure_time"], 0.01)
    assert e["beam"][0]["wavelength"] == pytest.approx(1.0)


def beam_meets(panel):
    o, f, s = (
        np.asarray(panel[k], float) for k in ("origin", "fast_axis", "slow_axis")
    )
    n = np.cross(f, s)
    t = o.dot(n) / np.array([0, 0, -1.0]).dot(n)
    hit = np.array([0, 0, -t])
    return (
        (hit - o).dot(f) / panel["pixel_size"][0],
        (hit - o).dot(s) / panel["pixel_size"][1],
        abs(o.dot(n / np.linalg.norm(n))),
    )


@needs_import
def test_the_overrides_do_what_they_say(tmp_path):
    e, out = run(
        tmp_path,
        "--beam-centre",
        "123.5,77.25",
        "--distance",
        "180",
        "--wavelength",
        "0.8",
        "--image-range",
        "3,7",
    )
    p = e["detector"][0]["panels"][0]
    bx, by, d = beam_meets(p)
    assert (bx, by) == pytest.approx((123.5, 77.25), abs=1e-6)
    assert d == pytest.approx(180.0, abs=1e-6)
    _, fast, slow = expected_panel()
    assert np.allclose(p["fast_axis"], fast) and np.allclose(
        p["slow_axis"], slow
    )  # moved, not turned
    assert e["beam"][0]["wavelength"] == pytest.approx(0.8)
    e0, _ = run(tmp_path)
    assert (
        p["mu"] < e0["detector"][0]["panels"][0]["mu"]
    )  # shorter wavelength, less absorption
    s = e["scan"][0]
    assert s["image_range"] == [3, 7] and len(s["properties"]["oscillation"]) == 5
    assert "the beam centre given" in out and "the distance given" in out


@needs_import
def test_the_insulin_master_against_dials_import(tmp_path):
    master = os.environ.get("MXI_TEST_IMAGES")
    reference = os.environ.get("MXI_DIALS_IMPORTED")
    if not (
        master and reference and os.path.exists(master) and os.path.exists(reference)
    ):
        pytest.skip("set MXI_TEST_IMAGES and MXI_DIALS_IMPORTED")
    subprocess.run(
        [IMPORT, master, "-o", str(tmp_path / "ours.expt")],
        check=True,
        capture_output=True,
    )
    ours = json.load(open(tmp_path / "ours.expt"))
    theirs = json.load(open(reference))
    assert ours["goniometer"] == theirs["goniometer"]
    p, q = ours["detector"][0]["panels"][0], theirs["detector"][0]["panels"][0]
    assert p["mu"] == pytest.approx(q["mu"], rel=1e-12)
    for k in (
        "fast_axis",
        "slow_axis",
        "pixel_size",
        "trusted_range",
        "thickness",
        "material",
    ):
        assert p[k] == pytest.approx(q[k]) if k != "material" else p[k] == q[k], k
    assert ours["scan"][0]["image_range"] == theirs["scan"][0]["image_range"]
    assert np.allclose(
        ours["scan"][0]["properties"]["oscillation"],
        theirs["scan"][0]["properties"]["oscillation"],
    )
