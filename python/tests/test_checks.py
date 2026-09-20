"""The checks, each against a perturbation whose answer is known in advance."""

from __future__ import annotations

import json

import numpy as np
import pytest

import fixtures
from mxeq import expt
from mxeq.checks import indexed, integrated, refined, scaled, strong
from mxeq.checks.refined import misorientation, rotation_angle


def values(report, section_title):
    for section in report.sections:
        if section.title == section_title:
            return section.values
    raise AssertionError(
        f"no section {section_title!r}; have {[s.title for s in report.sections]}"
    )


# --------------------------------------------------------------------------
# strong
# --------------------------------------------------------------------------


def test_strong_identical_tables_match_completely():
    table = fixtures.strong_table(n=300)
    report = strong.check(table, table)
    v = values(report, "spatial match")
    assert v["n_matched"] == 300
    assert v["fraction_matched"] == pytest.approx(1.0)
    assert not report.warnings
    offsets = values(report, "centroid offset, A minus B")
    assert offsets["offset x"]["median"] == pytest.approx(0.0)


def test_strong_finds_a_planted_centroid_shift():
    a = fixtures.strong_table(n=300)
    b = a.select(np.arange(a.nrows))
    b.columns["xyzobs.px.value"] = a["xyzobs.px.value"] + np.array([0.25, 0.0, 0.0])
    v = values(strong.check(a, b), "centroid offset, A minus B")
    assert v["offset x"]["median"] == pytest.approx(-0.25, abs=1e-9)
    assert v["offset y"]["median"] == pytest.approx(0.0, abs=1e-9)


def test_strong_counts_dropped_spots_as_unmatched():
    a = fixtures.strong_table(n=400)
    b = fixtures.drop(a, fraction=0.1, seed=7)
    v = values(strong.check(a, b), "spatial match")
    assert v["n_matched"] == b.nrows
    assert v["n_only_a"] == a.nrows - b.nrows
    assert v["n_only_b"] == 0


def test_strong_reports_the_strength_of_what_did_not_match():
    """The unmatched population must be describable, not just counted."""
    a = fixtures.strong_table(n=400)
    order = np.argsort(a["intensity.sum.value"])
    # Drop the hundred weakest, which is how two spot finders normally differ.
    b = a.select(order[100:])
    v = values(strong.check(a, b), "unmatched, by strength")
    weak = v["A intensity.sum.value unmatched"]["median"]
    kept = v["A intensity.sum.value matched"]["median"]
    assert weak < kept


def test_strong_warns_when_nothing_matches():
    a = fixtures.strong_table(n=50)
    b = a.select(np.arange(a.nrows))
    b.columns["xyzobs.px.value"] = a["xyzobs.px.value"] + 500.0
    report = strong.check(a, b)
    assert report.warnings
    assert values(report, "spatial match")["n_matched"] == 0


def test_strong_reports_spot_depth():
    table = fixtures.strong_table(n=200)
    v = values(strong.check(table, table), "spot depth in images, A")
    assert v["depth"]["n"] == 200
    assert 0.0 <= v["fraction_single_image"] <= 1.0


# --------------------------------------------------------------------------
# indexed
# --------------------------------------------------------------------------


def _indexable_pair(n=400, seed=11):
    """Two tables describing the same spots, with positions and indices."""
    rng = np.random.default_rng(seed)
    table = fixtures.integrated_table(n=n, seed=seed)
    table.columns["xyzobs.px.value"] = table["xyzcal.px"] + rng.normal(0, 0.02, (n, 3))
    table.types["xyzobs.px.value"] = "vec3<double>"
    return table


def test_indexed_identical_solutions_give_the_identity():
    table = _indexable_pair()
    report = indexed.check(table, table)
    v = values(report, "reindexing")
    assert v["is_identity"] is True
    assert v["agreement_best"] == pytest.approx(1.0)
    assert values(report, "index agreement after reindexing")["n_disagree"] == 0


def test_indexed_recovers_a_planted_reindexing():
    operator = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=np.int64)
    a = _indexable_pair()
    b = fixtures.reindex(a, operator)
    report = indexed.check(a, b)
    v = values(report, "reindexing")
    assert v["is_identity"] is False
    assert np.array_equal(np.array(v["operator"]), operator)
    assert v["agreement_best"] == pytest.approx(1.0)
    assert v["agreement_identity"] < 0.05
    assert values(report, "index agreement after reindexing")["n_disagree"] == 0


def test_indexed_reports_disagreement_when_indices_really_differ():
    a = _indexable_pair()
    b = a.select(np.arange(a.nrows))
    hkl = b["miller_index"].copy()
    hkl[:40, 0] += 1  # near misses, one unit in h
    b.columns["miller_index"] = hkl
    report = indexed.check(a, b)
    v = values(report, "index agreement after reindexing")
    assert v["n_disagree"] == 40
    assert v["disagreement |dh| max"]["median"] == pytest.approx(1.0)


def test_indexed_counts_partially_indexed_tables():
    a = _indexable_pair()
    b = a.select(np.arange(a.nrows))
    hkl = b["miller_index"].copy()
    hkl[:100] = 0
    b.columns["miller_index"] = hkl
    report = indexed.check(a, b)
    assert values(report, "how much was indexed")["n_indexed_b"] == a.nrows - 100
    assert values(report, "indexed on one side only")["n_indexed_only_a"] == 100


def test_indexed_flags_a_metric_incompatible_operator(tmp_path):
    """A tetragonal cell cannot admit a-c exchange, and must say so."""
    path = fixtures.write_experiments(tmp_path / "tet.expt")
    experiments = expt.load(path)
    experiments[0].crystal.a = np.diag([50.0, 50.0, 120.0])
    swap_ac = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]], dtype=np.int64)
    a = _indexable_pair()
    b = fixtures.reindex(a, swap_ac)
    report = indexed.check(a, b, experiments_a=experiments)
    # Reported, but no longer a warning on its own: an operator that explains
    # every pair is a change of basis whether or not it is a lattice symmetry,
    # and two different reduced cells of one lattice are related by exactly
    # such an operator. The warning is reserved for one that explains neither.
    assert values(report, "reindexing")["is_a_lattice_symmetry"] is False
    assert values(report, "reindexing")["agreement_best"] > 0.9
    assert not report.warnings


# --------------------------------------------------------------------------
# refined
# --------------------------------------------------------------------------


def test_rotation_angle_of_a_known_rotation():
    for degrees in (0.0, 1.0, 30.0, 90.0, 179.0):
        radians = np.radians(degrees)
        r = np.array(
            [
                [np.cos(radians), -np.sin(radians), 0.0],
                [np.sin(radians), np.cos(radians), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        assert rotation_angle(r) == pytest.approx(degrees, abs=1e-6)


def test_misorientation_ignores_a_reindexing():
    """A basis change is not a misorientation and must not be reported as one."""
    a = np.linalg.inv(np.eye(3) * 78.0)
    swap = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    angle, operator = misorientation(a, a @ swap)
    assert angle == pytest.approx(0.0, abs=1e-9)
    assert not np.array_equal(operator, np.eye(3, dtype=int))


def test_misorientation_recovers_a_planted_rotation():
    radians = np.radians(0.7)
    r = np.array(
        [
            [np.cos(radians), -np.sin(radians), 0.0],
            [np.sin(radians), np.cos(radians), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    a = np.linalg.inv(np.eye(3) * 78.0)
    angle, _ = misorientation(a, r @ a)
    assert angle == pytest.approx(0.7, abs=1e-6)


def test_refined_identical_models_agree(tmp_path):
    path = fixtures.write_experiments(tmp_path / "a.expt")
    a = expt.load(path)
    report = refined.check(a, expt.load(path))
    assert values(report, "orientation")["misorientation_deg"] == pytest.approx(
        0.0, abs=1e-9
    )
    assert values(report, "unit cell")["volume_relative_difference"] == pytest.approx(
        0.0
    )
    assert values(report, "detector")["distance_difference_mm"] == pytest.approx(0.0)


def test_refined_recovers_a_planted_cell_change(tmp_path):
    a = expt.load(fixtures.write_experiments(tmp_path / "a.expt", cell=78.0))
    b = expt.load(fixtures.write_experiments(tmp_path / "b.expt", cell=78.1))
    v = values(refined.check(a, b), "unit cell")
    assert v["volume_relative_difference"] == pytest.approx(
        (78.1 / 78.0) ** 3 - 1, rel=1e-6
    )


def test_refined_recovers_a_planted_misorientation(tmp_path):
    radians = np.radians(0.35)
    r = np.array(
        [
            [np.cos(radians), -np.sin(radians), 0.0],
            [np.sin(radians), np.cos(radians), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    a = expt.load(fixtures.write_experiments(tmp_path / "a.expt"))
    b = expt.load(fixtures.write_experiments(tmp_path / "b.expt", rotation=r))
    v = values(refined.check(a, b), "orientation")
    assert v["misorientation_deg"] == pytest.approx(0.35, abs=1e-6)


def test_refined_warns_on_a_different_image_range(tmp_path):
    a = expt.load(fixtures.write_experiments(tmp_path / "a.expt", image_range=(1, 600)))
    b = expt.load(fixtures.write_experiments(tmp_path / "b.expt", image_range=(0, 599)))
    report = refined.check(a, b)
    assert any("image ranges" in w for w in report.warnings)


def test_refined_handles_scan_varying_models(tmp_path):
    a = expt.load(fixtures.write_experiments(tmp_path / "a.expt", n_scan_points=30))
    b = expt.load(fixtures.write_experiments(tmp_path / "b.expt", n_scan_points=30))
    v = values(refined.check(a, b), "scan-varying model")
    assert v["n_scan_points_a"] == 30
    assert v["scan-point misorientation deg"]["max"] == pytest.approx(0.0, abs=1e-7)


def test_refined_survives_mismatched_scan_point_counts(tmp_path):
    a = expt.load(fixtures.write_experiments(tmp_path / "a.expt", n_scan_points=30))
    b = expt.load(fixtures.write_experiments(tmp_path / "b.expt", n_scan_points=12))
    v = values(refined.check(a, b), "scan-varying model")
    assert v["n_scan_points_a"] == 30 and v["n_scan_points_b"] == 12


# --------------------------------------------------------------------------
# integrated
# --------------------------------------------------------------------------


def test_integrated_identical_tables_correlate_perfectly(tmp_path):
    table = fixtures.integrated_table(n=600)
    experiments = expt.load(fixtures.write_experiments(tmp_path / "a.expt"))
    report = integrated.check(table, table, experiments=experiments)
    v = values(report, "intensity.sum.value")
    assert v["cc"] == pytest.approx(1.0)
    assert v["relative difference"]["median"] == pytest.approx(0.0)
    assert v["pull"]["median"] == pytest.approx(0.0)
    assert values(report, "keyed match")["fraction_matched"] == pytest.approx(1.0)


def test_integrated_recovers_a_planted_intensity_ratio(tmp_path):
    a = fixtures.integrated_table(n=600)
    b = a.select(np.arange(a.nrows))
    b.columns["intensity.sum.value"] = a["intensity.sum.value"] * 1.05
    experiments = expt.load(fixtures.write_experiments(tmp_path / "a.expt"))
    v = values(integrated.check(a, b, experiments=experiments), "intensity.sum.value")
    # 2(I - 1.05I)/(I + 1.05I) = -0.05/1.025
    assert v["relative difference"]["median"] == pytest.approx(-0.05 / 1.025, rel=1e-6)
    assert v["cc"] == pytest.approx(1.0)


def test_integrated_join_survives_row_reordering():
    a = fixtures.integrated_table(n=400)
    order = np.random.default_rng(5).permutation(a.nrows)
    b = a.select(order)
    v = values(integrated.check(a, b), "keyed match")
    assert v["n_matched"] == 400


def test_integrated_needs_the_operator_when_the_basis_differs():
    operator = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.int64)
    a = fixtures.integrated_table(n=400)
    b = fixtures.reindex(a, operator)

    without = integrated.check(a, b)
    assert values(without, "keyed match")["fraction_matched"] < 0.05
    assert without.warnings

    with_operator = integrated.check(a, b, operator=operator)
    assert values(with_operator, "keyed match")["fraction_matched"] == pytest.approx(
        1.0
    )


def test_integrated_shell_table_is_present_and_ordered(tmp_path):
    table = fixtures.integrated_table(n=800)
    experiments = expt.load(fixtures.write_experiments(tmp_path / "a.expt"))
    v = values(
        integrated.check(table, table, experiments=experiments), "intensity.sum.value"
    )
    shells = v["intensity.sum.value_shells"]
    assert shells["header"][:3] == ["d_max", "d_min", "n"]
    d_max = [float(row[0]) for row in shells["rows"]]
    assert d_max == sorted(d_max, reverse=True)


def test_integrated_reports_scalar_columns():
    table = fixtures.integrated_table(n=300)
    v = values(integrated.check(table, table), "scalar columns")
    assert "partiality A" in v
    assert v["partiality B-A"]["max"] == pytest.approx(0.0)


# --------------------------------------------------------------------------
# scaled
# --------------------------------------------------------------------------


def test_scaled_identical_data_agree(tmp_path):
    table = fixtures.integrated_table(n=900)
    experiments = expt.load(fixtures.write_experiments(tmp_path / "a.expt"))
    report = scaled.check(table, table, experiments=experiments)
    assert values(report, "unmerged match")["cc_unmerged"] == pytest.approx(1.0)
    v = values(report, "merged A against merged B")
    assert v["cc_merged"] == pytest.approx(1.0)
    assert v["n_only_a"] == 0 and v["n_only_b"] == 0


def test_scaled_reports_merging_statistics(tmp_path):
    table = fixtures.integrated_table(n=900)
    experiments = expt.load(fixtures.write_experiments(tmp_path / "a.expt"))
    report = scaled.check(table, table, experiments=experiments)
    merging = values(report, "merging statistics, each dataset alone")["merging"]
    assert [row[0] for row in merging["rows"]] == ["A", "B"]
    assert int(merging["rows"][0][1]) > 0


def test_scaled_an_overall_scale_is_not_a_disagreement(tmp_path):
    """A constant factor is unobservable after merging and must not show up."""
    a = fixtures.integrated_table(n=900)
    b = a.select(np.arange(a.nrows))
    b.columns["intensity.scale.value"] = a["intensity.scale.value"] * 3.0
    b.columns["intensity.scale.variance"] = a["intensity.scale.variance"] * 9.0
    experiments = expt.load(fixtures.write_experiments(tmp_path / "a.expt"))
    report = scaled.check(a, b, experiments=experiments)
    assert values(report, "merged A against merged B")["cc_merged"] == pytest.approx(
        1.0
    )


def test_scaled_without_a_space_group_stops_after_the_unmerged_part():
    table = fixtures.integrated_table(n=200)
    report = scaled.check(table, table)
    assert any("space group" in w for w in report.warnings)
    assert values(report, "unmerged match")["n_matched"] == 200


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------


def test_reports_render_as_text_and_json():
    table = fixtures.strong_table(n=50)
    report = strong.check(table, table)
    assert "mxeq strong" in report.text()
    assert "No thresholds are applied" in report.text()
    decoded = json.loads(report.json())
    assert decoded["boundary"] == "strong"
    assert decoded["sections"]["spatial match"]["n_matched"] == 50
