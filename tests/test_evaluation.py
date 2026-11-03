from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from calder.evaluation.calibration import (
    brier_score,
    expected_calibration_error,
    grouped_calibration,
    predictive_interval_coverage,
    reliability_bins,
)
from calder.evaluation.coverage import confirmed_group_coverage, expected_new_group_coverage
from calder.evaluation.hit_metrics import area_under_discovery, discovery_metrics
from calder.evaluation.leakage import assert_disjoint, assert_same_decision, assert_unique_selection
from calder.evaluation.reports import (
    _model_noise_from_config,
    build_context_comparison,
    build_report,
    calibration_summary,
    detailed_calibration_tables,
)
from calder.evaluation.statistics import paired_bootstrap, paired_policy_difference


def test_discovery_denominators_and_step_area():
    memberships = np.eye(3)
    metrics = discovery_metrics([0, 1], [2, 0], 1, memberships, total_hit_count=2, costs=[1, 2, 3])
    assert metrics == {
        "confirmed_hits": 1.0,
        "precision": 0.5,
        "recall": 0.5,
        "group_coverage": 1.0,
        "experiments": 2.0,
        "budget_spent": 3.0,
    }
    assert area_under_discovery([0, 2, 5], [1, 3, 6]) == 11
    with pytest.raises(ValueError, match="twice"):
        discovery_metrics([0, 0], [2, 2], 1, memberships, 2)


def test_coverage_independence_formula_and_missing_groups():
    memberships = np.array([[1.0, 0.0], [1.0, 1.0]])
    assert expected_new_group_coverage([0.5, 0.5], memberships) == pytest.approx(1.25)
    assert expected_new_group_coverage([0.5, 0.5], memberships, [True, False]) == 0.5
    assert confirmed_group_coverage(memberships, [False, True]) == 2
    assert confirmed_group_coverage(np.zeros((2, 0)), [True, True]) == 0


def test_calibration_exact_small_examples_and_boundaries():
    assert brier_score([0, 1], [0, 1]) == 0
    assert brier_score([0.2, 0.8], [0, 1]) == pytest.approx(0.04)
    assert expected_calibration_error([0.2, 0.8], [0, 1], 2) == pytest.approx(0.2)
    bins = reliability_bins([0, 1], [0, 1], 2)
    assert bins["count"].tolist() == [1, 1]
    assert np.isnan(brier_score([], []))
    assert grouped_calibration([0.5], [1], [[1, 0]]).iloc[1]["n"] == 0
    with pytest.raises(ValueError):
        brier_score([1.1], [1])
    with pytest.raises(ValueError):
        brier_score([0.5], [0.5])


def test_interval_targets_and_zero_variance():
    assert predictive_interval_coverage([0, 0], [0, 1], [0, 2]) == 0.5
    with pytest.raises(ValueError):
        predictive_interval_coverage([0], [-1], [0])


def test_measured_predictive_interval_includes_likelihood_floor_only_for_measured_targets(tmp_path):
    frame = pd.DataFrame(
        {
            "policy": ["calder"],
            "seed": [0],
            "probability": [0.8],
            "is_hit": [True],
            "latent_label": [True],
            "mean": [0.0],
            "observed_effect": [0.3],
            "latent_effect": [0.3],
            "noise_variance": [0.0],
            "epistemic_variance": [0.0],
        }
    )
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"model": {"noise_floor": 0.04}}))
    restored = _model_noise_from_config(frame, config)
    assert "model_noise_variance" not in frame
    assert restored.model_noise_variance.iloc[0] == 0.04
    scores = calibration_summary(restored).set_index("label_target")
    assert scores.loc["measured", "interval_target"] == "measured_predictive_model_assumed"
    assert scores.loc["measured", "interval_coverage_95"] == 1.0
    assert scores.loc["latent", "interval_coverage_95"] == 0.0


def test_paired_bootstrap_preserves_pairing_and_rejects_ambiguous_rows():
    result = paired_bootstrap([10, 20, 30], [8, 18, 28], seed=1)
    assert result == {"mean_difference": 2.0, "ci_lower": 2.0, "ci_upper": 2.0, "n_pairs": 3}
    frame = pd.DataFrame({"policy": ["a", "b"], "seed": [0, 1], "hits": [4, 3]})
    with pytest.raises(ValueError, match="identical paired"):
        paired_policy_difference(frame, "a", "b", "hits")
    duplicate = pd.DataFrame({"policy": ["a", "a", "b"], "seed": [0, 0, 0], "hits": [4, 4, 3]})
    with pytest.raises(ValueError, match="Duplicate"):
        paired_policy_difference(duplicate, "a", "b", "hits")


def test_leakage_assertions_fail_loudly():
    with pytest.raises(ValueError, match="Leakage"):
        assert_disjoint([1, 2], [2, 3])
    with pytest.raises(ValueError, match="repeatedly"):
        assert_unique_selection([1], [1])
    with pytest.raises(AssertionError, match="Hidden-label"):
        assert_same_decision([1], [2])


def test_report_uses_actual_outputs_and_marks_missing_biology(tmp_path):
    rows = []
    for seed in [0, 1]:
        for policy, hits in [("probability_of_hit", 3 + seed), ("calder", 2 + seed)]:
            for round_number in [0, 1]:
                rows.append(
                    {
                        "dataset": "synthetic",
                        "scenario": "A",
                        "threshold_quantile": 0.8,
                        "seed": seed,
                        "policy": policy,
                        "round": round_number,
                        "additional_budget": round_number * 2,
                        "budget_spent": 2 + round_number * 2,
                        "confirmed_hits": round_number * hits,
                        "group_coverage": round_number,
                        "auc_hits": 0.0,
                    }
                )
    pd.DataFrame(rows).to_csv(tmp_path / "rounds.csv", index=False)
    (tmp_path / "experiment_manifest.json").write_text(
        json.dumps({"config_hash": "abc", "git_commit": None})
    )
    result = build_report(tmp_path)
    assert (tmp_path / "results.csv").exists()
    differences = pd.read_csv(tmp_path / "paired_differences.csv")
    row = differences[
        (differences.policy == "calder") & (differences.metric == "confirmed_hits")
    ].iloc[0]
    assert row.mean_difference == -1
    assert result["figure_status"]["figure_6"]["status"] == "unavailable"
    assert "controlled synthetic simulation" in (tmp_path / "summary.md").read_text()
    report = (tmp_path / "research_report.md").read_text()
    assert "confirmed-hit mean differences ranged from -1 to -1" in report
    assert "Biological K562–RPE1 performance remains unmeasured" in report
    manifest = json.loads((tmp_path / "experiment_manifest.json").read_text())
    assert manifest["config_hash"] == "abc"
    assert manifest["git_commit"] is None
    assert "source_tree_hash" in manifest["report_generation_provenance"]
    from calder.utils.serialization import file_sha256

    assert manifest["files"]["research_report.md"] == file_sha256(tmp_path / "research_report.md")


def test_calibration_report_distinguishes_measured_and_latent_targets():
    data = pd.DataFrame(
        {
            "policy": ["calder", "calder"],
            "seed": [0, 0],
            "probability": [0.8, 0.2],
            "raw_probability": [0.9, 0.1],
            "is_hit": [1, 0],
            "latent_label": [0, 1],
            "mean": [0, 1],
            "true_log_effect": [0, 0.5],
            "observed_effect": [0.1, 1.2],
            "epistemic_variance": [0.1, 0.2],
            "noise_variance": [0.3, 0.3],
        }
    )
    result = calibration_summary(data)
    assert set(result.label_target) == {"latent", "measured"}
    raw = result[result.probability_type == "raw_probability"]
    assert (
        raw.loc[raw.label_target == "latent", "brier_score"].iloc[0]
        > raw.loc[raw.label_target == "measured", "brier_score"].iloc[0]
    )


def test_saved_group_calibration_and_bins_keep_empty_bins():
    data = pd.DataFrame(
        {
            "policy": ["calder", "calder"],
            "seed": [0, 0],
            "probability": [0.8, 0.2],
            "raw_probability": [0.9, 0.1],
            "calibrated_probability": [0.8, 0.2],
            "is_hit": [True, False],
            "groups": ["[0, 1]", "[1]"],
        }
    )
    bins, groups = detailed_calibration_tables(data)
    assert bins["count"].sum() == 6
    assert len(bins) == 30
    assert set(groups.group) == {"0", "1"}
    assert set(groups.probability_type) == {
        "probability",
        "raw_probability",
        "calibrated_probability",
    }


def test_context_comparison_rejects_retuned_hyperparameters(tmp_path):
    left, right = tmp_path / "k562", tmp_path / "rpe1"
    for path, scale in [(left, 1), (right, 2)]:
        path.mkdir()
        (path / "config.json").write_text(
            json.dumps({"model": {"length_scale": scale}, "acquisition": {}, "protocol": {}})
        )
    with pytest.raises(ValueError, match="frozen model"):
        build_context_comparison(left, right, tmp_path / "comparison")
