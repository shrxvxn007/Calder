#!/usr/bin/env python3
"""Run nine fixed development sensitivities and summarize only within-scenario pairs."""

from __future__ import annotations

import argparse
import json
import warnings
from collections import Counter
from pathlib import Path

import pandas as pd

from calder.config import load_config
from calder.evaluation.statistics import paired_policy_difference
from calder.experiments.runner import run_benchmark
from calder.utils.serialization import content_hash, file_sha256

GRID = (
    ("pilot_fraction", (0.1, 0.2, 0.3)),
    ("batch_size", (5, 10, 20)),
    ("noise_scale", (0.05, 0.2, 0.5)),
)


def markdown(frame: pd.DataFrame) -> str:
    shown = frame.copy()
    for column in shown:
        if pd.api.types.is_numeric_dtype(shown[column]):
            shown[column] = shown[column].map(lambda value: f"{value:.4g}")
    rows = [
        "| " + " | ".join(shown.columns) + " |",
        "| " + " | ".join("---" for _ in shown.columns) + " |",
    ]
    rows.extend("| " + " | ".join(map(str, row)) + " |" for row in shown.to_numpy())
    return "\n".join(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("configs/sensitivity.yaml"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/reference_sensitivity"))
    args = parser.parse_args()
    base = load_config(args.config)
    if base.dataset.name != "synthetic" or base.model.name != "random_features":
        raise ValueError("This development grid requires synthetic data and the frozen RFF model")
    if base.protocol.policies != ["probability_of_hit", "calder"]:
        raise ValueError("This comparison requires the declared PoH and CALDER policy pair")
    if (
        base.protocol.pilot_fraction != 0.2
        or base.protocol.batch_size != 10
        or base.synthetic.noise_scale != 0.2
        or base.protocol.seeds != [0, 1, 2]
        or base.protocol.threshold_quantiles != [0.8]
        or base.synthetic.scenarios != ["A", "B", "C", "D"]
        or base.synthetic.n_candidates != 250
        or base.synthetic.n_features != 24
        or base.protocol.include_oracle
    ):
        raise ValueError(
            "Use the declared development base: n250/d24, A–D, seeds0–2, q0.8, pilot0.2/batch10/noise0.2, no oracle"
        )
    root = args.output_dir
    root.mkdir(parents=True, exist_ok=True)
    frames, source_runs, warning_summary = [], [], []
    for axis, values in GRID:
        for value in values:
            variant = f"{axis}_{value:g}"
            config = base.model_copy(deep=True)
            config.experiment_name = f"development_sensitivity_{variant}"
            if axis == "noise_scale":
                config.synthetic.noise_scale = value
            else:
                setattr(config.protocol, axis, value)
            target = root / variant
            expected_hash = content_hash(config.model_dump(mode="json"))
            manifest_path = target / "experiment_manifest.json"
            if not manifest_path.exists():
                print(f"Running {variant}", flush=True)
                with warnings.catch_warnings(record=True) as captured:
                    warnings.simplefilter("always")
                    try:
                        run_benchmark(config, target, make_report=False)
                    finally:
                        target.mkdir(parents=True, exist_ok=True)
                        records = [
                            {
                                "category": item.category.__name__,
                                "message": str(item.message),
                                "filename": item.filename,
                                "lineno": item.lineno,
                            }
                            for item in captured
                        ]
                        (target / "warnings.json").write_text(json.dumps(records, indent=2) + "\n")
                print(f"Finished {variant}: {len(captured)} captured warnings", flush=True)
            manifest = json.loads(manifest_path.read_text())
            if manifest["status"] != "complete" or manifest["config_hash"] != expected_hash:
                raise ValueError(
                    f"Existing run is incomplete or has a different configuration: {target}"
                )
            frame = pd.read_csv(target / "seed_level_results.csv")
            frame["variant"], frame["axis"], frame["factor_value"] = variant, axis, value
            frames.append(frame)
            source_runs.append(
                {
                    "variant": variant,
                    "axis": axis,
                    "factor_value": value,
                    "path": str(target.resolve()),
                    "config_hash": manifest["config_hash"],
                    "git_commit": manifest.get("git_commit"),
                    "git_dirty": manifest.get("git_dirty"),
                    "source_tree_hash": manifest.get("source_tree_hash"),
                    "n_trajectories": manifest["n_trajectories"],
                }
            )
            records = json.loads((target / "warnings.json").read_text())
            for (category, message), count in Counter(
                (r["category"], r["message"]) for r in records
            ).items():
                warning_summary.append(
                    {"variant": variant, "category": category, "message": message, "count": count}
                )
    endpoints = pd.concat(frames, ignore_index=True)
    expected_trajectories = (
        9
        * len(base.synthetic.scenarios)
        * len(base.protocol.seeds)
        * len(base.protocol.threshold_quantiles)
        * 2
    )
    if len(endpoints) != expected_trajectories:
        raise RuntimeError("Unexpected number of sensitivity trajectories")
    pairing = ["variant", "scenario", "seed", "threshold_quantile"]
    if endpoints.groupby(pairing).threshold.nunique().max() != 1:
        raise RuntimeError("Sensitivity policies did not use identical paired thresholds")
    if endpoints.groupby(pairing).additional_budget.nunique().max() != 1:
        raise RuntimeError("Sensitivity policies did not use identical realized budgets")
    # Three independently executed variants reproduce the same base setting.
    scientific_columns = [
        "scenario",
        "seed",
        "threshold_quantile",
        "policy",
        "selection_order",
        "gene_id",
        "observed_effect",
        "is_hit",
    ]
    base_variants = ["pilot_fraction_0.2", "batch_size_10", "noise_scale_0.2"]
    base_selections = [
        pd.read_csv(root / variant / "selections.csv")[scientific_columns]
        .sort_values(scientific_columns[:-2])
        .reset_index(drop=True)
        for variant in base_variants
    ]
    for replicate in base_selections[1:]:
        pd.testing.assert_frame_equal(base_selections[0], replicate)
    endpoints.to_csv(root / "seed_level_results.csv", index=False)
    keys = ["variant", "axis", "factor_value", "scenario", "model", "policy"]
    metrics = [
        "confirmed_hits",
        "active_hits",
        "latent_hits",
        "group_coverage",
        "precision",
        "recall",
        "auc_hits",
        "brier_score",
        "mae",
        "calibration_fitted_fraction",
        "n_selected",
        "additional_budget",
    ]
    summary = endpoints.groupby(keys, sort=True)[metrics].mean().reset_index()
    summary["n_seeds"] = len(base.protocol.seeds)
    summary.to_csv(root / "summary.csv", index=False)
    paired = []
    for identity, part in endpoints.groupby(
        ["variant", "axis", "factor_value", "scenario"], sort=True
    ):
        variant, axis, value, scenario = identity
        for metric in [
            "confirmed_hits",
            "active_hits",
            "latent_hits",
            "group_coverage",
            "auc_hits",
            "brier_score",
        ]:
            contrast = paired_policy_difference(part, "calder", "probability_of_hit", metric)
            paired.append(
                {
                    "variant": variant,
                    "axis": axis,
                    "factor_value": value,
                    "scenario": scenario,
                    "metric": metric,
                    **contrast,
                }
            )
    differences = pd.DataFrame(paired)
    differences.to_csv(root / "paired_differences.csv", index=False)
    warnings_frame = pd.DataFrame(
        warning_summary, columns=["variant", "category", "message", "count"]
    )
    warnings_frame.to_csv(root / "warning_summary.csv", index=False)
    primary_differences = differences[differences.metric == "confirmed_hits"]
    n_negative = int((primary_differences.mean_difference < 0).sum())
    n_zero = int((primary_differences.mean_difference == 0).sum())
    n_positive = int((primary_differences.mean_difference > 0).sum())
    text = [
        "# Exploratory one-factor sensitivity experiments\n",
        f"These are {len(endpoints)} actual controlled synthetic trajectories from nine predeclared development variants. The frozen RFF model, features and acquisition weights are unchanged. Four scenarios, seeds 0–2, threshold quantile 0.8 and the policies Probability-of-Hit and CALDER were used. No main or biological hyperparameter was selected from these outcomes.\n",
        "The one-factor grid varies pilot fraction (0.1, 0.2, 0.3), maximum batch size (5, 10, 20), and simulated observation-noise scale (0.05, 0.2, 0.5). All other fields retain the supplied base config. The base setting appears once in each factor axis and is intentionally retained, not counted as an independent population. Its three independent reruns produced exactly identical scientific selections and observed outcomes, checked automatically.\n",
        "Compare policies within the same variant and scenario. Pilot-size changes also change total measured candidates and the pilot-derived threshold. Noise changes affect noisy labels and their pilot-derived threshold. Consequently, unadjusted hit counts across variants do not measure one common fixed-threshold estimand. Confirmed hits include the pilot; active-only and latent-hit counts are saved separately.\n",
        f"For final measured hits, CALDER−PoH mean differences were negative in {n_negative}, zero in {n_zero} and positive in {n_positive} of the 36 scenario/variant cells. These descriptive cell counts are not independent statistical replications. No cross-scenario confidence interval is computed.\n",
        "## Means by scenario and variant\n",
        markdown(
            summary[
                [
                    "variant",
                    "scenario",
                    "policy",
                    "confirmed_hits",
                    "active_hits",
                    "group_coverage",
                    "brier_score",
                    "calibration_fitted_fraction",
                    "n_selected",
                ]
            ]
        ),
        "\n## Paired measured-hit differences\n",
        "Positive values favour CALDER. Intervals bootstrap whole paired seeds separately in each scenario and variant; three seeds give a limited exploratory interval, conditional on this fixed synthetic pool. Other metric contrasts are in paired_differences.csv. Multiple contrasts are descriptive, with no adjusted superiority claim.\n",
        markdown(
            primary_differences[
                ["variant", "scenario", "mean_difference", "ci_lower", "ci_upper", "n_pairs"]
            ]
        ),
        "\n## Warnings and limitations\n",
        markdown(warnings_frame),
        "\nSmall calibration reserves fall back to raw probabilities and retain explicit warnings. Calibration is fitted against measured pilot labels, whereas raw epistemic probabilities concern latent effects; prequential acquired-label diagnostics remain selection biased. Repeated seeds are algorithmic runs, not biological replications.\n",
        "\n## Reproduction\n",
        "From a frozen installation at repository root, run `uv run --frozen --extra dev python scripts/run_sensitivity.py --config configs/sensitivity.yaml --output-dir results/reference_sensitivity`. Existing completed variants are reused only when their exact config hashes match. Use another output directory for fresh replay. Each variant retains config.json, source hashes, Git metadata, full seed endpoints, selections, predictions and audit logs. make_report=False avoids redundant per-variant figure generation; use `calder report --run-dir <variant-directory>` when needed. No absent figure is replaced by invented values.\n",
    ]
    (root / "summary.md").write_text("\n".join(text))
    manifest = {
        "schema_version": 1,
        "status": "complete",
        "experiment_type": "exploratory-one-factor-sensitivity",
        "base_config": str(args.config.resolve()),
        "base_config_sha256": file_sha256(args.config),
        "grid": {axis: list(values) for axis, values in GRID},
        "source_runs": source_runs,
        "n_trajectories": len(endpoints),
        "checks": {
            "paired_thresholds": True,
            "paired_realized_budgets": True,
            "three_base_setting_trajectories_identical": True,
        },
        "files": {
            name: file_sha256(root / name)
            for name in [
                "summary.csv",
                "seed_level_results.csv",
                "paired_differences.csv",
                "warning_summary.csv",
                "summary.md",
            ]
        },
    }
    (root / "experiment_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"Saved {len(endpoints)} trajectories; CALDER−PoH hit differences: {n_negative} negative, {n_zero} zero, {n_positive} positive",
        flush=True,
    )


if __name__ == "__main__":
    main()
