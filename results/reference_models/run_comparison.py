"""Reproduce the declared three-model development comparison from repository root."""

from __future__ import annotations

import json
import warnings
from collections import Counter
from pathlib import Path

import pandas as pd

from calder.config import load_config
from calder.experiments.runner import run_benchmark
from calder.utils.serialization import file_sha256

base = load_config("configs/model_comparison.yaml")
root = Path("results/reference_models")
models = ["bayesian_ridge", "random_features", "ensemble"]
for model in models:
    target = root / model
    if not (target / "experiment_manifest.json").exists():
        config = base.model_copy(deep=True)
        config.model.name = model
        config.experiment_name = f"development_model_comparison_{model}"
        print(f"Running {model}", flush=True)
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            try:
                run_benchmark(config, target)
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
                (target / "warnings.txt").write_text(
                    "\n".join(
                        f"{r['category']}: {r['message']} ({r['filename']}:{r['lineno']})"
                        for r in records
                    )
                    + "\n"
                )
        print(f"Finished {model}: {len(captured)} captured warnings", flush=True)
    else:
        manifest = json.loads((target / "experiment_manifest.json").read_text())
        if manifest["status"] != "complete":
            raise RuntimeError(f"Existing incomplete run: {target}")

frames = [pd.read_csv(root / model / "seed_level_results.csv") for model in models]
all_seeds = pd.concat(frames, ignore_index=True)
assert len(all_seeds) == 108
assert all_seeds.groupby(["scenario", "seed", "threshold_quantile"]).threshold.nunique().max() == 1
random_selections = {}
comparison_keys = ["scenario", "seed", "threshold_quantile", "policy", "selection_order", "gene_id"]
for model in models:
    selected = pd.read_csv(root / model / "selections.csv")
    random_selections[model] = (
        selected[selected.policy == "random"][comparison_keys]
        .sort_values(comparison_keys)
        .reset_index(drop=True)
    )
for model in ["bayesian_ridge", "ensemble"]:
    pd.testing.assert_frame_equal(random_selections["random_features"], random_selections[model])
all_seeds.to_csv(root / "seed_level_results.csv", index=False)
metrics = [
    "confirmed_hits",
    "latent_hits",
    "group_coverage",
    "precision",
    "auc_hits",
    "brier_score",
    "mae",
    "runtime_seconds",
]
summary = (
    all_seeds.groupby(["scenario", "model", "policy"], sort=True)[metrics].mean().reset_index()
)
summary.to_csv(root / "results.csv", index=False)
contrasts = summary.pivot(
    index=["scenario", "model"], columns="policy", values=["confirmed_hits", "group_coverage"]
)
contrasts = pd.DataFrame(
    {
        "hit_difference_calder_minus_poh": contrasts["confirmed_hits"]["calder"]
        - contrasts["confirmed_hits"]["probability_of_hit"],
        "group_difference_calder_minus_poh": contrasts["group_coverage"]["calder"]
        - contrasts["group_coverage"]["probability_of_hit"],
    }
).reset_index()
contrasts.to_csv(root / "policy_differences.csv", index=False)
warning_records = []
for model in models:
    records = json.loads((root / model / "warnings.json").read_text())
    for (category, message), count in Counter(
        (r["category"], r["message"]) for r in records
    ).items():
        warning_records.append(
            {"model": model, "category": category, "message": message, "count": count}
        )
pd.DataFrame(warning_records, columns=["model", "category", "message", "count"]).to_csv(
    root / "warning_summary.csv", index=False
)


def markdown(frame):
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


text = [
    "# Exploratory development model comparison\n",
    "These are actual controlled synthetic runs, not biological results. Three predeclared models were evaluated without hyperparameter tuning on 250 candidates, 24 features, four scenarios, seeds 0–2, threshold quantile 0.8, and policies random, Probability-of-Hit and CALDER. The method hyperparameters of the main benchmark remain unchanged. Three seeds are a limited development check, not final confirmatory evidence.\n",
    "The random-features model uses 128 Fourier features; Bayesian ridge uses the original features. The ensemble uses five bootstrap scikit-learn MLPs with one hidden layer of 32 units and max_iter=300; disagreement is heuristic epistemic uncertainty and does not incorporate per-observation noise weights. This is a five-member MLP ensemble, not a claim to reproduce an external deep-ensemble protocol. All methods receive the same permitted candidate information and pilot split.\n",
    "Means below summarize three seeds separately within each scenario/model/policy. Confirmed hits use noisy measured log effects; latent hits are a separate ground-truth diagnostic. Group coverage uses the original fixed groups. AUC is the left-step hit × cost integral beyond the pilot, including pilot hits as the initial discovery count. Every trajectory spent 50 additional cost units and selected 100 candidates including the common pilot.\n",
    markdown(
        summary[
            [
                "scenario",
                "model",
                "policy",
                "confirmed_hits",
                "latent_hits",
                "group_coverage",
                "auc_hits",
                "brier_score",
                "mae",
            ]
        ]
    ),
    "\n## Observed development comparisons\n",
    "Probability-of-Hit had more mean confirmed hits than CALDER in all 12 scenario/model combinations. CALDER sometimes traded fewer hits for greater group coverage, while several scenarios had coverage already saturated at five groups. In scenario A, PoH with Bayesian ridge found 48.0 mean confirmed hits, the ensemble 46.0 and RFF 37.67; the model-family ranking therefore differed from simply assuming that the nonlinear surrogate must be best. These small comparisons are exploratory and do not change the predeclared main model or its settings.\n",
    markdown(contrasts),
    "\n## Warning record\n",
    markdown(pd.DataFrame(warning_records, columns=["model", "category", "message", "count"])),
    "\nConvergence warnings are retained rather than suppressed or resolved by tuning these runs. A warning means the MLP optimizer exhausted its iteration limit; ensemble comparisons therefore reflect this fixed training budget. The full warning sequence is retained in each model directory as warnings.json and warnings.txt. No claim of ensemble convergence or model-family superiority follows from these small runs.\n",
    "\n## Reproduction and audit\n",
    "Run `.venv/bin/python results/reference_models/run_comparison.py` from the repository root after the frozen installation. The script refuses existing incomplete runs and reuses completed model runs. To create independent replicas, use the saved per-model config.json snapshots with `calder benchmark --config <snapshot> --run-dir <new-empty-directory>`. Each model has its own configuration hash, Git metadata, environment information, all seed-level outcomes, pre-query prediction snapshots, audit logs and figures. Aggregate files here concatenate those real outputs.\n",
    "Repeated seeds describe algorithmic variation on one fixed synthetic pool; they are not biological replications. Calibration diagnostics on actively selected examples are selection biased. Negative/mixed results are retained; results must not be used to retune the final RPE1 benchmark.\n",
]
(root / "summary.md").write_text("\n".join(text))
sources = []
for model in models:
    source = json.loads((root / model / "experiment_manifest.json").read_text())
    sources.append(
        {
            "model": model,
            "path": model,
            "status": source["status"],
            "config_hash": source["config_hash"],
            "git_commit": source.get("git_commit"),
            "git_dirty": source.get("git_dirty"),
            "n_trajectories": source["n_trajectories"],
        }
    )
(root / "experiment_manifest.json").write_text(
    json.dumps(
        {
            "schema_version": 1,
            "status": "complete",
            "experiment_type": "exploratory-development-model-comparison",
            "source_runs": sources,
            "checks": {
                "trajectories": 108,
                "identical_random_trajectories_across_models": True,
                "paired_thresholds_across_models_and_policies": True,
            },
            "files": {
                name: file_sha256(root / name)
                for name in [
                    "seed_level_results.csv",
                    "results.csv",
                    "policy_differences.csv",
                    "warning_summary.csv",
                    "summary.md",
                ]
            },
        },
        indent=2,
    )
    + "\n"
)
print(
    summary[["scenario", "model", "policy", "confirmed_hits", "group_coverage"]].to_string(
        index=False
    ),
    flush=True,
)
print("Saved comparison summary and all seed-level rows", flush=True)
