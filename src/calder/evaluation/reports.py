"""Audit-friendly summaries derived from immutable saved trajectory outputs."""

from __future__ import annotations

import datetime as dt
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from calder.evaluation.calibration import (
    predictive_interval_coverage,
    probability_calibration,
    reliability_bins,
    uncertainty_error_association,
)
from calder.evaluation.plots import STRATA, plot_saved_results
from calder.evaluation.statistics import bootstrap_interval, paired_policy_difference
from calder.utils.serialization import environment_provenance, file_sha256, write_json


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size < 2:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _markdown_table(frame: pd.DataFrame) -> str:
    """Keep reporting independent of the optional pandas/tabulate dependency."""
    if frame.empty:
        return "No matching measured results are available.\n"
    display = frame.copy()
    for column in display:
        if pd.api.types.is_numeric_dtype(display[column]):
            display[column] = display[column].map(lambda x: f"{x:.4g}" if pd.notna(x) else "NA")
    records = [list(map(str, display.columns)), ["---"] * len(display.columns)]
    records += display.astype(str).values.tolist()
    return (
        "\n".join(
            "| " + " | ".join(value.replace("|", "\\|") for value in row) + " |" for row in records
        )
        + "\n"
    )


def summarize_endpoints(endpoints: pd.DataFrame) -> pd.DataFrame:
    """CI across independent random seeds within each context and threshold."""
    if endpoints.empty:
        return pd.DataFrame()
    keys = [k for k in (*STRATA, "policy") if k in endpoints]
    metric_names = [
        k
        for k in (
            "confirmed_hits",
            "latent_hits",
            "group_coverage",
            "precision",
            "recall",
            "auc_hits",
            "runtime_seconds",
            "peak_memory_mb",
        )
        if k in endpoints
    ]
    records = []
    for group, part in endpoints.groupby(keys, dropna=False, sort=True):
        group = group if isinstance(group, tuple) else (group,)
        record = dict(zip(keys, group, strict=True))
        record["n_seeds"] = int(part.seed.nunique()) if "seed" in part else len(part)
        for metric in metric_names:
            values = pd.to_numeric(part[metric], errors="coerce").dropna().to_numpy()
            if len(values):
                mean, lo, hi = bootstrap_interval(values)
                record.update(
                    {f"{metric}_mean": mean, f"{metric}_ci_lower": lo, f"{metric}_ci_upper": hi}
                )
        records.append(record)
    return pd.DataFrame(records)


def calibration_summary(prequential: pd.DataFrame) -> pd.DataFrame:
    """Prequential diagnostics, reported separately for latent/measured targets.

    Active examples are selected adaptively. These descriptive scores do not
    establish exchangeability, population calibration, or biological replication.
    """
    if prequential.empty or not {"probability", "is_hit"}.issubset(prequential.columns):
        return pd.DataFrame()
    keys = [k for k in (*STRATA, "policy", "seed") if k in prequential]
    records = []
    for group, part in prequential.groupby(keys, dropna=False, sort=True):
        group = group if isinstance(group, tuple) else (group,)
        identity = dict(zip(keys, group, strict=True))
        targets = [("measured", "is_hit")]
        if "latent_label" in part and part.latent_label.notna().all():
            targets.append(("latent", "latent_label"))
        for target_name, target_column in targets:
            for prediction_column in ("raw_probability", "calibrated_probability", "probability"):
                if prediction_column not in part:
                    continue
                valid = part[[prediction_column, target_column]].notna().all(axis=1)
                filtered = part[valid]
                scores = probability_calibration(
                    filtered[prediction_column].to_numpy(float),
                    filtered[target_column].to_numpy(int),
                )
                record = {
                    **identity,
                    "probability_type": prediction_column,
                    "label_target": target_name,
                    **scores,
                }
                if {"mean", "observed_effect", "epistemic_variance"}.issubset(filtered.columns):
                    latent_column = next(
                        (k for k in ("latent_effect", "true_log_effect") if k in filtered), None
                    )
                    response_column = (
                        latent_column
                        if target_name == "latent" and latent_column
                        else "observed_effect"
                    )
                    errors = np.abs(
                        filtered["mean"].to_numpy(float) - filtered[response_column].to_numpy(float)
                    )
                    record["mae"] = float(errors.mean())
                    record["uncertainty_error_spearman"] = uncertainty_error_association(
                        filtered.epistemic_variance, errors
                    )
                    variance = filtered.epistemic_variance.to_numpy(float)
                    if target_name == "measured" and "noise_variance" in filtered:
                        variance = variance + filtered.noise_variance.to_numpy(float)
                        if "model_noise_variance" in filtered:
                            variance = variance + filtered.model_noise_variance.to_numpy(float)
                            record["interval_target"] = "measured_predictive_model_assumed"
                        else:
                            record["interval_target"] = "measured_predictive_sampling_only"
                    elif target_name == "latent":
                        record["interval_target"] = "latent_epistemic"
                    else:
                        record["interval_target"] = (
                            "measured_without_observation_noise_not_comparable"
                        )
                    record["interval_coverage_95"] = predictive_interval_coverage(
                        filtered["mean"], variance, filtered[response_column]
                    )
                records.append(record)
    return pd.DataFrame(records)


def _model_noise_from_config(prequential: pd.DataFrame, config_path: Path) -> pd.DataFrame:
    """Restore a likelihood-discrepancy floor for older saved runs, in memory only."""
    if prequential.empty or not config_path.exists():
        return prequential
    config = json.loads(config_path.read_text())
    frame = prequential.copy()
    if "model_noise_variance" not in frame:
        frame["model_noise_variance"] = np.nan
    if "model" in config and "noise_floor" in config["model"]:
        frame["model_noise_variance"] = frame.model_noise_variance.fillna(
            float(config["model"]["noise_floor"])
        )
    elif "source_configs" in config and "dataset" in frame:
        for source in config["source_configs"]:
            name = source.get("dataset", {}).get("name")
            floor = source.get("model", {}).get("noise_floor")
            if name is not None and floor is not None:
                mask = (frame.dataset == name) & frame.model_noise_variance.isna()
                frame.loc[mask, "model_noise_variance"] = float(floor)
    if frame.model_noise_variance.isna().all():
        return frame.drop(columns="model_noise_variance")
    if not np.isfinite(frame.model_noise_variance).all() or (frame.model_noise_variance < 0).any():
        raise ValueError("Configured model-discrepancy variances must be finite and nonnegative")
    return frame


def detailed_calibration_tables(prequential: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Save bin estimates and observed-group diagnostics without reading hidden labels."""
    if prequential.empty or not {"probability", "is_hit"}.issubset(prequential.columns):
        return pd.DataFrame(), pd.DataFrame()
    keys = [k for k in (*STRATA, "policy", "seed") if k in prequential]
    bin_tables, group_records = [], []
    for group, part in prequential.groupby(keys, dropna=False, sort=True):
        values = group if isinstance(group, tuple) else (group,)
        identity = dict(zip(keys, values, strict=True))
        targets = [("measured", "is_hit")]
        if "latent_label" in part and part.latent_label.notna().all():
            targets.append(("latent", "latent_label"))
        groups = None
        if "groups" in part:
            try:
                groups = part["groups"].map(
                    lambda value: json.loads(value) if isinstance(value, str) else []
                )
            except (TypeError, ValueError):
                # Malformed group annotations must not be silently reinterpreted.
                raise ValueError(
                    "Cannot parse saved acquired-candidate group annotations"
                ) from None
        for target_name, target_column in targets:
            for prediction_column in ("raw_probability", "calibrated_probability", "probability"):
                if prediction_column not in part:
                    continue
                stamp = {
                    **identity,
                    "probability_type": prediction_column,
                    "label_target": target_name,
                }
                table = reliability_bins(
                    part[prediction_column].to_numpy(float), part[target_column].to_numpy(int)
                )
                for key, value in stamp.items():
                    table[key] = value
                bin_tables.append(table)
                if groups is not None:
                    names = sorted({str(name) for membership in groups for name in membership})
                    for name in names:
                        mask = groups.map(
                            lambda membership: name in {str(value) for value in membership}
                        )
                        metrics = probability_calibration(
                            part.loc[mask, prediction_column],
                            part.loc[mask, target_column].astype(int),
                        )
                        group_records.append({**stamp, "group": name, **metrics})
    return (
        pd.concat(bin_tables, ignore_index=True) if bin_tables else pd.DataFrame(),
        pd.DataFrame(group_records),
    )


def build_report(run_dir: Path | str) -> dict[str, object]:
    """Regenerate tables, figures and narrative from genuine saved run outputs."""
    root = Path(run_dir)
    if not root.is_dir():
        raise ValueError(f"Run directory does not exist: {root}")
    rounds, endpoints = _read_csv(root / "rounds.csv"), _read_csv(root / "seed_level_results.csv")
    if endpoints.empty:
        if rounds.empty:
            raise ValueError(f"No saved experiment outputs found in {root}; run a benchmark first")
        keys = [k for k in (*STRATA, "policy", "seed") if k in rounds]
        endpoints = rounds.sort_values("round").groupby(keys, dropna=False, sort=True).tail(1)
        endpoints.to_csv(root / "seed_level_results.csv", index=False)
    strata = [k for k in STRATA if k in endpoints]
    if (
        "seed" in endpoints
        and "policy" in endpoints
        and endpoints.duplicated([*strata, "policy", "seed"]).any()
    ):
        raise ValueError(
            "Duplicate seed-level endpoints: the run identifiers do not distinguish experimental conditions"
        )
    summary = summarize_endpoints(endpoints)
    summary.to_csv(root / "results.csv", index=False)
    prequential = _model_noise_from_config(
        _read_csv(root / "prequential.csv"), root / "config.json"
    )
    calibration = calibration_summary(prequential)
    calibration.to_csv(root / "calibration.csv", index=False)
    reliability, grouped = detailed_calibration_tables(prequential)
    reliability.to_csv(root / "reliability.csv", index=False)
    grouped.to_csv(root / "grouped_calibration.csv", index=False)
    comparison_records, warnings = [], []
    if "policy" in endpoints:
        for group, part in (
            endpoints.groupby(strata, dropna=False, sort=True) if strata else [((), endpoints)]
        ):
            values = group if isinstance(group, tuple) else (group,)
            identity = dict(zip(strata, values, strict=True))
            policies = sorted(part.policy.unique())
            reference = (
                "probability_of_hit"
                if "probability_of_hit" in policies
                else "random"
                if "random" in policies
                else policies[0]
            )
            for policy in policies:
                if policy == reference or policy.startswith("oracle"):
                    continue
                for metric in ("confirmed_hits", "group_coverage", "auc_hits"):
                    if metric not in part:
                        continue
                    try:
                        difference = paired_policy_difference(part, policy, reference, metric)
                    except ValueError as error:
                        warnings.append(f"{identity}: {policy} vs {reference} {metric}: {error}")
                        continue
                    comparison_records.append(
                        {
                            **identity,
                            "policy": policy,
                            "reference": reference,
                            "metric": metric,
                            **difference,
                        }
                    )
    comparisons = pd.DataFrame(comparison_records)
    comparisons.to_csv(root / "paired_differences.csv", index=False)
    figure_status = plot_saved_results(root)
    manifest_path = root / "experiment_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    fields = [
        k
        for k in (
            *STRATA,
            "policy",
            "n_seeds",
            "confirmed_hits_mean",
            "group_coverage_mean",
            "auc_hits_mean",
        )
        if k in summary
    ]
    biological = (
        "dataset" in endpoints
        and endpoints.dataset.astype(str).str.contains("replogle22", regex=False).any()
    )
    text = [
        "# CALDER experiment summary\n",
        "This report was generated from saved experiment outputs. No real laboratory experiments were performed by CALDER.\n",
        f"Experiment type: {'retrospective replay of public perturbation measurements' if biological else 'controlled synthetic simulation'}.\n",
        f"Configuration hash: `{manifest.get('config_hash', 'not recorded')}`. Git commit: `{manifest.get('git_commit', 'not recorded')}`.\n",
        "## Measured endpoint results\n",
        "Confirmed hits use observed log effects at the pilot-derived threshold. Synthetic latent-hit counts, when present, are a distinct ground-truth diagnostic. Group coverage counts predefined groups containing at least one confirmed measured hit.\n",
        _markdown_table(summary[fields]),
        "Full bootstrap intervals and secondary metrics are in `results.csv`; every run endpoint is in `seed_level_results.csv`. Seeds describe algorithmic variation on a fixed experimental pool, not independent biological replications.\n",
        "## Paired policy differences\n",
        "Positive differences favour the named policy. Intervals resample whole paired seeds within each dataset, scenario and threshold. They quantify Monte Carlo variation conditional on this pool. Multiple contrasts are descriptive; no multiplicity-adjusted superiority claim is made.\n",
        _markdown_table(comparisons),
        "## Calibration and prediction\n",
        "`calibration.csv` reports raw, pilot-calibrated and acquisition-used probability diagnostics: Brier score, fixed-bin ECE, prediction error, uncertainty–error association and interval coverage. Bin counts and group-specific diagnostics are in `reliability.csv` and `grouped_calibration.csv`. These are pre-query predictions on adaptively acquired examples. Selection bias prevents interpreting them as exchangeable population-level validation. Latent probabilities evaluated against noisy measured hit labels also reflect label noise. Calibration fitted to measured pilot hit labels is a measured-label correction, not guaranteed latent-effect calibration.\n",
        "Measured predictive intervals add the recorded sampling variance and configured likelihood-discrepancy variance floor to epistemic variance. Their target is labelled `measured_predictive_model_assumed`; this is a model assumption and is heuristic for ridge or ensemble uncertainty. Older saved runs restore the floor from their resolved configuration in memory, preserving original prediction tables. Latent-effect intervals use epistemic variance alone.\n",
        "## Figure availability\n",
    ]
    for figure, item in figure_status.items():
        text.append(
            f"- {figure}: {item['status']}. "
            + (
                ", ".join(f"[view](figures/{file})" for file in item["files"])
                if item["files"]
                else str(item.get("reason", ""))
            )
            + "\n"
        )
    if warnings:
        text += ["\n## Analysis limitations\n", *[f"- {warning}\n" for warning in warnings]]
    text += [
        "\n## Interpretation limits\n",
        "The policy is not guaranteed to outperform Probability-of-Hit or other baselines. Negative and mixed outcomes remain in the saved tables. A high transcriptional-response score does not establish therapeutic usefulness, and perturbation-induced expression changes do not identify direct causal regulatory edges. STRING communities are prior association groups; control-feature clusters are artificial feature clusters. Projected feature geometry has no direct causal interpretation.\n",
        "K562–RPE1 comparison requires measured outputs from both contexts under frozen method hyperparameters; missing context results are not replaced with synthetic values. Offline oracle rankings are full-information upper-bound references and are not deployable policies.\n",
        "\nSources and protocol are documented in the repository's `reports/bibliography.md` and `reports/research_report.md`.\n",
    ]
    (root / "summary.md").write_text("\n".join(text))
    _write_research_report(root, endpoints, summary, comparisons, biological)
    _update_report_manifest(root, manifest)
    return {
        "results": str(root / "results.csv"),
        "summary": str(root / "summary.md"),
        "figure_status": figure_status,
        "warnings": warnings,
    }


def _update_report_manifest(root: Path, manifest: dict) -> None:
    """Refresh derived-file hashes while preserving original experimental provenance."""
    names = (
        "results.csv",
        "calibration.csv",
        "reliability.csv",
        "grouped_calibration.csv",
        "paired_differences.csv",
        "summary.md",
        "research_report.md",
        "bibliography.md",
        "figure_status.json",
    )
    files = manifest.setdefault("files", {})
    for name in names:
        path = root / name
        if path.is_file():
            files[name] = file_sha256(path)
    for path in sorted((root / "figures").glob("*.png")):
        files[str(path.relative_to(root))] = file_sha256(path)
    provenance = environment_provenance()
    template = Path(__file__).resolve().parents[3] / "reports" / "protocol_template.md"
    bibliography = template.with_name("bibliography.md")
    manifest["report_generation_provenance"] = {
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_tree_hash": provenance.get("source_tree_hash"),
        "source_files": provenance.get("source_files", {}),
        "git_commit": provenance.get("git_commit"),
        "git_dirty": provenance.get("git_dirty"),
        "git_tracked_dirty": provenance.get("git_tracked_dirty"),
        "report_template_sha256": file_sha256(template) if template.is_file() else None,
        "bibliography_template_sha256": file_sha256(bibliography)
        if bibliography.is_file()
        else None,
    }
    write_json(root / "experiment_manifest.json", manifest)


def _write_research_report(
    root: Path,
    endpoints: pd.DataFrame,
    summary: pd.DataFrame,
    comparisons: pd.DataFrame,
    biological: bool,
) -> None:
    """Write an empirical abstract only after endpoint measurements exist."""
    n_policies = endpoints.policy.nunique() if "policy" in endpoints else 0
    n_seeds = endpoints.seed.nunique() if "seed" in endpoints else 0
    abstract = (
        f"We evaluated budgeted perturbation selection using {len(endpoints)} saved "
        f"trajectories, {n_policies} selection/reference rules and {n_seeds} random seeds "
        f"in {'retrospective biological replay' if biological else 'controlled synthetic simulation'}. "
    )
    if not comparisons.empty:
        calder = comparisons[
            (comparisons.policy == "calder") & (comparisons.reference == "probability_of_hit")
        ]
        clauses = []
        for metric, label in [
            ("confirmed_hits", "confirmed-hit"),
            ("group_coverage", "group-coverage"),
        ]:
            values = calder.loc[calder.metric == metric, "mean_difference"]
            if not values.empty:
                clauses.append(
                    f"{label} mean differences ranged from {values.min():.3g} to {values.max():.3g}"
                )
        if clauses:
            abstract += (
                "Across separately analyzed dataset/scenario/threshold strata, CALDER minus Probability-of-Hit "
                + "; ".join(clauses)
                + ". These ranges describe separate effects rather than a pooled estimate. "
            )
    abstract += (
        "Seed-paired bootstrap intervals quantify algorithmic variation conditional on the fixed pools. "
        "No laboratory experiments were performed, and these results do not establish universal policy superiority."
    )
    if not biological:
        abstract += " Biological K562–RPE1 performance remains unmeasured in this run."
    protocol_path = Path(__file__).resolve().parents[3] / "reports" / "protocol_template.md"
    bibliography_path = protocol_path.with_name("bibliography.md")
    if protocol_path.exists():
        document = protocol_path.read_text()
        document = re.sub(
            r"## Abstract\n.*?(?=\n## 1\.)",
            "## Abstract\n\n" + abstract + "\n",
            document,
            flags=re.DOTALL,
        )
        fields = [
            k
            for k in (
                *STRATA,
                "policy",
                "n_seeds",
                "confirmed_hits_mean",
                "group_coverage_mean",
                "auc_hits_mean",
            )
            if k in summary
        ]
        results_section = (
            f"## 6. Results\n\nRun: `{root.name}`. Input configuration and full provenance are saved in `config.json` and `experiment_manifest.json`.\n\n"
            + _markdown_table(summary[fields])
            + "\nPaired differences (policy minus reference), conditional on each fixed pool:\n\n"
            + _markdown_table(comparisons)
            + "\nFull confidence intervals are in `results.csv`; seed endpoints are in `seed_level_results.csv`. Prequential latent/measured diagnostics remain separate in `calibration.csv`. Figure availability and links are recorded in [summary.md](summary.md).\n"
        )
        document = re.sub(
            r"## 6\. Results\n.*?(?=\n## 7\.)", lambda _: results_section, document, flags=re.DOTALL
        )
        ablations = _read_csv(root / "ablations.csv")
        if ablations.empty:
            ablation_text = "## 7. Ablation studies\n\nThis run contains no completed ablation trajectories. No ablation effect is claimed. The predeclared ablation configuration must be executed separately.\n"
        else:
            ablation_fields = [
                k
                for k in (
                    *STRATA,
                    "policy",
                    "n_seeds",
                    "confirmed_hits_mean",
                    "group_coverage_mean",
                )
                if k in summary
            ]
            ablation_text = (
                "## 7. Ablation studies\n\nCompleted controlled variants are preserved in `ablations.csv`; component differences must be assessed within each fixed scenario and threshold.\n\n"
                + _markdown_table(summary[ablation_fields])
            )
        document = re.sub(
            r"## 7\. Ablation [Ss]tudies\n.*?(?=\n## 8\.)",
            lambda _: ablation_text,
            document,
            flags=re.DOTALL,
        )
        document = re.sub(
            r"## 10\. Conclusion\n.*?(?=\n## 11\.)",
            lambda _: "## 10. Conclusion\n\n" + abstract + "\n",
            document,
            flags=re.DOTALL,
        )
        if biological:
            document = re.sub(
                r"At implementation-time verification on 3 November 2026,.*?(?=\n\n)",
                "The original implementation encountered remote download failures. This run's actual processed biological inputs and checksums are recorded in its experiment manifest; historical source failures do not substitute for those measured inputs.",
                document,
                flags=re.DOTALL,
            )
    else:
        document = (
            "# CALDER empirical research report\n\n## Abstract\n\n"
            + abstract
            + "\n\n"
            + (root / "summary.md").read_text()
        )
    (root / "research_report.md").write_text(document)
    if bibliography_path.exists():
        (root / "bibliography.md").write_text(bibliography_path.read_text())


def build_context_comparison(
    k562_run: Path | str, rpe1_run: Path | str, output_dir: Path | str
) -> dict[str, object]:
    """Compare completed biological reruns with equal method hyperparameters.

    This is context re-evaluation (new pilot and sequential observations), not
    transfer of a K562-trained predictive model. No hidden labels are retuned.
    """
    sources = [Path(k562_run), Path(rpe1_run)]
    configs = [json.loads((source / "config.json").read_text()) for source in sources]
    for key in ("model", "acquisition"):
        if configs[0].get(key) != configs[1].get(key):
            raise ValueError(f"Context comparison requires frozen {key} hyperparameters")
    protocol_keys = (
        "pilot_fraction",
        "calibration_fraction",
        "calibration_min_samples",
        "batch_size",
        "budget_fractions",
        "threshold_quantiles",
        "policies",
    )
    for key in protocol_keys:
        if configs[0].get("protocol", {}).get(key) != configs[1].get("protocol", {}).get(key):
            raise ValueError(f"Context comparison requires the same preregistered protocol: {key}")
    endpoints = [_read_csv(source / "seed_level_results.csv") for source in sources]
    expected = ("replogle22k562", "replogle22rpe1")
    manifests = []
    for source, frame, dataset in zip(sources, endpoints, expected, strict=True):
        if frame.empty or "dataset" not in frame or set(frame.dataset.astype(str)) != {dataset}:
            raise ValueError(f"Expected measured {dataset} endpoints in {source}")
        manifest = json.loads((source / "experiment_manifest.json").read_text())
        if manifest.get("status") != "complete":
            raise ValueError(f"Context source run did not complete: {source}")
        manifests.append(manifest)
    root = Path(output_dir)
    if root.exists() and any(root.iterdir()):
        raise ValueError(f"Refusing to overwrite comparison directory: {root}")
    root.mkdir(parents=True, exist_ok=True)
    for filename in (
        "rounds.csv",
        "seed_level_results.csv",
        "prequential.csv",
        "selections.csv",
        "candidates.csv",
    ):
        tables = [_read_csv(source / filename) for source in sources]
        nonempty = [frame for frame in tables if not frame.empty]
        if nonempty:
            pd.concat(nonempty, ignore_index=True).to_csv(root / filename, index=False)
    comparison_manifest = {
        "schema_version": 1,
        "status": "complete",
        "experiment_type": "independent-context-reruns-with-frozen-methods",
        "direct_predictive_transfer": False,
        "source_runs": [
            {
                "path": str(source.resolve()),
                "config_hash": manifest.get("config_hash"),
                "git_commit": manifest.get("git_commit"),
            }
            for source, manifest in zip(sources, manifests, strict=True)
        ],
        "limitation": "Repeated seeds measure algorithmic variation; K562/RPE1 are distinct fixed contexts, not independent biological replicates.",
    }
    (root / "experiment_manifest.json").write_text(json.dumps(comparison_manifest, indent=2) + "\n")
    (root / "config.json").write_text(json.dumps({"source_configs": configs}, indent=2) + "\n")
    return build_report(root)
