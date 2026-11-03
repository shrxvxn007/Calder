"""Paired experiments, evaluation-only access to hidden labels, and provenance."""

from __future__ import annotations

import datetime as dt
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from calder.config import ExperimentConfig, load_config
from calder.experiments.active_loop import Trajectory, run_trajectory
from calder.experiments.oracle import OfflineOracle
from calder.experiments.simulator import simulate
from calder.experiments.splits import pilot_split
from calder.types import CandidatePool, OutcomeTable
from calder.utils.serialization import content_hash, environment_provenance, file_sha256, write_json


def evaluate_trajectory(
    trajectory: Trajectory, outcomes: OutcomeTable, pool: CandidatePool
) -> None:
    """Full labels are used only AFTER an entire trajectory has been selected."""
    total_hits = int(np.sum(outcomes.measured_effects > trajectory.threshold))
    for row in trajectory.rounds:
        selected = [int(r["index"]) for r in trajectory.selections if r["round"] <= row["round"]]
        row["recall"] = row["confirmed_hits"] / total_hits if total_hits else 0.0
        row["total_pool_hits"] = total_hits
        if outcomes.latent_effects is not None:
            labels = outcomes.latent_effects[selected] > trajectory.threshold
            row["latent_hits"] = int(labels.sum())
            row["latent_group_coverage"] = int(
                np.any(pool.memberships[np.array(selected)[labels]] > 0, axis=0).sum()
            )
    for record in trajectory.prequential:
        if outcomes.latent_effects is not None:
            value = float(outcomes.latent_effects[record["index"]])
            record["latent_label"] = bool(value > trajectory.threshold)
            record["latent_effect"] = value
            record["latent_absolute_error"] = abs(value - record["mean"])


def oracle_reference(
    pool: CandidatePool,
    outcomes: OutcomeTable,
    config: ExperimentConfig,
    seed: int,
    quantile: float,
) -> Trajectory:
    """Observed-label full-information ranking; equal-cost hit-count bound only.

    This is evaluation code and never calls or impersonates an acquisition policy.
    Coverage is incidental; this is not a coverage upper bound.
    """
    train, cal = pilot_split(
        len(pool.gene_ids),
        config.protocol.pilot_fraction,
        config.protocol.calibration_fraction,
        seed,
    )
    pilot = np.concatenate([train, cal])
    threshold = float(np.quantile(outcomes.measured_effects[train], quantile))
    available = np.setdiff1d(np.arange(len(pool.gene_ids)), pilot)
    ranking = available[np.lexsort((available, -outcomes.measured_effects[available]))]
    selected = pilot.tolist()
    selections = [
        {
            "round": 0,
            "index": int(i),
            "gene_id": str(pool.gene_ids[i]),
            "split": "pilot_train" if i in train else "pilot_calibration",
        }
        for i in pilot
    ]
    rows = []
    spent = 0.0
    pilot_cost = float(pool.costs[pilot].sum())

    def append(round_id: int, fraction: float) -> None:
        hits = outcomes.measured_effects[selected] > threshold
        rows.append(
            {
                "round": round_id,
                "threshold": threshold,
                "budget_spent": pilot_cost + spent,
                "additional_budget": spent,
                "budget_fraction": spent / len(pool.gene_ids),
                "checkpoint_fraction": fraction,
                "confirmed_hits": int(hits.sum()),
                "group_coverage": int(
                    np.any(pool.memberships[np.array(selected)[hits]] > 0, axis=0).sum()
                ),
                "precision": float(hits.mean()),
                "n_selected": len(selected),
                "active_hits": int(
                    (outcomes.measured_effects[selected[len(pilot) :]] > threshold).sum()
                ),
                "runtime_seconds": 0.0,
            }
        )

    append(0, 0.0)
    r = 0
    for fraction in config.protocol.budget_fractions:
        while spent < fraction * len(pool.gene_ids) - 1e-9:
            batch = []
            batch_cost = 0.0
            for i in ranking:
                if (
                    int(i) not in selected
                    and spent + batch_cost + pool.costs[i] <= fraction * len(pool.gene_ids) + 1e-9
                ):
                    batch.append(int(i))
                    batch_cost += float(pool.costs[i])
                    if len(batch) == config.protocol.batch_size:
                        break
            if not batch:
                break
            r += 1
            selected.extend(batch)
            spent += batch_cost
            selections.extend(
                [
                    {"round": r, "index": i, "gene_id": str(pool.gene_ids[i]), "split": "reference"}
                    for i in batch
                ]
            )
            append(
                r, fraction if abs(spent - fraction * len(pool.gene_ids)) < 1e-9 else float("nan")
            )
    return Trajectory(rows, selections, [], threshold, train, cal)


def _dataset_pairs(config: ExperimentConfig) -> list[tuple[str, CandidatePool, OutcomeTable]]:
    if config.dataset.name == "synthetic":
        return [
            (scenario, *simulate(config.synthetic, scenario))
            for scenario in config.synthetic.scenarios
        ]
    from calder.data.adapters import load_prepared

    if config.dataset.prepared_path is None:
        raise ValueError(
            "Real benchmark needs dataset.prepared_path. Run 'calder data prepare' first; no synthetic substitution is performed."
        )
    pool, outcomes = load_prepared(config.dataset.prepared_path)
    return [(config.dataset.name, pool, outcomes)]


def run_benchmark(
    config: ExperimentConfig | str | Path,
    run_dir: str | Path | None = None,
    make_report: bool = True,
    ablation_mode: bool = False,
) -> Path:
    if not isinstance(config, ExperimentConfig):
        config = load_config(config)
    payload = config.model_dump(mode="json")
    config_hash = content_hash(payload)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = Path(
        run_dir or Path(config.output_dir) / f"{config.experiment_name}_{stamp}_{config_hash[:8]}"
    )
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Refusing to overwrite existing run: {output}")
    # Validate/load input before creating a success-looking output directory.
    pairs = _dataset_pairs(config)
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "config.json", payload)
    started = time.perf_counter()
    provenance = environment_provenance()
    manifest = {
        "schema_version": 1,
        "status": "running",
        "config_hash": config_hash,
        "created_utc": stamp,
        **provenance,
        "datasets": [],
        "evaluation": "retrospective replay; repeated seeds are not independent biological replications",
        "probability_target": "raw epistemic p targets latent log effect; isotonic mapping targets measured-hit labels",
        "threshold_protocol": "quantile of pilot training responses only; fixed across paired policies",
        "calibration_protocol": "fixed random disjoint pilot reserve; refit mapping per round; never surrogate training",
        "primary_metrics_include_pilot": True,
    }
    write_json(output / "experiment_manifest.json", manifest)
    all_rounds, all_selections, all_prequential, endpoints, candidates = [], [], [], [], []
    try:
        with threadpool_limits(limits=1):
            for scenario, original_pool, outcomes in pairs:
                for i, gene in enumerate(original_pool.gene_ids):
                    candidates.append(
                        {
                            "dataset": config.dataset.name,
                            "scenario": scenario,
                            "index": i,
                            "gene_id": gene,
                            "feature_0": original_pool.features[i, 0],
                            "feature_1": original_pool.features[i, 1]
                            if original_pool.features.shape[1] > 1
                            else 0,
                        }
                    )
                manifest["datasets"].append(
                    {
                        "scenario": scenario,
                        "metadata": original_pool.metadata,
                        "candidate_hash": content_hash(
                            {
                                "ids": original_pool.gene_ids,
                                "features": original_pool.features,
                                "groups": original_pool.memberships,
                            }
                        ),
                        "evaluation_outcome_hash": content_hash(
                            {
                                "values": outcomes.measured_effects,
                                "variances": outcomes.noise_variances,
                                "latent": outcomes.latent_effects,
                            }
                        ),
                    }
                )
                for seed in config.protocol.seeds:
                    for quantile in config.protocol.threshold_quantiles:
                        policies = list(config.protocol.policies)
                        if config.protocol.include_oracle:
                            policies.append("oracle_reference")
                        for policy in policies:
                            pool = original_pool
                            if policy == "calder_random_groups":
                                shuffled = np.random.default_rng(seed + 6021).permutation(
                                    len(pool.gene_ids)
                                )
                                pool = CandidatePool(
                                    pool.gene_ids,
                                    pool.features,
                                    pool.costs,
                                    pool.memberships[shuffled],
                                    pool.group_names,
                                    pool.metadata,
                                )
                            elif policy == "calder_uninformative_features":
                                features = np.random.default_rng(seed + 9011).normal(
                                    size=pool.features.shape
                                )
                                pool = CandidatePool(
                                    pool.gene_ids,
                                    features,
                                    pool.costs,
                                    pool.memberships,
                                    pool.group_names,
                                    pool.metadata,
                                )
                            base = {
                                "run_id": output.name,
                                "dataset": config.dataset.name,
                                "scenario": scenario,
                                "policy": policy,
                                "seed": seed,
                                "threshold_quantile": quantile,
                                "model": config.model.name,
                                "config_hash": config_hash,
                            }
                            trace_id = f"{scenario}_seed{seed}_q{quantile:g}_{policy}"
                            if policy == "oracle_reference":
                                trace = oracle_reference(pool, outcomes, config, seed, quantile)
                            else:
                                trace = run_trajectory(
                                    pool,
                                    OfflineOracle(outcomes),
                                    config,
                                    seed,
                                    quantile,
                                    policy,
                                    output / "trajectories" / trace_id,
                                )
                            evaluate_trajectory(trace, outcomes, original_pool)
                            # Random group ablation decisions use corrupted groups, evaluation uses fixed original groups.
                            for row in trace.rounds:
                                selected = np.array(
                                    [
                                        r["index"]
                                        for r in trace.selections
                                        if r["round"] <= row["round"]
                                    ],
                                    dtype=int,
                                )
                                selected_hits = selected[
                                    outcomes.measured_effects[selected] > trace.threshold
                                ]
                                row["group_coverage"] = int(
                                    np.any(
                                        original_pool.memberships[selected_hits] > 0, axis=0
                                    ).sum()
                                )
                            all_rounds.extend([{**base, **r} for r in trace.rounds])
                            all_selections.extend(
                                [
                                    {**base, "threshold": trace.threshold, **r}
                                    for r in trace.selections
                                ]
                            )
                            all_prequential.extend([{**base, **r} for r in trace.prequential])
                            rows = trace.rounds
                            from calder.evaluation.hit_metrics import area_under_discovery

                            auc = area_under_discovery(
                                [r["additional_budget"] for r in rows],
                                [r["confirmed_hits"] for r in rows],
                            )
                            final = {
                                **base,
                                **rows[-1],
                                "auc_hits": auc,
                                "auc_hits_normalized": auc / rows[-1]["additional_budget"]
                                if rows[-1]["additional_budget"]
                                else 0.0,
                                "threshold": trace.threshold,
                                "pilot_train_count": len(trace.pilot_train),
                                "pilot_calibration_count": len(trace.pilot_calibration),
                            }
                            if trace.prequential:
                                final["brier_score"] = float(
                                    np.mean(
                                        [
                                            (r["probability"] - r["is_hit"]) ** 2
                                            for r in trace.prequential
                                        ]
                                    )
                                )
                                final["mae"] = float(
                                    np.mean([r["absolute_error"] for r in trace.prequential])
                                )
                                final["calibration_fitted_fraction"] = float(
                                    np.mean([r["calibration_fitted"] for r in trace.prequential])
                                )
                            endpoints.append(final)
        pd.DataFrame(all_rounds).to_csv(output / "rounds.csv", index=False)
        pd.DataFrame(all_selections).to_csv(output / "selections.csv", index=False)
        pd.DataFrame(all_prequential).to_csv(output / "prequential.csv", index=False)
        pd.DataFrame(endpoints).to_csv(output / "seed_level_results.csv", index=False)
        pd.DataFrame(candidates).to_csv(output / "candidates.csv", index=False)
        if ablation_mode:
            pd.DataFrame(endpoints).to_csv(output / "ablations.csv", index=False)
        else:
            pd.DataFrame(columns=["policy", "seed", "confirmed_hits", "group_coverage"]).to_csv(
                output / "ablations.csv", index=False
            )
        manifest["status"] = "complete"
        manifest["runtime_seconds"] = time.perf_counter() - started
        # ru_maxrss is a process high-water mark, not isolated per-trajectory usage.
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        manifest["process_peak_memory_mb"] = peak / (1024**2 if sys.platform == "darwin" else 1024)
        manifest["memory_definition"] = (
            "process lifetime resident high-water mark, including imports; not incremental trajectory allocation"
        )
        manifest["n_trajectories"] = len(endpoints)
        manifest["files"] = {p.name: file_sha256(p) for p in sorted(output.glob("*.csv"))}
        write_json(output / "experiment_manifest.json", manifest)
        if make_report:
            from calder.evaluation.reports import build_report

            build_report(output)
        manifest["files"] = {
            str(p.relative_to(output)): file_sha256(p)
            for p in sorted(output.glob("**/*"))
            if p.is_file() and p.name != "experiment_manifest.json"
        }
        write_json(output / "experiment_manifest.json", manifest)
        return output
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        write_json(output / "experiment_manifest.json", manifest)
        raise
