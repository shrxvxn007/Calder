"""Outcome-blind active loop. This module has no full-label evaluator interface."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from calder.acquisition import make_policy
from calder.config import ExperimentConfig
from calder.experiments.oracle import OfflineOracle
from calder.experiments.splits import pilot_split
from calder.models import make_model
from calder.models.calibration import MonotoneCalibrator, normal_hit_probability
from calder.types import CandidatePool, ObservedHistory, PredictiveDistribution
from calder.utils.logging import AuditLog
from calder.utils.seeds import derived_seed
from calder.utils.serialization import content_hash


@dataclass
class Trajectory:
    rounds: list[dict]
    selections: list[dict]
    prequential: list[dict]
    threshold: float
    pilot_train: np.ndarray
    pilot_calibration: np.ndarray


def run_trajectory(
    pool: CandidatePool,
    oracle: OfflineOracle,
    config: ExperimentConfig,
    seed: int,
    quantile: float,
    policy_name: str,
    artifact_dir: Path | None = None,
) -> Trajectory:
    """No hidden responses enter this API; oracle only supports selected queries.

    Threshold uses pilot TRAIN responses only. Held-out pilot calibration labels
    remain excluded from every surrogate fit. A fresh isotonic map is fitted to
    the same fixed reserve after each refit, never to adaptive acquisitions.
    """
    started = time.perf_counter()
    protocol = config.protocol
    train, cal = pilot_split(
        len(pool.gene_ids), protocol.pilot_fraction, protocol.calibration_fraction, seed
    )
    pilot = np.concatenate([train, cal])
    pilot_y, pilot_noise = oracle.query(pilot)
    threshold = float(np.quantile(pilot_y[: len(train)], quantile))
    y_by_index = dict(zip(pilot.tolist(), pilot_y.tolist(), strict=True))
    noise_by_index = dict(zip(pilot.tolist(), pilot_noise.tolist(), strict=True))
    observed = pilot.tolist()
    training = train.tolist()
    hit_groups = np.any(pool.memberships[pilot[pilot_y > threshold]] > 0, axis=0)
    audit = AuditLog(artifact_dir / "events.jsonl") if artifact_dir is not None else None
    if artifact_dir is not None:
        artifact_dir.mkdir(parents=True, exist_ok=True)
    if audit:
        audit.append(
            {
                "phase": "pilot",
                "train": train,
                "calibration": cal,
                "threshold": threshold,
                "seed": seed,
            }
        )
    selections: list[dict] = []
    prequential: list[dict] = []
    rounds: list[dict] = []
    for order, (index, y, variance) in enumerate(zip(pilot, pilot_y, pilot_noise, strict=True)):
        selections.append(
            {
                "round": 0,
                "index": int(index),
                "gene_id": str(pool.gene_ids[index]),
                "selection_order": order,
                "observed_effect": float(y),
                "noise_variance": float(variance),
                "is_hit": bool(y > threshold),
                "cost": float(pool.costs[index]),
                "groups": json.dumps(np.flatnonzero(pool.memberships[index] > 0).tolist()),
                "split": "pilot_train" if index in train else "pilot_calibration",
            }
        )
    pilot_cost = float(pool.costs[pilot].sum())
    additional_spent = 0.0
    round_id = 0

    def record_round(checkpoint: float | None = None) -> None:
        hits = sum(y_by_index[i] > threshold for i in observed)
        rounds.append(
            {
                "round": round_id,
                "threshold": threshold,
                "budget_spent": pilot_cost + additional_spent,
                "additional_budget": additional_spent,
                "budget_fraction": additional_spent / len(pool.gene_ids),
                "checkpoint_fraction": checkpoint,
                "confirmed_hits": hits,
                "group_coverage": int(hit_groups.sum()),
                "precision": hits / len(observed),
                "n_selected": len(observed),
                "active_hits": sum(y_by_index[i] > threshold for i in observed[len(pilot) :]),
                "runtime_seconds": time.perf_counter() - started,
            }
        )

    record_round(0.0)
    parameters = config.acquisition.model_dump()
    actual_policy = policy_name
    if policy_name.startswith("calder_"):
        actual_policy = "calder"
        ablation = policy_name.removeprefix("calder_")
        for name, key in (
            ("no_diversity", "diversity_weight"),
            ("no_uncertainty", "uncertainty_weight"),
            ("no_similarity", "similarity_weight"),
        ):
            if ablation == name:
                parameters[key] = 0.0
    policy = make_policy(actual_policy, seed=seed, **parameters)
    for fraction in protocol.budget_fractions:
        target = float(len(pool.gene_ids) * fraction)
        while additional_spent < target - 1e-9:
            available = np.setdiff1d(np.arange(len(pool.gene_ids)), observed)
            if not len(available):
                break
            model_params = config.model.model_dump(exclude={"name", "noise_floor"})
            model_params["observation_noise"] = np.sqrt(config.model.noise_floor)
            model = make_model(config.model.name, seed=seed, **model_params)
            train_ids = np.array(training)
            model.fit(
                pool.features[train_ids],
                np.array([y_by_index[i] for i in training]),
                np.array([noise_by_index[i] for i in training]),
            )
            means, variances = model.predict(pool.features[available])
            raw = normal_hit_probability(means, variances, threshold)
            calibrator = MonotoneCalibrator(min_samples=protocol.calibration_min_samples)
            if len(cal):
                cm, cv = model.predict(pool.features[cal])
                calibrator.fit(
                    normal_hit_probability(cm, cv, threshold),
                    np.array([y_by_index[i] > threshold for i in cal], dtype=int),
                )
            calibrated = calibrator.predict(raw)
            apply_calibration = actual_policy == "calder" and policy_name != "calder_no_calibration"
            probabilities = calibrated if apply_calibration else raw
            draw = (
                model.sample_posterior(
                    pool.features[available], seed=derived_seed(seed, f"thompson:{round_id}")
                )
                if actual_policy == "thompson"
                else None
            )
            prediction = PredictiveDistribution(means, variances, probabilities, raw, draw)
            history = ObservedHistory(
                np.array(observed),
                np.array([y_by_index[i] for i in observed]),
                np.array([noise_by_index[i] for i in observed]),
                hit_groups.copy(),
            )
            batch = policy.select_batch(
                candidate_features=pool.features[available],
                predictive_distribution=prediction,
                observed_history=history,
                remaining_budget=target - additional_spent,
                batch_size=protocol.batch_size,
                costs=pool.costs[available],
                group_assignments=pool.memberships[available],
            )
            local = np.asarray(batch.indices, dtype=int)
            if not len(local):
                break
            if (
                np.any(local < 0)
                or np.any(local >= len(available))
                or len(set(local.tolist())) != len(local)
                or len(local) > protocol.batch_size
            ):
                raise RuntimeError("Acquisition violated candidate/batch constraints")
            indices = available[local]
            if pool.costs[indices].sum() > target - additional_spent + 1e-8:
                raise RuntimeError("Acquisition exceeded budget")
            # Persist full-pool predictions BEFORE querying this round's selected labels.
            if artifact_dir is not None:
                np.savez_compressed(
                    artifact_dir / f"predictions_{round_id + 1:03d}.npz",
                    indices=available,
                    means=means,
                    epistemic_variances=variances,
                    raw_probabilities=raw,
                    calibrated_probabilities=calibrated,
                    policy_probabilities=probabilities,
                )
            if audit:
                audit.append(
                    {
                        "phase": "pre_query",
                        "round": round_id + 1,
                        "selected": indices,
                        "training_indices": train_ids,
                        "calibration_indices": cal,
                        "predictions_hash": content_hash(
                            {
                                "indices": available,
                                "mean": means,
                                "variance": variances,
                                "p": probabilities,
                            }
                        ),
                    }
                )
            values, noise = oracle.query(indices)
            round_id += 1
            for pos, (index, j, y, v) in enumerate(zip(indices, local, values, noise, strict=True)):
                additional_spent += float(pool.costs[index])
                record = {
                    "round": round_id,
                    "index": int(index),
                    "gene_id": str(pool.gene_ids[index]),
                    "selection_order": len(observed),
                    "mean": float(means[j]),
                    "epistemic_variance": float(variances[j]),
                    "raw_probability": float(raw[j]),
                    "calibrated_probability": float(calibrated[j]),
                    "probability": float(probabilities[j]),
                    "calibration_fitted": calibrator.is_fitted,
                    "observed_effect": float(y),
                    "noise_variance": float(v),
                    "model_noise_variance": config.model.noise_floor,
                    "is_hit": bool(y > threshold),
                    "cost": float(pool.costs[index]),
                    "remaining_budget": max(
                        0.0, len(pool.gene_ids) * protocol.budget_fractions[-1] - additional_spent
                    ),
                    "groups": json.dumps(np.flatnonzero(pool.memberships[index] > 0).tolist()),
                    "split": "active",
                    "threshold": threshold,
                    "measured_log_effect": float(y),
                    "absolute_error": float(abs(means[j] - y)),
                }
                if len(batch.components) > pos:
                    record.update(batch.components[pos])
                selections.append(record)
                prequential.append(record.copy())
                observed.append(int(index))
                training.append(int(index))
                y_by_index[int(index)] = float(y)
                noise_by_index[int(index)] = float(v)
                if y > threshold:
                    hit_groups |= pool.memberships[index] > 0
            if audit:
                audit.append(
                    {
                        "phase": "revealed",
                        "round": round_id,
                        "indices": indices,
                        "values": values,
                        "noise_variances": noise,
                    }
                )
            record_round(fraction if abs(additional_spent - target) < 1e-8 else None)
        if not rounds or rounds[-1].get("checkpoint_fraction") != fraction:
            # Retain every seed at each declared checkpoint, including early stopping.
            record_round(fraction)
    return Trajectory(rounds, selections, prequential, threshold, train, cal)
