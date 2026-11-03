"""Reproduce a descriptive probability-plateau audit from saved snapshots.

Example: python scripts/audit_probability_plateaus.py --run-dir
results/reference_synthetic --output reports/synthetic_probability_plateau_audit.csv
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from calder.config import ExperimentConfig
from calder.experiments.simulator import simulate
from calder.utils.serialization import content_hash, file_sha256


def audit(run_dir: Path, output: Path, quantile: float = 0.8) -> pd.DataFrame:
    config = ExperimentConfig.model_validate(json.loads((run_dir / "config.json").read_text()))
    manifest = json.loads((run_dir / "experiment_manifest.json").read_text())
    if config.dataset.name != "synthetic" or config.synthetic.heterogeneous_costs:
        raise ValueError("This audit supports the equal-cost synthetic reference run only")
    selections = pd.read_csv(run_dir / "selections.csv")
    rows = []
    for scenario in config.synthetic.scenarios:
        pool, _ = simulate(config.synthetic, scenario)
        candidate_hash = content_hash(
            {"ids": pool.gene_ids, "features": pool.features, "groups": pool.memberships}
        )
        recorded_hash = next(
            item["candidate_hash"] for item in manifest["datasets"] if item["scenario"] == scenario
        )
        if candidate_hash != recorded_hash:
            raise ValueError("Regenerated pre-outcome features do not match the source run")
        for seed in config.protocol.seeds:
            records = selections[
                (selections.scenario == scenario)
                & (selections.seed == seed)
                & (selections.policy == "calder")
                & (selections.threshold_quantile == quantile)
            ]
            snapshots = run_dir / "trajectories" / f"{scenario}_seed{seed}_q{quantile:g}_calder"
            for path in sorted(snapshots.glob("predictions_*.npz")):
                round_id = int(path.stem.split("_")[-1])
                with np.load(path) as saved:
                    raw = saved["raw_probabilities"]
                    p = saved["calibrated_probabilities"]
                    indices = saved["indices"]
                    sd = np.sqrt(saved["epistemic_variances"])
                _, counts = np.unique(p, return_counts=True)
                sd /= max(float(sd.max()), 1e-12)
                seen = np.zeros(len(pool.group_names), bool)
                for previous in records[records["round"] < round_id].itertuples():
                    if previous.is_hit:
                        seen |= pool.memberships[previous.index] > 0
                remaining_groups = (~seen).astype(float)
                memberships = pool.memberships[indices]
                x = pool.features[indices]
                norm = np.linalg.norm(x, axis=1)
                unit = np.divide(x, norm[:, None], out=np.zeros_like(x), where=norm[:, None] > 0)
                similarity = np.zeros(len(p))
                available = np.ones(len(p), bool)
                batch = records[records["round"] == round_id].sort_values("selection_order")
                tied, matched, first_ties = 0, 0, None
                for position, selected in enumerate(batch.itertuples()):
                    gain = (
                        p
                        + config.acquisition.diversity_weight
                        * np.sum(memberships * p[:, None] * remaining_groups[None, :], axis=1)
                        + config.acquisition.uncertainty_weight * sd
                        - config.acquisition.similarity_weight * similarity
                    )
                    gain[~available] = -np.inf
                    i = int(np.flatnonzero(indices == selected.index)[0])
                    max_ties = int(np.sum(gain == gain.max()))
                    if position == 0:
                        first_ties = max_ties
                    tied += max_ties > 1
                    matched += gain[i] == gain.max()
                    available[i] = False
                    remaining_groups *= 1 - memberships[i] * p[i]
                    similarity += np.clip(unit @ unit[i], 0, 1)
                rows.append(
                    {
                        "scenario": scenario,
                        "seed": seed,
                        "round": round_id,
                        "n_candidates": len(raw),
                        "raw_unique_fraction": len(np.unique(raw)) / len(raw),
                        "mapped_unique_fraction": len(counts) / len(p),
                        "raw_unique_count": len(np.unique(raw)),
                        "mapped_unique_count": len(counts),
                        "mapped_zero_fraction": np.mean(p == 0),
                        "mapped_one_fraction": np.mean(p == 1),
                        "largest_plateau_fraction": counts.max() / len(p),
                        "top_plateau_fraction": np.mean(p == p.max()),
                        "all_zero": bool(np.all(p == 0)),
                        "all_one": bool(np.all(p == 1)),
                        "first_step_max_ties": first_ties,
                        "batch_steps": len(batch),
                        "steps_with_exact_max_tie": tied,
                        "selected_objective_maxima": matched,
                    }
                )
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError("No matching saved CALDER snapshots exist")
    if not (frame.batch_steps == frame.selected_objective_maxima).all():
        raise AssertionError("Saved selection disagrees with the reconstructed batch objective")
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    metadata = {
        "source_run": str(run_dir),
        "source_manifest_sha256": file_sha256(run_dir / "experiment_manifest.json"),
        "source_config_json_sha256": file_sha256(run_dir / "config.json"),
        "source_selections_sha256": file_sha256(run_dir / "selections.csv"),
        "source_snapshot_inventory_hash": content_hash(
            {
                name: value
                for name, value in manifest.get("files", {}).items()
                if name.startswith("trajectories/")
                and f"_q{quantile:g}_calder/predictions_" in name
            }
        ),
        "source_config_hash": manifest["config_hash"],
        "source_git_commit": manifest["git_commit"],
        "threshold_quantile": quantile,
        "policy": "calder",
        "n_snapshots": len(frame),
        "csv_sha256": file_sha256(output),
        "interpretation": "Descriptive snapshot/greedy-objective audit; no causal mechanism or tuning claim",
    }
    output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    return frame


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quantile", type=float, default=0.8)
    args = parser.parse_args()
    result = audit(args.run_dir, args.output, args.quantile)
    print(result.groupby("scenario")[["batch_steps", "steps_with_exact_max_tie"]].sum().to_string())
