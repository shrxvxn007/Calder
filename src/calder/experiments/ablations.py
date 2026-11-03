"""Predeclared one-factor ablations; no winner-driven retuning."""

from pathlib import Path

from calder.config import ExperimentConfig
from calder.experiments.runner import run_benchmark


def run_ablations(config: ExperimentConfig, run_dir: Path | None = None) -> Path:
    variants = ["calder", *[f"calder_{name}" for name in config.ablations]]
    allowed = {
        "no_calibration",
        "no_diversity",
        "no_uncertainty",
        "no_similarity",
        "random_groups",
        "uninformative_features",
    }
    if not set(config.ablations) <= allowed:
        raise ValueError(f"Unknown component ablation; supported: {sorted(allowed)}")
    changed = config.model_copy(deep=True)
    changed.protocol.policies = variants
    changed.protocol.include_oracle = False
    return run_benchmark(changed, run_dir=run_dir, ablation_mode=True)


def sensitivity_grid(config: ExperimentConfig) -> list[ExperimentConfig]:
    """Small declared grid, for development only; report every variant, no auto-tuning."""
    variants = []
    for field, values in {"pilot_fraction": [0.1, 0.2, 0.3], "batch_size": [5, 10, 20]}.items():
        for value in values:
            variant = config.model_copy(deep=True)
            setattr(variant.protocol, field, value)
            variant.experiment_name = f"sensitivity_{field}_{value}"
            variants.append(variant)
    for noise in [0.05, 0.2, 0.5]:
        variant = config.model_copy(deep=True)
        variant.synthetic.noise_scale = noise
        variant.experiment_name = f"sensitivity_noise_{noise}"
        variants.append(variant)
    for model in ["bayesian_ridge", "random_features", "ensemble"]:
        variant = config.model_copy(deep=True)
        variant.model.name = model
        variant.experiment_name = f"sensitivity_model_{model}"
        variants.append(variant)
    return variants
