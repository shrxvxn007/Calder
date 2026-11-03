"""Strict, versioned experiment configuration; unknown keys fail loudly."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class SyntheticConfig(StrictModel):
    n_candidates: int = Field(250, ge=20)
    n_features: int = Field(24, ge=3)
    n_groups: int = Field(5, ge=2)
    scenarios: list[Literal["A", "B", "C", "D"]] = ["A", "B", "C", "D"]
    data_seed: int = 317
    noise_scale: float = Field(0.2, ge=0)
    distribution_shift: float = Field(0.8, ge=0)
    group_corruption: float = Field(0.25, ge=0, le=1)
    group_overlap: float = Field(0.1, ge=0, le=1)
    outlier_fraction: float = Field(0.02, ge=0, le=0.5)
    heterogeneous_costs: bool = False


class DatasetConfig(StrictModel):
    name: str = "synthetic"
    prepared_path: str | None = None
    local_path: str | None = None
    normalization: Literal["log_normalized", "counts"] = "log_normalized"
    perturbation_key: str = "perturbation"
    control_label: str = "control"
    min_cells: int = Field(12, ge=2)
    min_genes_per_cell: int = Field(200, ge=0)
    min_cells_per_gene: int = Field(3, ge=1)
    n_hvg: int = Field(1024, ge=2)
    n_groups: int = Field(5, ge=1)
    n_components: int = Field(16, ge=1)
    chunk_size: int = Field(2048, ge=1)
    max_embedding_controls: int = Field(2048, ge=2)
    effect_weighting: Literal["inverse_variance", "uniform"] = "inverse_variance"
    bootstrap_replicates: int = Field(50, ge=2)
    feature_mode: Literal["control", "network", "combined", "uninformative"] = "control"
    network_path: str | None = None
    network_version: str | None = None
    network_mapping: str | None = None


class ModelConfig(StrictModel):
    name: str = "random_features"
    n_features: int = Field(128, ge=2)
    prior_precision: float = Field(1.0, gt=0)
    noise_floor: float = Field(0.04, gt=0)
    length_scale: float = Field(1.0, gt=0)


class AcquisitionConfig(StrictModel):
    diversity_weight: float = Field(1.0, ge=0)
    uncertainty_weight: float = Field(0.1, ge=0)
    similarity_weight: float = Field(0.1, ge=0)
    ucb_beta: float = Field(2.0, ge=0)


class ProtocolConfig(StrictModel):
    seeds: list[int] = list(range(10))
    pilot_fraction: float = Field(0.2, gt=0, lt=1)
    calibration_fraction: float = Field(0.35, ge=0, lt=1)
    calibration_min_samples: int = Field(12, ge=2)
    batch_size: int = Field(10, ge=1)
    budget_fractions: list[float] = [0.05, 0.1, 0.2]
    threshold_quantiles: list[float] = [0.7, 0.8, 0.9]
    policies: list[str] = [
        "random",
        "greedy_mean",
        "uncertainty",
        "ucb",
        "thompson",
        "probability_of_hit",
        "poh_diversity",
        "calder",
    ]
    include_oracle: bool = True

    @model_validator(mode="after")
    def valid_protocol(self) -> ProtocolConfig:
        if not self.seeds or len(set(self.seeds)) != len(self.seeds):
            raise ValueError("Seeds must be nonempty and unique")
        if not self.policies or len(set(self.policies)) != len(self.policies):
            raise ValueError("Policies must be nonempty and unique")
        if not self.budget_fractions or any(x <= 0 or x > 1 for x in self.budget_fractions):
            raise ValueError("Additional budget fractions must be in (0,1]")
        if self.budget_fractions != sorted(set(self.budget_fractions)):
            raise ValueError("Budget fractions must be strictly increasing")
        if not self.threshold_quantiles or any(not 0 < q < 1 for q in self.threshold_quantiles):
            raise ValueError("Threshold quantiles must be in (0,1)")
        return self


class ExperimentConfig(StrictModel):
    schema_version: Literal[1] = 1
    experiment_name: str = "synthetic"
    output_dir: str = "results"
    dataset: DatasetConfig = Field(default_factory=DatasetConfig)
    synthetic: SyntheticConfig = Field(default_factory=SyntheticConfig)
    model: ModelConfig = Field(default_factory=ModelConfig)
    acquisition: AcquisitionConfig = Field(default_factory=AcquisitionConfig)
    protocol: ProtocolConfig = Field(default_factory=ProtocolConfig)
    ablations: list[str] = [
        "no_calibration",
        "no_diversity",
        "no_uncertainty",
        "no_similarity",
        "random_groups",
        "uninformative_features",
    ]


def load_config(path: str | Path) -> ExperimentConfig:
    with Path(path).open() as handle:
        raw = yaml.safe_load(handle)
    return ExperimentConfig.model_validate(raw)
