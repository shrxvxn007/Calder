"""Shared data contracts. CandidatePool deliberately contains no outcome labels."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class CandidatePool:
    gene_ids: NDArray
    features: NDArray
    costs: NDArray
    memberships: NDArray
    group_names: tuple[str, ...]
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        n = len(self.gene_ids)
        if n == 0 or len(set(self.gene_ids.tolist())) != n:
            raise ValueError("Candidate IDs must be nonempty and unique")
        if self.features.ndim != 2 or self.features.shape[0] != n:
            raise ValueError("Features must have one row per candidate")
        if self.costs.shape != (n,) or np.any(self.costs <= 0):
            raise ValueError("Costs must be positive and aligned")
        if self.memberships.shape != (n, len(self.group_names)):
            raise ValueError("Group memberships must align with candidates and group names")
        for a in (self.features, self.costs, self.memberships):
            if not np.isfinite(a).all():
                raise ValueError("Candidate metadata contains nonfinite values")
        if np.any((self.memberships < 0) | (self.memberships > 1)):
            raise ValueError("Membership values must lie in [0, 1]")


@dataclass(frozen=True)
class OutcomeTable:
    """Evaluator/oracle input ONLY. Effects and variances are in log space."""

    measured_effects: NDArray
    noise_variances: NDArray
    latent_effects: NDArray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.measured_effects.ndim != 1:
            raise ValueError("Outcomes must be a vector")
        if self.noise_variances.shape != self.measured_effects.shape:
            raise ValueError("Noise variance shape does not match outcomes")
        if (
            not np.isfinite(self.measured_effects).all()
            or not np.isfinite(self.noise_variances).all()
        ):
            raise ValueError("Outcomes and noise variances must be finite")
        if np.any(self.noise_variances < 0):
            raise ValueError("Noise variances cannot be negative")
        if self.latent_effects is not None:
            if (
                self.latent_effects.shape != self.measured_effects.shape
                or not np.isfinite(self.latent_effects).all()
            ):
                raise ValueError("Latent effects must be finite and aligned")


@dataclass(frozen=True)
class PredictiveDistribution:
    means: NDArray
    epistemic_variances: NDArray
    probabilities: NDArray
    raw_probabilities: NDArray | None = None
    posterior_draw: NDArray | None = None


@dataclass(frozen=True)
class ObservedHistory:
    indices: NDArray
    values: NDArray
    noise_variances: NDArray
    hit_groups: NDArray


@dataclass(frozen=True)
class BatchSelection:
    indices: NDArray
    components: list[dict[str, float]] = field(default_factory=list)
