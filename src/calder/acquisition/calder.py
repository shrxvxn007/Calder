"""Greedy marginal utility per cost for the CALDER and standard decision rules.

PoH is an independent implementation of threshold probability ranking, not a
reproduction of the 2026 paper's experimental protocol. Thompson uses a supplied
joint latent posterior function draw when available, and explicitly falls back to
independent marginal Gaussian draws otherwise. CALDER's pair penalty
is the sum of nonnegative cosine similarities over unordered batch pairs. This
penalty can make utility nonmonotone; greedy has no stated optimality guarantee.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from calder.acquisition.budget import validate_batch
from calder.types import BatchSelection, ObservedHistory, PredictiveDistribution


class GreedyPolicy:
    def __init__(
        self,
        name: str = "calder",
        seed: int = 0,
        diversity_weight: float = 1.0,
        uncertainty_weight: float = 0.1,
        similarity_weight: float = 0.1,
        ucb_beta: float = 2.0,
        **_: object,
    ) -> None:
        self.name = name
        self.rng = np.random.default_rng(seed)
        self.diversity_weight = float(diversity_weight)
        self.uncertainty_weight = float(uncertainty_weight)
        self.similarity_weight = float(similarity_weight)
        self.ucb_beta = float(ucb_beta)
        if (
            min(
                self.diversity_weight,
                self.uncertainty_weight,
                self.similarity_weight,
                self.ucb_beta,
            )
            < 0
        ):
            raise ValueError("Acquisition weights must be nonnegative")

    def select_batch(
        self,
        candidate_features: NDArray,
        predictive_distribution: PredictiveDistribution,
        observed_history: ObservedHistory,
        remaining_budget: float,
        batch_size: int,
        costs: NDArray,
        group_assignments: NDArray,
    ) -> BatchSelection:
        x, c, memberships = (
            np.asarray(a, dtype=float) for a in (candidate_features, costs, group_assignments)
        )
        pred = predictive_distribution
        mean = np.asarray(pred.means, dtype=float)
        var = np.asarray(pred.epistemic_variances, dtype=float)
        p = np.asarray(pred.probabilities, dtype=float)
        n = len(c)
        if x.ndim != 2 or len(x) != n or memberships.ndim != 2 or len(memberships) != n:
            raise ValueError("Candidate features, costs and groups must align")
        if any(a.shape != (n,) for a in (mean, var, p)) or np.asarray(
            observed_history.hit_groups
        ).shape != (memberships.shape[1],):
            raise ValueError("Prediction or observed-group shape does not match candidate pool")
        if not all(np.isfinite(a).all() for a in (x, c, memberships, mean, var, p)):
            raise ValueError("Acquisition arrays must be finite")
        if (
            np.any(c <= 0)
            or np.any(var < 0)
            or np.any((p < 0) | (p > 1))
            or np.any((memberships < 0) | (memberships > 1))
        ):
            raise ValueError("Invalid costs, variance, probability or membership")
        if batch_size < 0 or remaining_budget < 0 or not np.isfinite(remaining_budget):
            raise ValueError("Budget and batch size must be finite and nonnegative")
        if n == 0 or batch_size == 0:
            return BatchSelection(np.array([], dtype=int))
        std = np.sqrt(var)
        normalized_std = std / max(float(np.max(std)), 1e-12)
        remaining_group_probability = (~np.asarray(observed_history.hit_groups, dtype=bool)).astype(
            float
        )
        norms = np.linalg.norm(x, axis=1)
        unit = np.divide(x, norms[:, None], out=np.zeros_like(x), where=norms[:, None] > 0)
        similarity_sum = np.zeros(n)
        if self.name == "random":
            # A continuous priority induces a uniformly random candidate ordering.
            basic = self.rng.random(n)
        elif self.name == "greedy_mean":
            basic = mean
        elif self.name == "uncertainty":
            basic = std
        elif self.name == "ucb":
            basic = mean + self.ucb_beta * std
        elif self.name == "thompson":
            if pred.posterior_draw is not None:
                basic = np.asarray(pred.posterior_draw, dtype=float)
                if basic.shape != (n,) or not np.isfinite(basic).all():
                    raise ValueError("Posterior function draw must be finite and aligned")
            else:
                basic = mean + std * self.rng.standard_normal(n)
        else:
            basic = p
        selected: list[int] = []
        records: list[dict[str, float]] = []
        available = np.ones(n, dtype=bool)
        budget = float(remaining_budget)
        for _ in range(min(batch_size, n)):
            feasible = available & (c <= budget + 1e-10)
            if not np.any(feasible):
                break
            cover_gain = np.sum(
                memberships * p[:, None] * remaining_group_probability[None, :], axis=1
            )
            diversity = (
                self.diversity_weight * cover_gain
                if self.name in {"calder", "poh_diversity"}
                else np.zeros(n)
            )
            uncertainty = (
                self.uncertainty_weight * normalized_std if self.name == "calder" else np.zeros(n)
            )
            penalty = (
                self.similarity_weight * similarity_sum if self.name == "calder" else np.zeros(n)
            )
            gain = basic + diversity + uncertainty - penalty
            ratio = gain if self.name == "random" else gain / c
            # argmax uses first row for deterministic ties. Candidate IDs are sorted upstream.
            i = int(np.argmax(np.where(feasible, ratio, -np.inf)))
            # CALDER utility is nonmonotone; allow an empty/truncated batch if all
            # additions would lower the objective. Standard policies fill their batch.
            if self.name == "calder" and gain[i] < 0:
                break
            selected.append(i)
            records.append(
                {
                    "hit_probability": float(p[i]),
                    "coverage_gain": float(cover_gain[i]),
                    "diversity_bonus": float(diversity[i]),
                    "uncertainty_bonus": float(uncertainty[i]),
                    "similarity_penalty": float(penalty[i]),
                    "marginal_utility": float(gain[i]),
                    "utility_per_cost": float(ratio[i]),
                }
            )
            available[i] = False
            budget -= c[i]
            remaining_group_probability *= 1 - memberships[i] * p[i]
            similarity_sum += np.clip(unit @ unit[i], 0, 1)
        result = BatchSelection(np.asarray(selected, dtype=int), records)
        validate_batch(result.indices, c, remaining_budget, batch_size)
        return result
