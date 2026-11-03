"""Controlled finite populations; latent response is never a candidate feature."""

from __future__ import annotations

import numpy as np

from calder.config import SyntheticConfig
from calder.types import CandidatePool, OutcomeTable


def simulate(config: SyntheticConfig, scenario: str = "A") -> tuple[CandidatePool, OutcomeTable]:
    if scenario not in "ABCD" or len(scenario) != 1:
        raise ValueError("Scenario must be A, B, C or D")
    rng = np.random.default_rng(np.random.SeedSequence([config.data_seed, ord(scenario)]))
    n, d, k = config.n_candidates, config.n_features, config.n_groups
    if k > n:
        raise ValueError("Cannot have more groups than candidates")
    labels = np.arange(n) % k
    rng.shuffle(labels)
    centers = rng.normal(0, 1.3, (k, d))
    features = centers[labels] + rng.normal(0, 0.65, (n, d))
    # A is an exactly linear log-response sanity case; B-D are misspecified for linear models.
    t = (features[:, 0] - centers[labels, 0]) / 0.65
    u = (features[:, 1] - centers[labels, 1]) / 0.65
    if scenario == "A":
        latent = 0.3 + 0.45 * features[:, 0] - 0.25 * features[:, 1] + 0.15 * features[:, 2]
    else:
        modes = 1.8 * np.exp(-0.5 * ((t - 1.0) ** 2 + (u + 0.3) ** 2))
        modes += 1.6 * np.exp(-0.5 * ((t + 1.5) ** 2 + (u - 0.5) ** 2))
        latent = -0.4 + modes + 0.25 * np.sin(features[:, 2])
        latent += 0.15 * np.cos(labels * t) + 0.1 * features[:, 3]
    shift_mask = np.zeros(n, dtype=bool)
    if scenario == "D":
        shift_mask = labels == k - 1
        features[shift_mask, :3] += config.distribution_shift
        latent[shift_mask] += config.distribution_shift * (0.5 - np.tanh(t[shift_mask]))
    outliers = rng.choice(n, size=int(config.outlier_fraction * n), replace=False)
    if scenario != "A":
        latent[outliers] += rng.normal(0, 1.5, len(outliers))
    noise_sd = np.full(n, config.noise_scale)
    if scenario in {"C", "D"}:
        noise_sd *= 0.3 + 2.5 / (1 + np.exp(-t))
    measured = latent + rng.normal(size=n) * noise_sd
    memberships = np.eye(k)[labels]
    if config.group_overlap:
        overlap = rng.random(n) < config.group_overlap
        second = (labels + rng.integers(1, k, n)) % k
        memberships[np.flatnonzero(overlap), second[overlap]] = 1
    if scenario == "D":
        corrupt = rng.random(n) < config.group_corruption
        memberships[corrupt] = np.eye(k)[rng.integers(k, size=corrupt.sum())]
    costs = rng.uniform(0.5, 1.5, n) if config.heterogeneous_costs else np.ones(n)
    metadata = {
        "dataset": "synthetic",
        "scenario": scenario,
        "data_seed": config.data_seed,
        "group_kind": "synthetic predefined groups",
        "feature_provenance": "seeded pre-outcome features",
        "cost_kind": "simulated heterogeneous" if config.heterogeneous_costs else "equal unit cost",
        "scenario_description": {
            "A": "well-specified smooth linear log-response",
            "B": "separated nonlinear hit modes",
            "C": "nonlinear modes with heteroscedastic measurements",
            "D": "corrupted groups and context shift",
        }[scenario],
    }
    pool = CandidatePool(
        np.array([f"gene_{i:04d}" for i in range(n)]),
        features,
        costs,
        memberships,
        tuple(f"group_{i}" for i in range(k)),
        metadata,
    )
    return pool, OutcomeTable(
        measured, noise_sd**2, latent, {"scale": "log_effect", "known_truth": True}
    )
