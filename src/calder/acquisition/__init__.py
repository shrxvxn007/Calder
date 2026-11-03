"""Standard acquisition policies share a common outcome-free implementation."""

from __future__ import annotations

from calder.acquisition.calder import GreedyPolicy

POLICIES = (
    "random",
    "greedy_mean",
    "uncertainty",
    "ucb",
    "thompson",
    "probability_of_hit",
    "poh_diversity",
    "calder",
)


def make_policy(name: str, seed: int = 0, **params: object) -> GreedyPolicy:
    if name not in POLICIES:
        raise ValueError(f"Unknown policy {name!r}; choose from {POLICIES}")
    return GreedyPolicy(name=name, seed=seed, **params)
