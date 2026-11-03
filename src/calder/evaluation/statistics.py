"""Paired seed-level bootstrap: Monte Carlo runs are not biological replicates."""

from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike


def bootstrap_interval(
    values: ArrayLike, *, seed: int = 0, n_resamples: int = 2000, confidence: float = 0.95
) -> tuple[float, float, float]:
    x = np.asarray(values, float)
    if x.ndim != 1 or not len(x) or not np.isfinite(x).all():
        raise ValueError("Bootstrap requires a nonempty finite vector")
    if n_resamples < 1 or not 0 < confidence < 1:
        raise ValueError("Invalid bootstrap configuration")
    rng = np.random.default_rng(seed)
    # Bounded memory for large run collections.
    draws = np.empty(n_resamples)
    for start in range(0, n_resamples, 500):
        size = min(500, n_resamples - start)
        draws[start : start + size] = x[rng.integers(0, len(x), (size, len(x)))].mean(axis=1)
    alpha = (1 - confidence) / 2
    lo, hi = np.quantile(draws, [alpha, 1 - alpha])
    return float(x.mean()), float(lo), float(hi)


def paired_bootstrap(
    left: ArrayLike,
    right: ArrayLike,
    *,
    seed: int = 0,
    n_resamples: int = 2000,
    confidence: float = 0.95,
) -> dict[str, float]:
    """CI for mean(left - right), resampling whole independent paired runs."""
    a, b = np.asarray(left, float), np.asarray(right, float)
    if a.shape != b.shape or a.ndim != 1:
        raise ValueError("Paired bootstrap vectors must align")
    mean, low, high = bootstrap_interval(
        a - b, seed=seed, n_resamples=n_resamples, confidence=confidence
    )
    return {"mean_difference": mean, "ci_lower": low, "ci_upper": high, "n_pairs": len(a)}


def paired_policy_difference(
    table: pd.DataFrame,
    policy: str,
    reference: str,
    metric: str,
    pair_columns: tuple[str, ...] = ("seed",),
    **kwargs: object,
) -> dict[str, float]:
    """Require complete one-to-one pairs; never silently compare unmatched seeds."""
    columns = [*pair_columns, metric]
    a, b = table.loc[table.policy == policy, columns], table.loc[table.policy == reference, columns]
    if a.duplicated(list(pair_columns)).any() or b.duplicated(list(pair_columns)).any():
        raise ValueError("Duplicate seed-level endpoints make policy pairing ambiguous")
    merged = a.merge(b, on=list(pair_columns), how="outer", suffixes=("_a", "_b"), indicator=True)
    if merged.empty or (merged._merge != "both").any():
        raise ValueError("Policy comparison requires identical paired run identifiers")
    return paired_bootstrap(merged[f"{metric}_a"], merged[f"{metric}_b"], **kwargs)
