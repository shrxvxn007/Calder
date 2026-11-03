"""Discovery metrics with explicit experimental and evaluation denominators."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def area_under_discovery(budgets: ArrayLike, hits: ArrayLike) -> float:
    """Left-continuous step integral: hits count only after their experiment.

    The first point is the pilot endpoint. Integration therefore excludes pilot
    acquisition. Units are hit × cost; dividing by the final budget span gives
    the average number of hits available during sequential acquisition.
    """
    x, y = np.asarray(budgets, float), np.asarray(hits, float)
    if x.ndim != 1 or y.shape != x.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Budget and discovery vectors must be finite and aligned")
    if np.any(np.diff(x) < 0):
        raise ValueError("Budgets must be nondecreasing")
    return float(np.sum(np.diff(x) * y[:-1]))


def discovery_metrics(
    selected_indices: ArrayLike,
    observed_values: ArrayLike,
    threshold: float,
    memberships: ArrayLike,
    total_hit_count: int,
    costs: ArrayLike | None = None,
) -> dict[str, float]:
    """Confirmed hits use measured effects; recall needs an offline denominator.

    Group coverage counts any positive fixed membership among confirmed hits.
    Fractional memberships are interpreted as group inclusion here; expected
    probabilistic group coverage belongs to the acquisition objective instead.
    """
    indices = np.asarray(selected_indices, dtype=int)
    values = np.asarray(observed_values, dtype=float)
    groups = np.asarray(memberships, dtype=float)
    if indices.ndim != 1 or values.shape != indices.shape:
        raise ValueError("Observed values must align with selected indices")
    if len(np.unique(indices)) != len(indices):
        raise ValueError("A candidate cannot count twice in discovery metrics")
    if groups.ndim != 2 or np.any(indices < 0) or np.any(indices >= len(groups)):
        raise ValueError("Membership rows do not cover selected indices")
    if not np.isfinite(values).all() or not np.isfinite(threshold):
        raise ValueError("Effects and threshold must be finite")
    if total_hit_count < 0:
        raise ValueError("Total hit count cannot be negative")
    hits = values > threshold
    count = int(hits.sum())
    if count > total_hit_count:
        raise ValueError("Selected confirmed hits exceed the offline denominator")
    covered = int(np.any(groups[indices[hits]] > 0, axis=0).sum()) if count else 0
    output = {
        "confirmed_hits": float(count),
        "precision": float(count / len(indices)) if len(indices) else 0.0,
        "recall": float(count / total_hit_count) if total_hit_count else float("nan"),
        "group_coverage": float(covered),
        "experiments": float(len(indices)),
    }
    if costs is not None:
        expense = np.asarray(costs, dtype=float)
        if expense.shape != (len(groups),) or np.any(expense <= 0):
            raise ValueError("Costs must be positive and align with the complete pool")
        output["budget_spent"] = float(expense[indices].sum())
    return output
