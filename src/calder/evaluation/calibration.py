"""Descriptive calibration metrics; no exchangeability guarantee is assumed."""

from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from scipy.stats import norm, spearmanr


def _validated(probabilities: ArrayLike, labels: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    p, y = np.asarray(probabilities, float), np.asarray(labels, float)
    if p.ndim != 1 or p.shape != y.shape:
        raise ValueError("Probabilities and labels must be aligned vectors")
    if not np.isfinite(p).all() or not np.isfinite(y).all():
        raise ValueError("Calibration inputs must be finite")
    if np.any((p < 0) | (p > 1)) or np.any((y != 0) & (y != 1)):
        raise ValueError("Probabilities must lie in [0,1]; labels must be binary")
    return p, y


def brier_score(probabilities: ArrayLike, labels: ArrayLike) -> float:
    p, y = _validated(probabilities, labels)
    return float(np.mean((p - y) ** 2)) if len(p) else float("nan")


def reliability_bins(probabilities: ArrayLike, labels: ArrayLike, n_bins: int = 10) -> pd.DataFrame:
    """Fixed equal-width bins. Empty bins have NaN estimates and zero count."""
    p, y = _validated(probabilities, labels)
    if n_bins < 1:
        raise ValueError("n_bins must be positive")
    bin_id = np.minimum((p * n_bins).astype(int), n_bins - 1)
    records = []
    for k in range(n_bins):
        mask = bin_id == k
        records.append(
            {
                "bin": k,
                "lower": k / n_bins,
                "upper": (k + 1) / n_bins,
                "count": int(mask.sum()),
                "mean_probability": float(p[mask].mean()) if mask.any() else float("nan"),
                "hit_frequency": float(y[mask].mean()) if mask.any() else float("nan"),
            }
        )
    return pd.DataFrame(records)


def expected_calibration_error(
    probabilities: ArrayLike, labels: ArrayLike, n_bins: int = 10
) -> float:
    p, y = _validated(probabilities, labels)
    if len(p) == 0:
        return float("nan")
    table = reliability_bins(p, y, n_bins)
    gap = (table.mean_probability - table.hit_frequency).abs().fillna(0)
    return float(np.sum(table["count"] / len(p) * gap))


def probability_calibration(
    probabilities: ArrayLike, labels: ArrayLike, n_bins: int = 10
) -> dict[str, float]:
    p, y = _validated(probabilities, labels)
    return {
        "n": len(p),
        "brier_score": brier_score(p, y),
        "ece": expected_calibration_error(p, y, n_bins),
        "hit_frequency": float(y.mean()) if len(y) else float("nan"),
        "mean_probability": float(p.mean()) if len(p) else float("nan"),
    }


def grouped_calibration(
    probabilities: ArrayLike,
    labels: ArrayLike,
    memberships: ArrayLike,
    group_names: list[str] | None = None,
) -> pd.DataFrame:
    p, y = _validated(probabilities, labels)
    m = np.asarray(memberships, float)
    if m.ndim != 2 or m.shape[0] != len(p):
        raise ValueError("Memberships must align with calibration examples")
    names = group_names or [str(k) for k in range(m.shape[1])]
    if len(names) != m.shape[1]:
        raise ValueError("Group name count must align with membership columns")
    return pd.DataFrame(
        [
            {"group": name, **probability_calibration(p[m[:, k] > 0], y[m[:, k] > 0])}
            for k, name in enumerate(names)
        ]
    )


def predictive_interval_coverage(
    means: ArrayLike, variances: ArrayLike, targets: ArrayLike, level: float = 0.95
) -> float:
    """Coverage of Gaussian intervals; caller must match target and variance type.

    Use epistemic variance with latent synthetic targets. For noisy measured
    targets, observation variance must be added explicitly by the caller.
    """
    mu, var, y = np.asarray(means, float), np.asarray(variances, float), np.asarray(targets, float)
    if mu.ndim != 1 or var.shape != mu.shape or y.shape != mu.shape:
        raise ValueError("Interval arrays must align")
    if not (0 < level < 1) or np.any(var < 0):
        raise ValueError("Invalid coverage level or negative variance")
    if not (np.isfinite(mu).all() and np.isfinite(var).all() and np.isfinite(y).all()):
        raise ValueError("Interval inputs must be finite")
    half_width = norm.ppf((1 + level) / 2) * np.sqrt(var)
    return float(np.mean(np.abs(mu - y) <= half_width)) if len(mu) else float("nan")


def uncertainty_error_association(variances: ArrayLike, absolute_errors: ArrayLike) -> float:
    var, err = np.asarray(variances, float), np.asarray(absolute_errors, float)
    if var.shape != err.shape or var.ndim != 1 or np.any(var < 0):
        raise ValueError("Uncertainty and error vectors must align")
    if len(var) < 3 or np.ptp(var) == 0 or np.ptp(err) == 0:
        return float("nan")
    return float(spearmanr(np.sqrt(var), err).statistic)
