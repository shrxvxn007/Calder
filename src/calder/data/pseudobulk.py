"""Weighted transcriptional response and deterministic cell bootstrap uncertainty."""

from __future__ import annotations

import numpy as np


def response_magnitude(perturbed_mean, control_mean, weights):
    delta = np.asarray(perturbed_mean) - np.asarray(control_mean)
    weights = np.asarray(weights, dtype=float)
    if weights.ndim != 1 or delta.shape[-1] != len(weights):
        raise ValueError("Response weights must align with the final expression dimension")
    if (
        not np.isfinite(delta).all()
        or not np.isfinite(weights).all()
        or np.any(weights < 0)
        or weights.sum() <= 0
    ):
        raise ValueError(
            "Response inputs must be finite and weights nonnegative with positive total"
        )
    return np.sqrt(np.maximum(np.sum(weights * delta**2, axis=-1), 0))


def bootstrap_means(matrix, replicates: int, rng: np.random.Generator):
    """Exact multinomial cell bootstrap without densifying expression matrices."""
    n = matrix.shape[0]
    if n == 0:
        raise ValueError("Cannot bootstrap zero cells")
    counts = rng.multinomial(n, np.full(n, 1 / n), size=replicates)
    return np.asarray(counts @ matrix) / n


def estimate_response(
    matrix, control_mean, weights, bootstrap_control_means, rng, epsilon: float = 1e-6
):
    mean = np.asarray(matrix.mean(axis=0)).ravel()
    effect = float(response_magnitude(mean, control_mean, weights))
    replicate_means = bootstrap_means(matrix, len(bootstrap_control_means), rng)
    replicates = response_magnitude(replicate_means, bootstrap_control_means, weights)
    log_replicates = np.log(replicates + epsilon)
    return {
        "effect_estimate": effect,
        "effect_standard_error": float(np.std(replicates, ddof=1)),
        "log_effect": float(np.log(effect + epsilon)),
        "log_noise_variance": float(np.var(log_replicates, ddof=1)),
    }


def energy_distance(x, y, max_cells: int = 256, seed: int = 0):
    """Optional bounded biased empirical E-distance; this is not scPerturb code."""
    from scipy import sparse
    from scipy.spatial.distance import cdist

    if max_cells < 2 or x.shape[0] == 0 or y.shape[0] == 0 or x.shape[1] != y.shape[1]:
        raise ValueError("E-distance needs nonempty aligned expression matrices and max_cells>=2")
    rng = np.random.default_rng(seed)
    x = x[rng.choice(x.shape[0], min(x.shape[0], max_cells), replace=False)]
    y = y[rng.choice(y.shape[0], min(y.shape[0], max_cells), replace=False)]
    x = x.toarray() if sparse.issparse(x) else np.asarray(x)
    y = y.toarray() if sparse.issparse(y) else np.asarray(y)
    return float(max(0, 2 * cdist(x, y).mean() - cdist(x, x).mean() - cdist(y, y).mean()))
