"""Numerically stable threshold-exceedance probabilities of latent Gaussians."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.special import ndtr


def normal_hit_probability(
    means: NDArray, epistemic_variances: NDArray, threshold: float
) -> NDArray:
    means, variance = np.asarray(means, dtype=float), np.asarray(epistemic_variances, dtype=float)
    if means.shape != variance.shape or np.any(variance < 0) or not np.isfinite(threshold):
        raise ValueError("Means/variances must align; variances must be nonnegative")
    if not np.isfinite(means).all() or not np.isfinite(variance).all():
        raise ValueError("Predictive distributions must be finite")
    result = ndtr((means - threshold) / np.sqrt(np.maximum(variance, 1e-300)))
    # A degenerate Gaussian has strict exceedance mass 0 at the threshold.
    return np.where(variance == 0, (means > threshold).astype(float), result)


probability_of_hit = normal_hit_probability
