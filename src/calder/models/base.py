"""Surrogates predict latent log effects, with epistemic variance only."""

from __future__ import annotations

from typing import Protocol

import numpy as np
from numpy.typing import NDArray


class Surrogate(Protocol):
    def fit(
        self, features: NDArray, observed_effects: NDArray, noise_estimates: NDArray | None = None
    ) -> Surrogate: ...
    def predict(self, features: NDArray) -> tuple[NDArray, NDArray]: ...

    def sample_posterior(self, features: NDArray, seed: int = 0) -> NDArray: ...


def validate_training(
    features: NDArray, values: NDArray, noise: NDArray | None
) -> tuple[NDArray, NDArray, NDArray]:
    x, y = np.asarray(features, dtype=float), np.asarray(values, dtype=float)
    if x.ndim != 2 or y.shape != (len(x),) or len(x) == 0:
        raise ValueError("Training features and outcomes must be aligned and nonempty")
    n = np.zeros(len(y)) if noise is None else np.asarray(noise, dtype=float)
    if n.shape != y.shape or np.any(n < 0):
        raise ValueError("Observation noise must be an aligned vector of nonnegative variances")
    if not all(np.isfinite(a).all() for a in (x, y, n)):
        raise ValueError("Training arrays must be finite")
    return x, y, n
