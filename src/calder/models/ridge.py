"""Deterministic ridge baseline; it does not claim posterior uncertainty."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from sklearn.linear_model import Ridge

from calder.models.base import validate_training


class RidgeModel:
    def __init__(self, alpha: float = 1.0, observation_noise: float = 0.1, **_: object) -> None:
        self.model = Ridge(alpha=alpha)
        self.observation_noise = observation_noise
        self.fitted = False

    def fit(
        self, features: NDArray, observed_effects: NDArray, noise_estimates: NDArray | None = None
    ) -> RidgeModel:
        x, y, noise = validate_training(features, observed_effects, noise_estimates)
        self.model.fit(x, y, sample_weight=1 / np.maximum(noise + self.observation_noise**2, 1e-10))
        self.fitted = True
        return self

    def predict(self, features: NDArray) -> tuple[NDArray, NDArray]:
        if not self.fitted:
            raise RuntimeError("Fit the surrogate before prediction")
        return self.model.predict(features), np.zeros(len(features))

    def sample_posterior(self, features: NDArray, seed: int = 0) -> NDArray:
        """The deterministic ridge baseline has a degenerate predictive distribution."""
        return self.predict(features)[0]
