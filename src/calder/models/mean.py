"""Exchangeable global-mean baseline with uncertainty about its mean."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from calder.models.base import validate_training


class MeanModel:
    def __init__(self, observation_noise: float = 0.1, **_: object) -> None:
        self.observation_noise = observation_noise
        self.mean_: float | None = None

    def fit(
        self, features: NDArray, observed_effects: NDArray, noise_estimates: NDArray | None = None
    ) -> MeanModel:
        _, y, noise = validate_training(features, observed_effects, noise_estimates)
        self.mean_ = float(np.mean(y))
        self.variance_ = float(
            max(
                np.var(y, ddof=1) if len(y) > 1 else self.observation_noise**2,
                np.mean(noise),
                1e-10,
            )
            / len(y)
        )
        return self

    def predict(self, features: NDArray) -> tuple[NDArray, NDArray]:
        if self.mean_ is None:
            raise RuntimeError("Fit the surrogate before prediction")
        return np.full(len(features), self.mean_), np.full(len(features), self.variance_)

    def sample_posterior(self, features: NDArray, seed: int = 0) -> NDArray:
        if self.mean_ is None:
            raise RuntimeError("Fit the surrogate before posterior sampling")
        value = np.random.default_rng(seed).normal(self.mean_, np.sqrt(self.variance_))
        return np.full(len(features), value)
