"""Gaussian-prior Bayesian regression on seeded random Fourier features."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from calder.models.bayesian_ridge import BayesianRidgeModel


class RandomFourierBayesianModel(BayesianRidgeModel):
    def __init__(
        self, n_features: int = 128, length_scale: float = 2.0, seed: int = 0, **kwargs: object
    ) -> None:
        super().__init__(**kwargs)
        if n_features < 1 or length_scale <= 0:
            raise ValueError("Random-feature dimension and length scale must be positive")
        self.n_features = n_features
        self.length_scale = length_scale
        self.seed = seed

    def _transform(self, x: NDArray) -> NDArray:
        return np.sqrt(2 / self.n_features) * np.cos(x @ self.frequencies_ + self.phases_)

    def fit(
        self, features: NDArray, observed_effects: NDArray, noise_estimates: NDArray | None = None
    ) -> RandomFourierBayesianModel:
        x = np.asarray(features, dtype=float)
        if x.ndim != 2:
            raise ValueError("Features must be a matrix")
        self.raw_dimension_ = x.shape[1]
        rng = np.random.default_rng(self.seed)
        self.frequencies_ = rng.normal(size=(x.shape[1], self.n_features)) / self.length_scale
        self.phases_ = rng.uniform(0, 2 * np.pi, size=self.n_features)
        super().fit(self._transform(x), observed_effects, noise_estimates)
        return self

    def predict(self, features: NDArray) -> tuple[NDArray, NDArray]:
        if self._factor is None:
            raise RuntimeError("Fit the surrogate before prediction")
        x = np.asarray(features, dtype=float)
        if x.ndim != 2 or x.shape[1] != self.raw_dimension_:
            raise ValueError("Prediction feature dimension differs from training")
        return super().predict(self._transform(x))

    def sample_posterior(self, features: NDArray, seed: int | np.random.Generator = 0) -> NDArray:
        if self._factor is None:
            raise RuntimeError("Fit the surrogate before posterior sampling")
        x = np.asarray(features, dtype=float)
        if x.ndim != 2 or x.shape[1] != self.raw_dimension_:
            raise ValueError("Prediction feature dimension differs from training")
        return super().sample_posterior(self._transform(x), seed)


RandomFeaturesModel = RandomFourierBayesianModel
