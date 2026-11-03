"""Optional exact GP for small candidate libraries, without test-label tuning."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel

from calder.models.base import validate_training


class GaussianProcessModel:
    def __init__(
        self,
        length_scale: float = 2.0,
        prior_precision: float = 1.0,
        observation_noise: float = 0.1,
        seed: int = 0,
        **_: object,
    ) -> None:
        self.length_scale, self.prior_precision = length_scale, prior_precision
        self.observation_noise, self.seed = observation_noise, seed
        self.model = None

    def fit(
        self, features: NDArray, observed_effects: NDArray, noise_estimates: NDArray | None = None
    ) -> GaussianProcessModel:
        x, y, noise = validate_training(features, observed_effects, noise_estimates)
        if len(x) > 3000:
            raise ValueError("Exact GP is limited to 3,000 observations; use random features")
        self.model = GaussianProcessRegressor(
            kernel=ConstantKernel(1 / self.prior_precision, constant_value_bounds="fixed")
            * RBF(self.length_scale, length_scale_bounds="fixed"),
            alpha=np.maximum(noise + self.observation_noise**2, 1e-10),
            optimizer=None,
            normalize_y=False,
            random_state=self.seed,
        )
        self.model.fit(x, y)
        return self

    def predict(self, features: NDArray) -> tuple[NDArray, NDArray]:
        if self.model is None:
            raise RuntimeError("Fit the surrogate before prediction")
        mean, std = self.model.predict(features, return_std=True)
        return mean, std**2

    def sample_posterior(self, features: NDArray, seed: int = 0) -> NDArray:
        """Joint latent GP draw; alpha is training noise and is not added here."""
        if self.model is None:
            raise RuntimeError("Fit the surrogate before posterior sampling")
        return self.model.sample_y(features, n_samples=1, random_state=seed).ravel()
