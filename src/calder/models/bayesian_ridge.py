"""Known-noise Gaussian linear regression using Cholesky posterior solves.

The name distinguishes this baseline from sklearn's evidence-optimized BayesianRidge:
CALDER uses explicit, configuration-fixed prior and noise parameters throughout.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import cho_factor, cho_solve, solve_triangular

from calder.models.base import validate_training


class BayesianRidgeModel:
    def __init__(
        self,
        prior_precision: float = 1.0,
        observation_noise: float = 0.1,
        fit_intercept: bool = True,
        **_: object,
    ) -> None:
        if prior_precision <= 0 or observation_noise < 0:
            raise ValueError("Prior precision must be positive and observation noise nonnegative")
        self.prior_precision = float(prior_precision)
        self.observation_noise = float(observation_noise)
        self.fit_intercept = fit_intercept
        self._factor = None

    def _design(self, x: NDArray) -> NDArray:
        return np.column_stack((np.ones(len(x)), x)) if self.fit_intercept else x

    def fit(
        self, features: NDArray, observed_effects: NDArray, noise_estimates: NDArray | None = None
    ) -> BayesianRidgeModel:
        x, y, noise = validate_training(features, observed_effects, noise_estimates)
        self.n_features_in_ = x.shape[1]
        phi = self._design(x)
        # A fixed residual/model-discrepancy floor is added to supplied sampling variance.
        variance = np.maximum(noise + self.observation_noise**2, 1e-10)
        weighted = phi / np.sqrt(variance[:, None])
        precision = weighted.T @ weighted + self.prior_precision * np.eye(phi.shape[1])
        self._factor = cho_factor(precision, lower=True, check_finite=False)
        self.coef_ = cho_solve(self._factor, phi.T @ (y / variance), check_finite=False)
        return self

    def predict(self, features: NDArray) -> tuple[NDArray, NDArray]:
        if self._factor is None:
            raise RuntimeError("Fit the surrogate before prediction")
        x = np.asarray(features, dtype=float)
        if x.ndim != 2 or x.shape[1] != self.n_features_in_ or not np.isfinite(x).all():
            raise ValueError("Prediction features have invalid shape or values")
        phi = self._design(x)
        solved = solve_triangular(self._factor[0], phi.T, lower=True, check_finite=False)
        return phi @ self.coef_, np.maximum(np.sum(solved**2, axis=0), 0.0)

    def sample_posterior(self, features: NDArray, seed: int | np.random.Generator = 0) -> NDArray:
        """Draw one coherent latent function from the Gaussian weight posterior."""
        if self._factor is None:
            raise RuntimeError("Fit the surrogate before posterior sampling")
        x = np.asarray(features, dtype=float)
        if x.ndim != 2 or x.shape[1] != self.n_features_in_ or not np.isfinite(x).all():
            raise ValueError("Prediction features have invalid shape or values")
        rng = np.random.default_rng(seed)
        perturbation = solve_triangular(
            self._factor[0],
            rng.standard_normal(len(self.coef_)),
            lower=True,
            trans="T",
            check_finite=False,
        )
        return self._design(x) @ (self.coef_ + perturbation)


BayesianRidge = BayesianRidgeModel
