"""Optional five-member CPU MLP ensemble; disagreement is heuristic uncertainty."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from sklearn.neural_network import MLPRegressor

from calder.models.base import validate_training


class EnsembleModel:
    def __init__(
        self,
        seed: int = 0,
        members: int = 5,
        hidden_units: int = 32,
        max_iter: int = 300,
        **_: object,
    ) -> None:
        if members < 2:
            raise ValueError("Ensemble requires at least two members")
        self.seed, self.members = seed, members
        self.hidden_units, self.max_iter = hidden_units, max_iter
        self.models: list[MLPRegressor] = []

    def fit(
        self, features: NDArray, observed_effects: NDArray, noise_estimates: NDArray | None = None
    ) -> EnsembleModel:
        x, y, _ = validate_training(features, observed_effects, noise_estimates)
        self.models = []
        rng = np.random.default_rng(self.seed)
        for member in range(self.members):
            take = rng.integers(0, len(x), len(x))
            model = MLPRegressor(
                hidden_layer_sizes=(self.hidden_units,),
                random_state=self.seed + member,
                max_iter=self.max_iter,
                early_stopping=False,
            )
            model.fit(x[take], y[take])
            self.models.append(model)
        return self

    def predict(self, features: NDArray) -> tuple[NDArray, NDArray]:
        if not self.models:
            raise RuntimeError("Fit the surrogate before prediction")
        p = np.stack([m.predict(features) for m in self.models])
        return np.mean(p, axis=0), np.var(p, axis=0, ddof=1)

    def sample_posterior(self, features: NDArray, seed: int = 0) -> NDArray:
        """One bootstrap member is a heuristic function draw, not a Bayesian posterior."""
        if not self.models:
            raise RuntimeError("Fit the surrogate before posterior sampling")
        member = int(np.random.default_rng(seed).integers(len(self.models)))
        return self.models[member].predict(features)
