"""Configuration-fixed surrogate factory."""

from __future__ import annotations

from calder.models.base import Surrogate
from calder.models.bayesian_ridge import BayesianRidgeModel
from calder.models.ensemble import EnsembleModel
from calder.models.gaussian_process import GaussianProcessModel
from calder.models.mean import MeanModel
from calder.models.random_features import RandomFourierBayesianModel
from calder.models.ridge import RidgeModel


def make_model(name: str, seed: int = 0, **params: object) -> Surrogate:
    models = {
        "mean": MeanModel,
        "ridge": RidgeModel,
        "bayesian_ridge": BayesianRidgeModel,
        "rff": RandomFourierBayesianModel,
        "random_features": RandomFourierBayesianModel,
        "gaussian_process": GaussianProcessModel,
        "ensemble": EnsembleModel,
    }
    if name not in models:
        raise ValueError(f"Unknown model {name!r}; choose from {sorted(models)}")
    return models[name](seed=seed, **params)
