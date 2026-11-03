"""Pilot-only isotonic calibration with an explicit small-sample fallback.

Labels based on noisy measured responses calibrate observed-hit events, not latent
biological truth. This limitation must be retained in reports.
"""

from __future__ import annotations

import warnings

import numpy as np
from numpy.typing import NDArray
from sklearn.isotonic import IsotonicRegression

from calder.models.uncertainty import normal_hit_probability

__all__ = ["MonotoneCalibrator", "normal_hit_probability"]


class MonotoneCalibrator:
    def __init__(self, min_samples: int = 12) -> None:
        if min_samples < 2:
            raise ValueError("At least two calibration samples are required")
        self.min_samples = min_samples
        self.is_fitted = False
        self.model = IsotonicRegression(y_min=0, y_max=1, out_of_bounds="clip")
        self.reason = "not fitted"

    def fit(self, probabilities: NDArray, labels: NDArray) -> MonotoneCalibrator:
        p, y = np.asarray(probabilities, dtype=float), np.asarray(labels, dtype=float)
        if (
            p.ndim != 1
            or p.shape != y.shape
            or not np.isfinite(p).all()
            or not np.isfinite(y).all()
            or np.any((p < 0) | (p > 1))
            or not np.isin(y, [0, 1]).all()
        ):
            raise ValueError("Calibration requires aligned finite probabilities and binary labels")
        self.is_fitted = False
        if len(p) < self.min_samples or len(np.unique(y)) < 2 or len(np.unique(p)) < 2:
            self.reason = "too few calibration examples or insufficient class/probability variation"
            warnings.warn(f"Using raw probabilities: {self.reason}", RuntimeWarning, stacklevel=2)
            return self
        self.model.fit(p, y)
        self.is_fitted = True
        self.reason = "pilot-only isotonic mapping; no exchangeability guarantee"
        return self

    def predict(self, probabilities: NDArray) -> NDArray:
        p = np.asarray(probabilities, dtype=float)
        if not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
            raise ValueError("Probabilities must be finite and lie in [0, 1]")
        return np.asarray(self.model.predict(p) if self.is_fitted else p.copy())
