"""Uniform random pilots shared across policies; calibration rows never train the surrogate."""

from __future__ import annotations

import numpy as np


def pilot_split(
    n: int, fraction: float, calibration_fraction: float, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    size = min(n - 1, max(3, int(np.ceil(n * fraction))))
    indices = np.random.default_rng(seed).permutation(n)[:size]
    n_cal = min(size - 2, int(np.floor(size * calibration_fraction)))
    return indices[n_cal:], indices[:n_cal]
