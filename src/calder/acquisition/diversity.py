"""Expected new-group coverage under a working independent-hit approximation."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def expected_new_group_coverage(
    probabilities: NDArray, memberships: NDArray, discovered_groups: NDArray | None = None
) -> float:
    p = np.asarray(probabilities, dtype=float)
    m = np.asarray(memberships, dtype=float)
    if (
        m.ndim != 2
        or m.shape[0] != len(p)
        or np.any((m < 0) | (m > 1))
        or np.any((p < 0) | (p > 1))
    ):
        raise ValueError("Coverage requires probabilities and memberships in [0, 1]")
    undiscovered = (
        np.ones(m.shape[1], dtype=bool)
        if discovered_groups is None
        else ~np.asarray(discovered_groups, dtype=bool)
    )
    if undiscovered.shape != (m.shape[1],):
        raise ValueError("Discovered-group vector has the wrong shape")
    return float(np.sum((1 - np.prod(1 - m * p[:, None], axis=0))[undiscovered]))


def nonnegative_cosine(features: NDArray) -> NDArray:
    x = np.asarray(features, dtype=float)
    norms = np.linalg.norm(x, axis=1)
    unit = np.divide(x, norms[:, None], out=np.zeros_like(x), where=norms[:, None] > 0)
    return np.clip(unit @ unit.T, 0, 1)
