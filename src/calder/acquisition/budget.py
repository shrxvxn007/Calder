"""Shared cost and uniqueness checks used before oracle access."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def validate_batch(
    indices: NDArray, costs: NDArray, remaining_budget: float, batch_size: int
) -> None:
    idx = np.asarray(indices)
    if idx.ndim != 1 or not np.issubdtype(idx.dtype, np.integer):
        raise ValueError("Batch indices must be an integer vector")
    if len(idx) > batch_size or len(np.unique(idx)) != len(idx):
        raise ValueError("Batch size exceeded or duplicate candidates selected")
    if np.any(idx < 0) or np.any(idx >= len(costs)):
        raise ValueError("Batch index out of bounds")
    if np.sum(costs[idx]) > remaining_budget + 1e-9:
        raise ValueError("Selected batch exceeds remaining experimental budget")
