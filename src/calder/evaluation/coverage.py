"""Fixed-group hit coverage; group labels are association proxies."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def confirmed_group_coverage(memberships: ArrayLike, hit_mask: ArrayLike) -> int:
    groups = np.asarray(memberships, float)
    hits = np.asarray(hit_mask, bool)
    if groups.ndim != 2 or hits.shape != (len(groups),):
        raise ValueError("Membership and hit arrays must align")
    return int(np.any(groups[hits] > 0, axis=0).sum()) if hits.any() else 0


def expected_new_group_coverage(
    probabilities: ArrayLike, memberships: ArrayLike, discovered: ArrayLike | None = None
) -> float:
    """Independent-hit approximation, summed over undiscovered fixed groups."""
    p, m = np.asarray(probabilities, float), np.asarray(memberships, float)
    if m.ndim != 2 or p.shape != (len(m),):
        raise ValueError("Probabilities and membership rows must align")
    if not np.isfinite(p).all() or not np.isfinite(m).all():
        raise ValueError("Probabilities and memberships must be finite")
    if np.any((p < 0) | (p > 1)) or np.any((m < 0) | (m > 1)):
        raise ValueError("Probabilities and memberships must lie in [0, 1]")
    seen = np.zeros(m.shape[1], bool) if discovered is None else np.asarray(discovered, bool)
    if seen.shape != (m.shape[1],):
        raise ValueError("Discovered group vector does not align")
    return float((1.0 - np.prod(1.0 - m * p[:, None], axis=0))[~seen].sum())
