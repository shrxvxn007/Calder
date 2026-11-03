"""Query-only retrospective oracle. Python privacy is an architectural, not security, boundary."""

from __future__ import annotations

import numpy as np

from calder.types import OutcomeTable


class OfflineOracle:
    __slots__ = ("__values", "__variances", "__queried")

    def __init__(self, outcomes: OutcomeTable) -> None:
        self.__values = outcomes.measured_effects.copy()
        self.__variances = outcomes.noise_variances.copy()
        self.__queried: set[int] = set()

    def query(self, indices: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ids = np.asarray(indices, dtype=int)
        if ids.ndim != 1 or len(set(ids.tolist())) != len(ids):
            raise ValueError("Oracle query must contain unique one-dimensional indices")
        if np.any(ids < 0) or np.any(ids >= len(self.__values)):
            raise IndexError("Oracle query outside candidate pool")
        if self.__queried.intersection(ids.tolist()):
            raise ValueError("Repeated oracle queries are forbidden")
        self.__queried.update(ids.tolist())
        return self.__values[ids].copy(), self.__variances[ids].copy()


def full_information_ranking(outcomes: OutcomeTable, candidate_indices: np.ndarray) -> np.ndarray:
    """Evaluation ONLY: full-response ranking; never exposed to a deployable policy."""
    truth = (
        outcomes.latent_effects
        if outcomes.latent_effects is not None
        else outcomes.measured_effects
    )
    return candidate_indices[np.lexsort((candidate_indices, -truth[candidate_indices]))]
