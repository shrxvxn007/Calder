"""Outcome-free acquisition interface.

Every input array contains only the currently unmeasured candidates. Returned
indices are local rows in those arrays; the runner maps them to library indices.
"""

from __future__ import annotations

from typing import Protocol

from numpy.typing import NDArray

from calder.types import BatchSelection, ObservedHistory, PredictiveDistribution


class AcquisitionPolicy(Protocol):
    def select_batch(
        self,
        candidate_features: NDArray,
        predictive_distribution: PredictiveDistribution,
        observed_history: ObservedHistory,
        remaining_budget: float,
        batch_size: int,
        costs: NDArray,
        group_assignments: NDArray,
    ) -> BatchSelection: ...
