"""Assertions at the boundary between permitted history and hidden evaluation."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def assert_disjoint(training_indices: ArrayLike, calibration_indices: ArrayLike) -> None:
    if np.intersect1d(training_indices, calibration_indices).size:
        raise ValueError("Leakage: surrogate training and calibrator reservation overlap")


def assert_unique_selection(selected_indices: ArrayLike, previous_indices: ArrayLike) -> None:
    selected = np.asarray(selected_indices, int)
    if len(np.unique(selected)) != len(selected) or np.intersect1d(selected, previous_indices).size:
        raise ValueError("Leakage/trajectory error: a candidate was selected repeatedly")


def assert_same_decision(before: ArrayLike, after: ArrayLike) -> None:
    if not np.array_equal(before, after):
        raise AssertionError("Hidden-label permutation changed an acquisition decision")
