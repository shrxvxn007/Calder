import numpy as np
import pytest

from calder.acquisition import POLICIES, make_policy
from calder.acquisition.budget import validate_batch
from calder.types import ObservedHistory, PredictiveDistribution


@pytest.mark.parametrize("policy", POLICIES)
def test_heterogeneous_cost_feasibility(policy):
    costs = np.array([2.0, 0.5, 9.0, 1.0, 3.0])
    args = dict(
        candidate_features=np.eye(5),
        predictive_distribution=PredictiveDistribution(
            np.arange(5.0), np.ones(5), np.linspace(0.1, 0.9, 5)
        ),
        observed_history=ObservedHistory(
            np.array([], dtype=int), np.array([]), np.array([]), np.array([False])
        ),
        remaining_budget=2.0,
        batch_size=4,
        costs=costs,
        group_assignments=np.ones((5, 1)),
    )
    out = make_policy(policy).select_batch(**args)
    assert costs[out.indices].sum() <= 2 + 1e-10
    assert len(set(out.indices)) == len(out.indices)
    args["remaining_budget"] = 0.1
    assert len(make_policy(policy).select_batch(**args).indices) == 0


def test_invalid_batch_rejected():
    for batch in [np.array([0, 0]), np.array([2]), np.array([0, 1])]:
        with pytest.raises(ValueError):
            validate_batch(batch, np.ones(2), 1, 2)
