"""Decision rules share distributions and cannot accept hidden labels."""

import inspect

import numpy as np
import pytest

from calder.acquisition import POLICIES, make_policy
from calder.acquisition.diversity import expected_new_group_coverage
from calder.types import ObservedHistory, PredictiveDistribution


def history(k=2):
    return ObservedHistory(
        np.array([], dtype=int), np.array([]), np.array([]), np.zeros(k, dtype=bool)
    )


def inputs():
    return dict(
        candidate_features=np.eye(4),
        predictive_distribution=PredictiveDistribution(
            np.array([1.0, 2.0, 3.0, 4.0]),
            np.array([4.0, 3.0, 2.0, 1.0]),
            np.array([0.9, 0.8, 0.7, 0.6]),
        ),
        observed_history=history(),
        remaining_budget=2.0,
        batch_size=2,
        costs=np.ones(4),
        group_assignments=np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 1.0]]),
    )


@pytest.mark.parametrize("name", POLICIES)
def test_policy_shared_interface_unique_and_budget_valid(name):
    policy = make_policy(name, seed=9)
    out = policy.select_batch(**inputs())
    assert len(out.indices) == len(set(out.indices)) == 2
    assert len(out.components) == 2
    assert "oracle" not in inspect.signature(policy.select_batch).parameters


def test_baseline_rankings_and_diversity_are_sensible():
    args = inputs()
    np.testing.assert_array_equal(make_policy("greedy_mean").select_batch(**args).indices, [3, 2])
    np.testing.assert_array_equal(make_policy("uncertainty").select_batch(**args).indices, [0, 1])
    np.testing.assert_array_equal(
        make_policy("probability_of_hit").select_batch(**args).indices, [0, 1]
    )
    np.testing.assert_array_equal(
        make_policy("poh_diversity", diversity_weight=2).select_batch(**args).indices, [0, 2]
    )


def test_previously_covered_groups_receive_no_new_group_reward():
    args = inputs()
    args["observed_history"] = ObservedHistory(
        np.array([9]), np.array([3.0]), np.array([0.1]), np.array([True, False])
    )
    chosen = make_policy("poh_diversity", diversity_weight=2).select_batch(**args)
    assert chosen.indices[0] == 2
    assert chosen.components[0]["coverage_gain"] == pytest.approx(0.7)


def test_coverage_formula_including_overlap_and_empty_batch():
    assert expected_new_group_coverage(np.array([0.5, 0.5]), np.ones((2, 1))) == pytest.approx(0.75)
    assert expected_new_group_coverage(np.array([0.5, 0.5]), np.ones((2, 1)), np.array([True])) == 0
    assert expected_new_group_coverage(np.array([]), np.empty((0, 2))) == 0
    assert expected_new_group_coverage(np.array([0.5]), np.array([[0.5, 1.0]])) == 0.75


def test_empty_pool_and_missing_group_assignments():
    args = inputs()
    args.update(
        candidate_features=np.empty((0, 4)),
        predictive_distribution=PredictiveDistribution(np.array([]), np.array([]), np.array([])),
        costs=np.array([]),
        group_assignments=np.empty((0, 2)),
    )
    assert len(make_policy("calder").select_batch(**args).indices) == 0
    args = inputs()
    args.update(group_assignments=np.empty((4, 0)), observed_history=history(0))
    assert len(make_policy("calder").select_batch(**args).indices) == 2


def test_random_reproducibility_and_deterministic_ties():
    a, b = make_policy("random", seed=22), make_policy("random", seed=22)
    np.testing.assert_array_equal(
        a.select_batch(**inputs()).indices, b.select_batch(**inputs()).indices
    )
    args = inputs()
    args["predictive_distribution"] = PredictiveDistribution(
        np.zeros(4), np.ones(4), np.ones(4) * 0.5
    )
    np.testing.assert_array_equal(
        make_policy("probability_of_hit").select_batch(**args).indices, [0, 1]
    )


def test_thompson_uses_joint_posterior_draw_when_supplied():
    args = inputs()
    pred = args["predictive_distribution"]
    args["predictive_distribution"] = PredictiveDistribution(
        pred.means,
        pred.epistemic_variances,
        pred.probabilities,
        posterior_draw=np.array([3.0, 1.0, 2.0, 4.0]),
    )
    np.testing.assert_array_equal(make_policy("thompson").select_batch(**args).indices, [3, 0])
