import numpy as np
import pytest

from calder.config import ExperimentConfig, SyntheticConfig
from calder.experiments.oracle import OfflineOracle
from calder.experiments.simulator import simulate


@pytest.mark.parametrize("scenario", list("ABCD"))
def test_synthetic_known_truth_and_groups(scenario):
    config = SyntheticConfig()
    pool, labels = simulate(config, scenario)
    assert pool.features.shape == (250, 24)
    assert pool.memberships.shape == (250, 5)
    assert np.all(pool.memberships.sum(axis=1) >= 1)
    assert not np.array_equal(labels.measured_effects, labels.latent_effects)
    assert np.isfinite(labels.measured_effects).all()
    p2, y2 = simulate(config, scenario)
    np.testing.assert_array_equal(pool.features, p2.features)
    np.testing.assert_array_equal(labels.measured_effects, y2.measured_effects)
    assert not hasattr(pool, "measured_effects")


def test_easy_landscape_exact_linear_ground_truth():
    pool, outcomes = simulate(SyntheticConfig(), "A")
    expected = (
        0.3 + 0.45 * pool.features[:, 0] - 0.25 * pool.features[:, 1] + 0.15 * pool.features[:, 2]
    )
    np.testing.assert_allclose(outcomes.latent_effects, expected)


def test_noise_and_group_corruption_configurable():
    pool, outcomes = simulate(SyntheticConfig(noise_scale=0, group_overlap=0), "C")
    np.testing.assert_array_equal(outcomes.measured_effects, outcomes.latent_effects)
    assert (pool.memberships.sum(axis=1) == 1).all()


def test_oracle_copies_and_checks_queries():
    _, outcomes = simulate(SyntheticConfig())
    oracle = OfflineOracle(outcomes)
    before = outcomes.measured_effects.copy()
    outcomes.measured_effects[:] = 111
    values, _ = oracle.query(np.array([0, 1]))
    np.testing.assert_array_equal(values, before[:2])
    with pytest.raises(ValueError, match="Repeated"):
        oracle.query(np.array([0]))
    with pytest.raises(ValueError, match="unique"):
        oracle.query(np.array([2, 2]))
    with pytest.raises(IndexError):
        oracle.query(np.array([-1]))


def test_config_rejects_typo_and_invalid_budget():
    with pytest.raises(ValueError):
        ExperimentConfig.model_validate({"protcol": {}})
    with pytest.raises(ValueError):
        ExperimentConfig.model_validate({"protocol": {"budget_fractions": [0.2, 0.1]}})
