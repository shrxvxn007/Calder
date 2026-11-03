"""Analytical checks distinguish epistemic variance from measurement noise."""

import numpy as np
import pytest

from calder.models import make_model
from calder.models.bayesian_ridge import BayesianRidgeModel
from calder.models.uncertainty import normal_hit_probability


def test_scalar_posterior_matches_closed_form():
    x = np.array([[1.0], [2.0]])
    y = np.array([1.0, 3.0])
    noise = np.array([0.25, 1.0])
    model = BayesianRidgeModel(prior_precision=2, observation_noise=0, fit_intercept=False)
    model.fit(x, y, noise)
    mean, variance = model.predict(np.array([[3.0]]))
    precision = 2 + 1 / 0.25 + 4
    expected_coef = (1 / 0.25 + 6) / precision
    np.testing.assert_allclose(mean, [3 * expected_coef])
    np.testing.assert_allclose(variance, [9 / precision])


def test_sampling_noise_is_not_added_to_latent_predictive_variance():
    model = BayesianRidgeModel(prior_precision=1, observation_noise=0, fit_intercept=False)
    model.fit(np.array([[1.0]]), np.array([2.0]), np.array([4.0]))
    _, variance = model.predict(np.array([[1.0]]))
    np.testing.assert_allclose(variance, [0.8])


@pytest.mark.parametrize("name", ["mean", "ridge", "bayesian_ridge", "rff", "gaussian_process"])
def test_model_stable_with_degenerate_and_near_zero_noise_data(name):
    x = np.ones((5, 3))
    model = make_model(name, seed=4, observation_noise=0)
    model.fit(x, np.full(5, 2.0), np.zeros(5))
    means, variance = model.predict(x)
    assert np.isfinite(means).all() and np.isfinite(variance).all()
    assert np.min(variance) >= 0
    means, variance = (
        model.predict(np.empty((0, 3))) if name not in {"ridge", "gaussian_process"} else ([], [])
    )
    assert len(means) == len(variance) == 0


def test_rff_seed_and_posterior_variance_contract():
    rng = np.random.default_rng(41)
    x = rng.normal(size=(15, 4))
    y = np.sin(x[:, 0])
    a, b = make_model("rff", seed=6), make_model("rff", seed=6)
    a.fit(x[:5], y[:5], np.full(5, 0.01))
    b.fit(x[:5], y[:5], np.full(5, 0.01))
    for left, right in zip(a.predict(x), b.predict(x), strict=True):
        np.testing.assert_array_equal(left, right)
    variance_before = a.predict(x)[1]
    a.fit(x, y, np.full(15, 0.01))
    assert np.all(a.predict(x)[1] <= variance_before + 1e-10)


def test_hit_probabilities_normal_and_zero_variance():
    p = normal_hit_probability(np.array([0.0, 1.0, -1.0]), np.ones(3), 0)
    np.testing.assert_allclose(p, [0.5, 0.841344746, 0.158655254])
    np.testing.assert_array_equal(
        normal_hit_probability(np.array([-1, 0, 1]), np.zeros(3), 0), [0, 0, 1]
    )
    with pytest.raises(ValueError):
        normal_hit_probability(np.zeros(2), np.array([-1, 0]), 0)


def test_unfit_and_invalid_data_fail_loudly():
    model = make_model("rff")
    with pytest.raises(RuntimeError):
        model.predict(np.zeros((2, 3)))
    with pytest.raises(ValueError):
        model.fit(np.zeros((2, 3)), np.zeros(2), np.array([0, -1]))


def test_joint_weight_draw_has_analytical_mean_variance_and_shared_values():
    model = BayesianRidgeModel(prior_precision=2, observation_noise=0, fit_intercept=False)
    model.fit(np.array([[1.0], [2.0]]), np.array([1.0, 3.0]), np.array([0.25, 1.0]))
    x = np.array([[1.0], [1.0], [2.0]])
    draws = np.stack([model.sample_posterior(x, seed=s) for s in range(3000)])
    # Identical feature vectors must share a function value, not independent draws.
    np.testing.assert_array_equal(draws[:, 0], draws[:, 1])
    np.testing.assert_allclose(draws[:, 2], 2 * draws[:, 0])
    means, variances = model.predict(x)
    np.testing.assert_allclose(np.mean(draws, axis=0), means, atol=0.04)
    np.testing.assert_allclose(np.var(draws, axis=0), variances, rtol=0.08)


@pytest.mark.parametrize("name", ["bayesian_ridge", "rff", "gaussian_process", "mean", "ridge"])
def test_joint_posterior_draw_is_seeded(name):
    x = np.array([[0.0], [1.0], [2.0]])
    model = make_model(name, seed=9)
    model.fit(x, np.array([1.0, 0.0, -1.0]), np.full(3, 0.02))
    np.testing.assert_array_equal(
        model.sample_posterior(x, seed=7), model.sample_posterior(x, seed=7)
    )
