from dataclasses import replace

import numpy as np

from calder.config import ExperimentConfig
from calder.experiments.active_loop import run_trajectory
from calder.experiments.oracle import OfflineOracle
from calder.experiments.simulator import simulate
from calder.experiments.splits import pilot_split
from calder.models import make_model


def small_config():
    config = ExperimentConfig()
    config.synthetic.n_candidates = 60
    config.model.n_features = 16
    config.protocol.pilot_fraction = 0.3
    config.protocol.budget_fractions = [0.1]
    config.protocol.batch_size = 6
    config.protocol.calibration_min_samples = 3
    return config


def test_shuffling_hidden_labels_cannot_change_next_selected_batch():
    config = small_config()
    pool, outcomes = simulate(config.synthetic, "B")
    train, cal = pilot_split(60, 0.3, config.protocol.calibration_fraction, 2)
    hidden = np.setdiff1d(np.arange(60), np.concatenate([train, cal]))
    changed = outcomes.measured_effects.copy()
    changed[hidden] = changed[hidden][::-1] * 100
    changed_table = replace(outcomes, measured_effects=changed)
    for name in config.protocol.policies:
        first = run_trajectory(pool, OfflineOracle(outcomes), config, 2, 0.8, name)
        second = run_trajectory(pool, OfflineOracle(changed_table), config, 2, 0.8, name)
        ids1 = [r["index"] for r in first.selections if r["round"] == 1]
        ids2 = [r["index"] for r in second.selections if r["round"] == 1]
        assert ids1 == ids2
        np.testing.assert_array_equal(
            [r["mean"] for r in first.prequential], [r["mean"] for r in second.prequential]
        )


def test_calibration_reserve_never_enters_training(monkeypatch):
    import calder.experiments.active_loop as loop

    config = small_config()
    pool, outcomes = simulate(config.synthetic, "A")
    train, cal = pilot_split(60, 0.3, config.protocol.calibration_fraction, 1)
    seen = []
    original = make_model

    def spy_factory(*args, **kwargs):
        model = original(*args, **kwargs)
        fit = model.fit

        def spy_fit(features, y, noise):
            for feature in features:
                matches = np.flatnonzero(np.all(pool.features == feature, axis=1))
                seen.extend(matches.tolist())
            return fit(features, y, noise)

        model.fit = spy_fit
        return model

    monkeypatch.setattr(loop, "make_model", spy_factory)
    run_trajectory(pool, OfflineOracle(outcomes), config, 1, 0.8, "calder")
    assert set(train) <= set(seen)
    assert not set(cal) & set(seen)
