import json

import numpy as np
import pandas as pd

from calder.config import ExperimentConfig
from calder.experiments.runner import run_benchmark


def test_complete_saved_replay_and_logging(tmp_path):
    config = ExperimentConfig()
    config.synthetic.n_candidates = 40
    config.synthetic.scenarios = ["A"]
    config.model.n_features = 12
    config.protocol.seeds = [3]
    config.protocol.threshold_quantiles = [0.8]
    config.protocol.policies = ["random", "calder"]
    config.protocol.budget_fractions = [0.1, 0.2]
    config.protocol.calibration_min_samples = 2
    left = run_benchmark(config, tmp_path / "left", make_report=False)
    right = run_benchmark(config, tmp_path / "right", make_report=False)
    a = pd.read_csv(left / "selections.csv").drop(columns="run_id")
    b = pd.read_csv(right / "selections.csv").drop(columns="run_id")
    pd.testing.assert_frame_equal(a, b)
    keys = ["policy", "seed", "threshold_quantile"]
    for _, frame in a.groupby(keys):
        assert frame["index"].is_unique
    for _, frame in pd.read_csv(left / "rounds.csv").groupby(keys):
        assert frame.threshold.nunique() == 1
        assert np.all(np.diff(frame.confirmed_hits) >= 0)
        assert frame.additional_budget.max() <= 8 + 1e-8
    events_path = next((left / "trajectories").glob("*/events.jsonl"))
    events = [json.loads(line) for line in events_path.read_text().splitlines()]
    phases = [r["phase"] for r in events]
    assert phases.index("pre_query") < phases.index("revealed")
    pre = events[phases.index("pre_query")]
    assert not set(pre["training_indices"]) & set(pre["calibration_indices"])
    manifest = json.loads((left / "experiment_manifest.json").read_text())
    assert manifest["status"] == "complete"
    assert (
        manifest["config_hash"]
        == json.loads((right / "experiment_manifest.json").read_text())["config_hash"]
    )
    assert "selections.csv" in manifest["files"]
