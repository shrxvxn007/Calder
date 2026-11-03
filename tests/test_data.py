"""Analytical and H5AD integration tests for control-only preparation."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from calder.data.controls import inverse_variance_weights
from calder.data.download import download_dataset, verify_remote
from calder.data.gene_features import feature_clusters
from calder.data.network import network_features
from calder.data.perturbations import eligible_perturbations
from calder.data.pseudobulk import bootstrap_means, response_magnitude


def test_response_and_bootstrap_analytical():
    assert response_magnitude([3, 4], [0, 0], [0.5, 0.5]) == pytest.approx(np.sqrt(12.5))
    matrix = sparse.csr_matrix(np.ones((8, 3)))
    np.testing.assert_allclose(bootstrap_means(matrix, 20, np.random.default_rng(0)), 1)
    weights = inverse_variance_weights(np.array([0, 0, 1e-12]))
    assert np.isfinite(weights).all()
    assert weights.sum() == pytest.approx(1)


def test_single_gene_rules_and_unknown_groups():
    retained, excluded = eligible_perturbations(
        ["A"] * 12 + ["A+B"] * 12 + ["Z"] * 12 + ["control"], ["A", "B"], "control", 12
    )
    assert retained == ["A"]
    assert excluded == {"A+B": "ambiguous_or_combination_label", "Z": "identifier_not_in_var_names"}
    memberships, groups = feature_clusters(np.zeros((1, 2)), n_groups=5)
    assert memberships.shape == (1, 1)
    assert groups == ("feature_cluster_0",)


def make_h5ad(tmp_path, perturb_scale=1.0, missing_control=False):
    ad = pytest.importorskip("anndata")
    rng = np.random.default_rng(31)
    genes = ["A", "B", "C", "D", "E", "F"]
    labels = np.array(
        ["control"] * 35 + ["A"] * 16 + ["B"] * 16 + ["C"] * 16 + ["A+B"] * 5, dtype=object
    )
    x = rng.poisson(2.5, (len(labels), len(genes))).astype(float)
    x[labels == "A", 0] += 5 * perturb_scale
    x[labels == "B", 1] += 4 * perturb_scale
    x[labels == "C", 2] += 3 * perturb_scale
    x = np.log1p(x)
    if missing_control:
        labels[labels == "control"] = "untreated"
    data = ad.AnnData(
        sparse.csr_matrix(x),
        obs=pd.DataFrame({"perturbation": labels}, index=[f"cell{i}" for i in range(len(labels))]),
        var=pd.DataFrame(index=genes),
    )
    path = tmp_path / f"input_{perturb_scale}_{missing_control}.h5ad"
    data.write_h5ad(path)
    return path


def test_backed_preparation_roundtrip_and_hidden_outcome_isolation(tmp_path):
    from calder.data.adapters import load_prepared, prepare_dataset
    from calder.data.inspect import inspect_dataset

    first = make_h5ad(tmp_path, 1)
    second = make_h5ad(tmp_path, 100)
    kwargs = {"hvg_n": 4, "n_components": 2, "n_groups": 2, "seed": 9, "bootstrap_replicates": 12}
    with pytest.warns(UserWarning, match="fewer than 30"):
        pool, outcomes = prepare_dataset(first, tmp_path / "one", **kwargs)
    with pytest.warns(UserWarning, match="fewer than 30"):
        pool_changed, outcomes_changed = prepare_dataset(second, tmp_path / "two", **kwargs)
    np.testing.assert_array_equal(pool.features, pool_changed.features)
    np.testing.assert_array_equal(pool.memberships, pool_changed.memberships)
    assert not np.allclose(outcomes.measured_effects, outcomes_changed.measured_effects)
    assert inspect_dataset(first)["control_cells"] == 35
    saved_pool, saved_outcomes = load_prepared(tmp_path / "one")
    np.testing.assert_array_equal(saved_pool.features, pool.features)
    np.testing.assert_array_equal(saved_outcomes.measured_effects, outcomes.measured_effects)
    assert np.all(outcomes.noise_variances >= 0)
    manifest = json.loads((tmp_path / "one" / "manifest.json").read_text())
    assert manifest["pool"]["normalization_status"] == "log_normalized"
    assert manifest["pool"]["excluded_perturbations"]["A+B"] == "ambiguous_or_combination_label"
    assert "latent_effects" not in np.load(tmp_path / "one" / "outcomes.npz").files


def test_actual_schema_is_validated(tmp_path):
    from calder.data.inspect import inspect_dataset

    path = make_h5ad(tmp_path, missing_control=True)
    with pytest.raises(ValueError, match="No untreated"):
        inspect_dataset(path)
    assert inspect_dataset(path, control_label="untreated")["control_cells"] == 35


def test_dense_count_input_is_normalized_once(tmp_path):
    ad = pytest.importorskip("anndata")
    from calder.data.preprocess import collect_rows

    counts = np.array([[1, 3], [2, 2]], dtype=float)
    data = ad.AnnData(counts)
    path = tmp_path / "dense.h5ad"
    data.write_h5ad(path)
    backed = ad.read_h5ad(path, backed="r")
    try:
        expected = np.log1p(np.array([[2500, 7500], [5000, 5000]]))
        transformed = collect_rows(backed, [0, 1], [0, 1], "counts", chunk_size=1)
        np.testing.assert_allclose(transformed, expected)
        unchanged = collect_rows(backed, [0, 1], [0, 1], "log_normalized")
        np.testing.assert_array_equal(unchanged, counts)
    finally:
        backed.file.close()


def test_prepared_artifact_tampering_is_rejected(tmp_path):
    from calder.data.adapters import load_prepared, save_prepared
    from calder.types import CandidatePool, OutcomeTable

    pool = CandidatePool(np.array(["A", "B"]), np.eye(2), np.ones(2), np.eye(2), ("x", "y"))
    table = OutcomeTable(np.array([1.0, 2.0]), np.ones(2))
    save_prepared(pool, table, tmp_path)
    with (tmp_path / "outcomes.npz").open("ab") as handle:
        handle.write(b"tamper")
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_prepared(tmp_path)


def test_network_mapping_is_explicit_and_seeded(tmp_path):
    edges = tmp_path / "network.tsv"
    edges.write_text("gene1\tgene2\tscore\nA\tB\t0.9\nB\tC\t0.8\n")
    f1, m1, names, metadata = network_features(
        edges, ["A", "B", "C", "Z"], seed=2, version="toy-v1"
    )
    f2, m2, _, _ = network_features(edges, ["A", "B", "C", "Z"], seed=2, version="toy-v1")
    np.testing.assert_array_equal(f1, f2)
    np.testing.assert_array_equal(m1, m2)
    assert metadata["unmapped_genes"] == ["Z"]
    assert len(names) >= 1
    with pytest.raises(ValueError, match="version"):
        network_features(edges, ["A", "B"])


def test_remote_unknown_dataset_is_actionable():
    with pytest.raises(ValueError, match="Unknown dataset"):
        verify_remote("invented")


def test_download_enforces_verified_size_before_body(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "calder.data.download.verify_remote", lambda *_: {"size_bytes": 100, "url": "unused"}
    )
    with pytest.raises(RuntimeError, match="exceeding"):
        download_dataset("replogle22k562", tmp_path, max_bytes=10)


def test_cli_local_schema_and_rpe1_config_routing(tmp_path):
    """A synthetic H5AD fixture checks routing, never a biological benchmark."""
    from typer.testing import CliRunner

    from calder.cli import app

    path = make_h5ad(tmp_path)
    config = tmp_path / "rpe1.yaml"
    config.write_text(
        "dataset:\n"
        "  name: replogle22rpe1\n"
        "  n_hvg: 4\n"
        "  n_groups: 2\n"
        "  n_components: 2\n"
        "  bootstrap_replicates: 6\n"
        "  effect_weighting: uniform\n"
    )
    runner = CliRunner()
    inspection = runner.invoke(app, ["data", "inspect", "--path", str(path)])
    assert inspection.exit_code == 0, inspection.output
    assert json.loads(inspection.output)["control_cells"] == 35
    output = tmp_path / "prepared"
    result = runner.invoke(
        app,
        [
            "data",
            "prepare",
            "--path",
            str(path),
            "--config",
            str(config),
            "--output-dir",
            str(output),
        ],
    )
    assert result.exit_code == 0, result.output
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["pool"]["dataset"] == "replogle22rpe1"
    assert manifest["pool"]["source_url"].endswith("replogle22rpe1_processed_complete.h5ad")
    assert manifest["outcomes"]["effect_weighting"] == "uniform"
    assert manifest["pool"]["local_input"] == str(path.resolve())
