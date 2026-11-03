"""Prepare outcome-free pools and separate, auditable offline response tables."""

from __future__ import annotations

import json
import logging
import warnings
from contextlib import ExitStack
from pathlib import Path

import numpy as np
import pandas as pd

from calder.data.controls import control_moments, inverse_variance_weights
from calder.data.download import DATASETS, file_sha256
from calder.data.gene_features import control_features, feature_clusters
from calder.data.network import network_features
from calder.data.perturbations import eligible_perturbations
from calder.data.preprocess import collect_rows, quality_mask, response_cache
from calder.data.pseudobulk import bootstrap_means, estimate_response
from calder.data.validation import validate_adata
from calder.types import CandidatePool, OutcomeTable

LOGGER = logging.getLogger(__name__)


def save_prepared(pool: CandidatePool, outcomes: OutcomeTable, output_dir: str | Path) -> None:
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        directory / "pool.npz",
        gene_ids=pool.gene_ids.astype(str),
        features=pool.features,
        costs=pool.costs,
        memberships=pool.memberships,
        group_names=np.asarray(pool.group_names),
    )
    np.savez_compressed(
        directory / "outcomes.npz",
        measured_effects=outcomes.measured_effects,
        noise_variances=outcomes.noise_variances,
    )
    (directory / "manifest.json").write_text(
        json.dumps(
            {
                "pool": pool.metadata,
                "outcomes": outcomes.metadata,
                "artifact_sha256": {
                    name: file_sha256(directory / name) for name in ("pool.npz", "outcomes.npz")
                },
            },
            indent=2,
        )
        + "\n"
    )


def load_prepared(path: str | Path) -> tuple[CandidatePool, OutcomeTable]:
    directory = Path(path)
    if not (directory / "pool.npz").is_file() or not (directory / "outcomes.npz").is_file():
        raise FileNotFoundError(
            f"Prepared pool/outcomes missing in {directory}; run calder data prepare first"
        )
    manifest = json.loads((directory / "manifest.json").read_text())
    for name, expected_hash in manifest.get("artifact_sha256", {}).items():
        if file_sha256(directory / name) != expected_hash:
            raise ValueError(f"Prepared artifact checksum mismatch: {name}; prepare the data again")
    with np.load(directory / "pool.npz", allow_pickle=False) as data:
        pool = CandidatePool(
            data["gene_ids"],
            data["features"],
            data["costs"],
            data["memberships"],
            tuple(map(str, data["group_names"])),
            manifest["pool"],
        )
    with np.load(directory / "outcomes.npz", allow_pickle=False) as data:
        outcomes = OutcomeTable(
            data["measured_effects"], data["noise_variances"], metadata=manifest["outcomes"]
        )
    if len(pool.gene_ids) != len(outcomes.measured_effects):
        raise ValueError("Prepared feature and outcome tables have inconsistent candidate counts")
    return pool, outcomes


def prepare_dataset(
    path: str | Path,
    output_dir: str | Path,
    *,
    dataset: str | None = None,
    source_url: str | None = None,
    perturbation_key: str = "perturbation",
    control_label: str = "control",
    normalization_status: str = "log_normalized",
    min_cells: int = 12,
    min_expressed_genes: int = 200,
    min_gene_cells: int = 3,
    hvg_n: int = 1024,
    n_components: int = 16,
    n_groups: int = 5,
    seed: int = 0,
    bootstrap_replicates: int = 30,
    epsilon: float = 1e-6,
    variance_stabilizer: float = 1e-3,
    weight_clip_quantile: float = 0.95,
    feature_mode: str = "control",
    network_path: str | Path | None = None,
    network_version: str | None = None,
    network_mapping: str | Path | None = None,
    max_embedding_controls: int = 2048,
    chunk_size: int = 2048,
    effect_weighting: str = "inverse_variance",
) -> tuple[CandidatePool, OutcomeTable]:
    """No transformation is fitted using intervention outcomes.

    Input counts are library-size normalized to 10,000 and log1p once. Input
    log_normalized X is left unchanged. Response genes are selected by variance
    among controls, after control-only detection filtering. Scores are weighted
    RMS expression changes, not clinical phenotypes or direct regulatory edges.
    """
    import anndata as ad

    if min_cells < 2 or bootstrap_replicates < 2 or hvg_n < 1 or epsilon <= 0:
        raise ValueError("Require min_cells>=2, bootstrap_replicates>=2, hvg_n>=1 and epsilon>0")
    if normalization_status not in {"counts", "log_normalized"}:
        raise ValueError("Declare normalization_status='counts' or 'log_normalized'")
    if effect_weighting not in {"inverse_variance", "uniform"}:
        raise ValueError("effect_weighting must be inverse_variance or uniform")
    feature_mode = {
        "control_only": "control",
        "network_only": "network",
        "none": "uninformative",
    }.get(feature_mode, feature_mode)
    if feature_mode not in {"control", "network", "combined", "uninformative"}:
        raise ValueError("feature_mode must be control, network, combined or uninformative")
    if feature_mode in {"network", "combined"} and network_path is None:
        raise ValueError(
            "Network/combined features require a local, versioned network TSV; use control otherwise"
        )
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"H5AD missing at {path}; download verified public data or supply a local path"
        )
    rng = np.random.default_rng(seed)
    data = ad.read_h5ad(path, backed="r")
    resources = ExitStack()
    try:
        validate_adata(data, perturbation_key, control_label)
        LOGGER.info(
            "Preparing %s: %d cells, %d genes; declared scale=%s",
            path,
            data.n_obs,
            data.n_vars,
            normalization_status,
        )
        labels = data.obs[perturbation_key].astype(str).to_numpy()
        LOGGER.info("Applying configured cell quality criteria in sequential chunks")
        keep = quality_mask(data, normalization_status, min_expressed_genes, chunk_size)
        controls = keep & (labels == control_label)
        LOGGER.info(
            "Cell QC retained %d/%d cells (%d controls); estimating control-only moments",
            int(keep.sum()),
            data.n_obs,
            int(controls.sum()),
        )
        means, variances, detection, n_controls = control_moments(
            data, controls, normalization_status, chunk_size
        )
        candidates, exclusions = eligible_perturbations(
            labels[keep], data.var_names, control_label, min_cells
        )
        if not candidates:
            raise ValueError(
                "No eligible single-gene perturbations remain; inspect IDs, labels and min_cells"
            )
        columns = data.var_names.get_indexer(candidates)
        eligible_genes = np.flatnonzero(detection * n_controls >= min_gene_cells - 1e-9)
        if not len(eligible_genes):
            raise ValueError("No expression genes pass control-only detection filtering")
        hvg_columns = eligible_genes[np.argsort(-variances[eligible_genes], kind="stable")[:hvg_n]]
        LOGGER.info(
            "Selected %d candidate perturbations and %d control-selected response genes; excluded %d perturbation labels",
            len(candidates),
            len(hvg_columns),
            len(exclusions),
        )
        # All dimensionality reduction uses untreated controls, including candidate genes.
        control_rows = np.flatnonzero(controls)
        embedding_rows = np.sort(
            rng.choice(control_rows, min(max_embedding_controls, n_controls), replace=False)
        )
        cached_responses, embedding_matrix = resources.enter_context(
            response_cache(
                data,
                hvg_columns,
                embedding_rows,
                columns,
                normalization_status,
                output_dir,
                chunk_size,
            )
        )
        LOGGER.info(
            "Fitting control-only candidate features and fixed groups from %d embedding controls",
            len(embedding_rows),
        )
        features = control_features(
            embedding_matrix, means, variances, detection, columns, seed, n_components
        )
        membership, group_names = feature_clusters(features, n_groups, seed)
        network_metadata = None
        if network_path is not None:
            net_features, membership, group_names, network_metadata = network_features(
                network_path, candidates, seed, network_version, network_mapping
            )
            if feature_mode == "network":
                features = net_features
            elif feature_mode == "combined":
                features = np.column_stack((features, net_features))
        if feature_mode == "uninformative":
            features = rng.normal(size=features.shape)
        control_matrix = (
            np.asarray(cached_responses[control_rows], dtype=float)
            if cached_responses is not None
            else collect_rows(data, control_rows, hvg_columns, normalization_status, chunk_size)
        )
        bootstrap_controls = bootstrap_means(control_matrix, bootstrap_replicates, rng)
        LOGGER.info(
            "Estimating perturbation responses with %d bootstrap replicates each",
            bootstrap_replicates,
        )
        response_mean = means[hvg_columns]
        if effect_weighting == "inverse_variance":
            weights = inverse_variance_weights(
                variances[hvg_columns], variance_stabilizer, weight_clip_quantile
            )
        else:
            weights = np.full(len(hvg_columns), 1 / len(hvg_columns))
        table = []
        for gene_index, gene in enumerate(candidates, start=1):
            rows = np.flatnonzero(keep & (labels == gene))
            matrix = (
                np.asarray(cached_responses[rows], dtype=float)
                if cached_responses is not None
                else collect_rows(data, rows, hvg_columns, normalization_status, chunk_size)
            )
            estimate = estimate_response(
                matrix, response_mean, weights, bootstrap_controls, rng, epsilon
            )
            table.append(
                {
                    "gene_id": gene,
                    "cell_count": len(rows),
                    **estimate,
                    "quality_status": "low_cell_count" if len(rows) < 30 else "pass",
                    "source": source_url or DATASETS.get(dataset, str(path.resolve())),
                }
            )
            if gene_index % 250 == 0 or gene_index == len(candidates):
                LOGGER.info(
                    "Response bootstrap: %d/%d perturbations complete", gene_index, len(candidates)
                )
        frame = pd.DataFrame(table)
        if (frame.cell_count < 30).any():
            warnings.warn(
                f"{int((frame.cell_count < 30).sum())} perturbations have fewer than 30 cells; bootstrap uncertainty may be unstable",
                stacklevel=2,
            )
        manifest = {
            "dataset": dataset or path.stem,
            "source_url": source_url or DATASETS.get(dataset),
            "local_input": str(path.resolve()),
            "sha256": file_sha256(path),
            "processing_version": "calder-0.1.0",
            "seed": seed,
            "original_cells": data.n_obs,
            "original_genes": data.n_vars,
            "retained_control_cells": n_controls,
            "retained_perturbation_cells": int(frame.cell_count.sum()),
            "excluded_cells": int(data.n_obs - n_controls - frame.cell_count.sum()),
            "retained_candidates": len(candidates),
            "excluded_perturbations": exclusions,
            "normalization_status": normalization_status,
            "transformation": "library size 10000 + log1p"
            if normalization_status == "counts"
            else "none: already log normalized",
            "qc": {
                "min_cells": min_cells,
                "min_gene_cells_in_controls": min_gene_cells,
                "min_expressed_genes_per_cell": min_expressed_genes
                if normalization_status == "counts"
                else None,
                "log_input_count_scale_qc": "skipped; raw counts unavailable"
                if normalization_status != "counts"
                else "applied",
            },
            "feature_selection": "highest control-only expression variance after control-only detection filtering",
            "response_genes": list(map(str, data.var_names[hvg_columns])),
            "response_gene_count": len(hvg_columns),
            "response_weights": weights.tolist(),
            "feature_mode": feature_mode,
            "feature_fit_source": "untreated controls and optional fixed external annotations only",
            "embedding_control_cells": len(embedding_rows),
            "embedding_components": n_components,
            "group_semantics": "artificial control-feature clusters"
            if network_metadata is None
            else "STRING association network modules",
            "network": network_metadata,
            "expression_access": "sequential dense-backed chunks with temporary float32 response-column memmap"
            if cached_responses is not None
            else "sparse backed row chunks; sparse storage preserved",
            "temporary_cache_bytes": data.n_obs * len(hvg_columns) * 4
            if cached_responses is not None
            else 0,
        }
        download_manifest = path.with_suffix(".download.json")
        if download_manifest.exists():
            manifest["download_provenance"] = json.loads(download_manifest.read_text())
        outcome_metadata = {
            "scale": "log(effect+epsilon)",
            "epsilon": epsilon,
            "effect_weighting": effect_weighting,
            "effect_definition": "sqrt(sum_j weights_j*(perturbation_mean_j-control_mean_j)^2)",
            "bootstrap_replicates": bootstrap_replicates,
            "bootstrap": "cell-level multinomial resampling of perturbation and control cells",
            "noise_estimates": "sample variance of log bootstrap response magnitudes",
            "shared_control_uncertainty": "common control bootstrap induces cross-candidate dependence; diagonal surrogate approximation ignores it",
            "effect_interpretation": "within-dataset transcriptional response, not therapeutic usefulness or direct regulatory edges",
            "sha256": manifest["sha256"],
            "dataset": manifest["dataset"],
        }
        pool = CandidatePool(
            np.asarray(candidates),
            features,
            np.ones(len(candidates)),
            membership,
            group_names,
            manifest,
        )
        outcomes = OutcomeTable(
            frame.log_effect.to_numpy(),
            frame.log_noise_variance.to_numpy(),
            metadata=outcome_metadata,
        )
        save_prepared(pool, outcomes, output_dir)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        frame.to_csv(Path(output_dir) / "perturbations.csv", index=False)
        LOGGER.info(
            "Prepared %d perturbation cells plus %d controls, %d response genes and %d candidates in %s; excluded %d cells",
            int(frame.cell_count.sum()),
            n_controls,
            len(hvg_columns),
            len(candidates),
            output_dir,
            manifest["excluded_cells"],
        )
        return pool, outcomes
    finally:
        try:
            resources.close()
        finally:
            data.file.close()
