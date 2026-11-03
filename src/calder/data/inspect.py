"""Read metadata in backed mode; no whole-matrix materialization."""

from __future__ import annotations

from pathlib import Path

from calder.data.validation import validate_adata


def inspect_dataset(
    path: str | Path, perturbation_key: str = "perturbation", control_label: str = "control"
) -> dict:
    import anndata as ad

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"H5AD not found at {path}. Download the dataset or supply --path.")
    data = ad.read_h5ad(path, backed="r")
    try:
        validate_adata(data, perturbation_key, control_label)
        labels = data.obs[perturbation_key].astype(str)
        return {
            "path": str(path.resolve()),
            "cells": data.n_obs,
            "genes": data.n_vars,
            "perturbations_excluding_control": int(labels[labels != control_label].nunique()),
            "control_cells": int((labels == control_label).sum()),
            "perturbation_key": perturbation_key,
            "control_label": control_label,
            "obs_columns": list(data.obs.columns),
            "var_columns": list(data.var.columns),
            "matrix_type": type(data.X).__name__,
            "backed": True,
            "normalization_status": "not inferred: explicitly declare counts or log_normalized at preparation",
            "perturbation_cell_counts": {str(k): int(v) for k, v in labels.value_counts().items()},
        }
    finally:
        data.file.close()
