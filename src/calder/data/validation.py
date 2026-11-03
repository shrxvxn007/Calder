"""Validate actual input schemas rather than assuming external conventions."""

from __future__ import annotations

import numpy as np


def validate_adata(
    adata, perturbation_key: str = "perturbation", control_label: str = "control"
) -> None:
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("AnnData must contain cells and expression genes")
    if perturbation_key not in adata.obs:
        raise ValueError(
            f"Missing obs[{perturbation_key!r}]; available columns: {list(adata.obs.columns)}"
        )
    labels = adata.obs[perturbation_key]
    if labels.isna().any():
        raise ValueError(
            "Perturbation annotations contain missing labels; resolve these explicitly"
        )
    if not (labels.astype(str) == control_label).any():
        raise ValueError(
            f"No untreated cells with label {control_label!r}; set control_label explicitly"
        )
    if not adata.var_names.is_unique:
        raise ValueError(
            "Expression gene identifiers must be unique; resolve duplicates before preparation"
        )
    if adata.X is None:
        raise ValueError("AnnData.X is missing; supply the documented expression representation")


def validate_block(matrix, normalization_status: str) -> None:
    from scipy import sparse

    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError(
            "Expression must be finite and nonnegative (counts or log-normalized expression)"
        )
    if normalization_status == "counts" and np.any(np.abs(values - np.rint(values)) > 1e-5):
        raise ValueError("normalization_status='counts' requires integer-valued nonnegative counts")
