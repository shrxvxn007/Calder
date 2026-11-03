"""Chunked expression access. Counts are transformed once; log input is unchanged."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from calder.data.validation import validate_block


def prepare_dataset(path, output_dir, **kwargs):
    """Public preprocessing entry point; implementation lives in adapters."""
    from calder.data.adapters import prepare_dataset as prepare

    return prepare(path, output_dir, **kwargs)


def iter_blocks(adata, chunk_size: int = 2048, normalization_status: str = "log_normalized"):
    if normalization_status not in {"log_normalized", "counts"}:
        raise ValueError("Set normalization_status explicitly to 'log_normalized' or 'counts'")
    for start in range(0, adata.n_obs, chunk_size):
        stop = min(start + chunk_size, adata.n_obs)
        block = adata.X[start:stop]
        validate_block(block, normalization_status)
        if sparse.issparse(block):
            block = block.tocsr().astype(np.float64, copy=True)
        else:
            block = np.asarray(block, dtype=np.float64).copy()
        if normalization_status == "counts":
            totals = np.asarray(block.sum(axis=1)).ravel()
            scale = np.divide(1e4, totals, out=np.zeros_like(totals), where=totals > 0)
            if sparse.issparse(block):
                block = block.multiply(scale[:, None]).tocsr()
                block.data = np.log1p(block.data)
            else:
                block = np.log1p(block * scale[:, None])
        yield start, stop, block


def quality_mask(
    adata, normalization_status: str, min_expressed_genes: int = 200, chunk_size: int = 2048
):
    """Count-scale QC only. Log input is not asserted to preserve raw-count QC."""
    keep = np.ones(adata.n_obs, dtype=bool)
    if normalization_status != "counts":
        return keep
    for start in range(0, adata.n_obs, chunk_size):
        stop = min(start + chunk_size, adata.n_obs)
        block = adata.X[start:stop]
        validate_block(block, normalization_status)
        expressed = np.asarray((block > 0).sum(axis=1)).ravel()
        keep[start:stop] = expressed >= min_expressed_genes
    return keep


def collect_rows(adata, rows, columns, normalization_status: str, chunk_size: int = 2048):
    """Materialize only requested rows/columns, retaining sparse storage if possible."""
    rows = np.asarray(rows, dtype=int)
    pieces = []
    for start in range(0, len(rows), chunk_size):
        # Sorted fancy row indexing is supported by both HDF5 and backed CSR.
        selected = np.sort(rows[start : start + chunk_size])
        block = adata.X[selected]
        validate_block(block, normalization_status)
        block = (
            block.astype(np.float64, copy=True)
            if sparse.issparse(block)
            else np.asarray(block, dtype=float).copy()
        )
        if normalization_status == "counts":
            totals = np.asarray(block.sum(axis=1)).ravel()
            scale = np.divide(1e4, totals, out=np.zeros_like(totals), where=totals > 0)
            if sparse.issparse(block):
                block = block.multiply(scale[:, None]).tocsr()
                block.data = np.log1p(block.data)
            else:
                block = np.log1p(block * scale[:, None])
        pieces.append(block[:, columns])
    if not pieces:
        return sparse.csr_matrix((0, len(columns)))
    if any(sparse.issparse(piece) for piece in pieces):
        return sparse.vstack(pieces, format="csr")
    return np.vstack(pieces)
