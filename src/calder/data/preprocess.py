"""Chunked expression access. Counts are transformed once; log input is unchanged."""

from __future__ import annotations

import logging
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from scipy import sparse

from calder.data.validation import validate_block

LOGGER = logging.getLogger(__name__)


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
        LOGGER.info("Count-scale cell QC skipped for declared log-normalized input")
        return keep
    next_progress = 20
    for start in range(0, adata.n_obs, chunk_size):
        stop = min(start + chunk_size, adata.n_obs)
        block = adata.X[start:stop]
        validate_block(block, normalization_status)
        expressed = np.asarray((block > 0).sum(axis=1)).ravel()
        keep[start:stop] = expressed >= min_expressed_genes
        percent = int(100 * stop / adata.n_obs)
        if percent >= next_progress:
            LOGGER.info("Cell QC: %d%% of cells scanned", percent)
            next_progress = (percent // 20 + 1) * 20
    return keep


@contextmanager
def response_cache(
    adata,
    response_columns,
    embedding_rows,
    embedding_columns,
    normalization_status: str,
    output_dir,
    chunk_size: int = 2048,
):
    """Sequentially cache selected response columns for dense compressed HDF5.

    Dense HDF5 row fancy indexing repeatedly decompresses whole large chunks.
    A temporary float32 disk memmap avoids that, with bounded RAM and predictable
    scratch use. Sparse files retain sparse row access and are never densified.
    Candidate embeddings are collected from controls during the same pass.
    """
    import h5py

    if not isinstance(adata.X, h5py.Dataset):
        yield (
            None,
            collect_rows(
                adata, embedding_rows, embedding_columns, normalization_status, chunk_size
            ),
        )
        return
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    needed = adata.n_obs * len(response_columns) * np.dtype("float32").itemsize
    LOGGER.info(
        "Building sequential response cache: %d columns, %.2f GB temporary disk",
        len(response_columns),
        needed / 1e9,
    )
    if shutil.disk_usage(directory).free < needed + 100_000_000:
        raise RuntimeError(
            f"Dense-backed preparation needs {needed:,} temporary bytes for selected response genes. Free disk space, reduce n_hvg, or choose an output directory on a larger disk."
        )
    with tempfile.TemporaryDirectory(prefix="calder-expression-", dir=directory) as temporary:
        cache = np.memmap(
            Path(temporary) / "responses.f32",
            mode="w+",
            dtype="float32",
            shape=(adata.n_obs, len(response_columns)),
        )
        try:
            selected = np.zeros(adata.n_obs, dtype=bool)
            selected[embedding_rows] = True
            embeddings = []
            next_progress = 20
            for start, stop, block in iter_blocks(adata, chunk_size, normalization_status):
                cache[start:stop] = block[:, response_columns]
                mask = selected[start:stop]
                if mask.any():
                    embeddings.append(block[mask][:, embedding_columns])
                percent = int(100 * stop / adata.n_obs)
                if percent >= next_progress:
                    LOGGER.info("Response cache: %d%% of cells written", percent)
                    next_progress = (percent // 20 + 1) * 20
            cache.flush()
            yield cache, np.vstack(embeddings)
        finally:
            cache.flush()
            cache._mmap.close()


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
