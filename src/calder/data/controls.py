"""Control-only moments and stable, normalized inverse-variance weights."""

from __future__ import annotations

import logging

import numpy as np
from scipy import sparse

from calder.data.preprocess import iter_blocks

LOGGER = logging.getLogger(__name__)


def control_moments(adata, control_mask, normalization_status: str, chunk_size: int = 2048):
    total = np.zeros(adata.n_vars)
    squares = np.zeros(adata.n_vars)
    detected = np.zeros(adata.n_vars)
    n = 0
    next_progress = 20
    for start, stop, block in iter_blocks(adata, chunk_size, normalization_status):
        percent = int(100 * stop / adata.n_obs)
        if percent >= next_progress:
            LOGGER.info("Control-only moments: %d%% of sequential input scan", percent)
            next_progress = (percent // 20 + 1) * 20
        x = block[np.asarray(control_mask[start:stop], dtype=bool)]
        if not x.shape[0]:
            continue
        n += x.shape[0]
        total += np.asarray(x.sum(axis=0)).ravel()
        squares += np.asarray(
            x.multiply(x).sum(axis=0) if sparse.issparse(x) else (x * x).sum(axis=0)
        ).ravel()
        detected += np.asarray((x > 0).sum(axis=0)).ravel()
    if n < 2:
        raise ValueError("At least two retained untreated control cells are required")
    mean = total / n
    variance = np.maximum((squares - total * total / n) / (n - 1), 0)
    return mean, variance, detected / n, n


def inverse_variance_weights(variances, stabilizer: float = 1e-3, clip_quantile: float = 0.95):
    variances = np.asarray(variances, dtype=float)
    if stabilizer <= 0 or not 0 < clip_quantile <= 1:
        raise ValueError("Positive stabilizer and clip_quantile in (0,1] are required")
    inverse = 1.0 / (variances + stabilizer)
    inverse = np.minimum(inverse, np.quantile(inverse, clip_quantile))
    return inverse / inverse.sum()
