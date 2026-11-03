"""Outcome-free candidate representations from untreated expression alone."""

from __future__ import annotations

import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler


def control_features(
    control_gene_matrix,
    means,
    variances,
    detection,
    candidate_columns,
    seed: int = 0,
    n_components: int = 16,
):
    """Rows of control_gene_matrix are cells; embeddings use genes as observations."""
    candidate_columns = np.asarray(candidate_columns, dtype=int)
    base = np.column_stack(
        (
            means[candidate_columns],
            np.log1p(variances[candidate_columns]),
            detection[candidate_columns],
        )
    )
    matrix = control_gene_matrix.T
    components = min(n_components, matrix.shape[0] - 1, matrix.shape[1] - 1)
    if components >= 1:
        embedding = TruncatedSVD(n_components=components, random_state=seed).fit_transform(matrix)
        base = np.column_stack((base, embedding))
    return StandardScaler().fit_transform(base)


def feature_clusters(features, n_groups: int = 5, seed: int = 0):
    n_groups = min(n_groups, len(features))
    if n_groups < 1:
        raise ValueError("At least one feature cluster is required")
    labels = KMeans(n_clusters=n_groups, n_init=10, random_state=seed).fit_predict(features)
    # Duplicate feature vectors can leave empty clusters. Remove them deterministically.
    occupied = sorted(set(labels.tolist()))
    remap = {value: i for i, value in enumerate(occupied)}
    labels = np.array([remap[value] for value in labels])
    return np.eye(len(occupied))[labels], tuple(
        f"feature_cluster_{i}" for i in range(len(occupied))
    )
