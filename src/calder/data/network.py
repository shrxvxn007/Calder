"""Versioned, local human STRING association networks with explicit ID mapping."""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler


def network_features(
    path: str | Path,
    gene_ids,
    seed: int = 0,
    version: str | None = None,
    mapping_path: str | Path | None = None,
    min_score: float = 0.7,
):
    """TSV columns: protein1/protein2/combined_score or gene1/gene2/score.

    STRING protein IDs require a supplied two-column protein_id/gene_id mapping.
    Gene-symbol edges require explicit gene1/gene2 headers. Network modules are
    association proxies, never asserted to be causal pathways.
    """
    if not version:
        raise ValueError("A documented STRING/network release version is required")
    edges = pd.read_csv(path, sep=r"\s+", compression="infer")
    if {"gene1", "gene2", "score"} <= set(edges.columns):
        first, second, score = "gene1", "gene2", "score"
    elif {"protein1", "protein2", "combined_score"} <= set(edges.columns):
        first, second, score = "protein1", "protein2", "combined_score"
        for column in (first, second):
            if not edges[column].astype(str).str.startswith("9606.").all():
                raise ValueError(
                    "STRING protein edges must use human taxon 9606 identifiers; supply a human-only file"
                )
        if mapping_path is None:
            raise ValueError(
                "STRING protein IDs require mapping_path with protein_id and gene_id columns"
            )
        mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)
        if not {"protein_id", "gene_id"} <= set(mapping.columns):
            raise ValueError("Mapping file needs protein_id and gene_id columns")
        if mapping.protein_id.duplicated().any():
            raise ValueError("Ambiguous duplicate protein IDs in network mapping")
        remap = dict(zip(mapping.protein_id, mapping.gene_id, strict=True))
        edges[first] = edges[first].map(remap)
        edges[second] = edges[second].map(remap)
    else:
        raise ValueError("Network file needs gene1/gene2/score or protein1/protein2/combined_score")
    scores = pd.to_numeric(edges[score], errors="raise")
    if scores.max() > 1:
        scores = scores / 1000
    if not np.isfinite(scores).all() or (scores < 0).any() or (scores > 1).any():
        raise ValueError("Association scores must be in [0,1] or STRING's [0,1000] scale")
    graph = nx.Graph()
    ids = list(map(str, gene_ids))
    available = set(edges[first].dropna()) | set(edges[second].dropna())
    graph.add_nodes_from(ids)
    for a, b, weight in zip(edges[first], edges[second], scores, strict=True):
        if a in graph and b in graph and a != b and weight >= min_score:
            graph.add_edge(a, b, weight=float(weight))
    if graph.number_of_edges():
        communities = list(nx.community.louvain_communities(graph, weight="weight", seed=seed))
    else:
        communities = [{gene} for gene in ids]
    # Missing mappings receive no group, rather than fictitious singleton modules.
    communities = [set(c) & available for c in communities]
    communities = sorted((sorted(c) for c in communities if c), key=lambda c: c[0])
    labels = {gene: i for i, community in enumerate(communities) for gene in community}
    membership = np.zeros((len(ids), len(communities)))
    for i, gene in enumerate(ids):
        if gene in labels:
            membership[i, labels[gene]] = 1
    degree = np.array(
        [[graph.degree(gene), graph.degree(gene, weight="weight")] for gene in ids], dtype=float
    )
    # Sparse adjacency SVD scales to the candidate-induced human network.
    adjacency = nx.to_scipy_sparse_array(
        graph, nodelist=ids, weight="weight", format="csr", dtype=float
    )
    from sklearn.decomposition import TruncatedSVD

    if len(ids) > 2 and graph.number_of_edges():
        embedding = TruncatedSVD(
            n_components=min(16, len(ids) - 1), random_state=seed
        ).fit_transform(adjacency)
        degree = np.column_stack((degree, embedding))
    metadata = {
        "source": str(Path(path).resolve()),
        "version": version,
        "mapping": str(mapping_path) if mapping_path else None,
        "min_score": min_score,
        "edges": graph.number_of_edges(),
        "unmapped_genes": sorted(set(ids) - available),
        "group_semantics": "STRING association network modules (candidate-induced graph)",
    }
    from calder.data.download import file_sha256

    metadata["sha256"] = file_sha256(path)
    if mapping_path is not None:
        metadata["mapping_sha256"] = file_sha256(mapping_path)
    return (
        StandardScaler().fit_transform(degree),
        membership,
        tuple(f"network_module_{i}" for i in range(len(communities))),
        metadata,
    )
