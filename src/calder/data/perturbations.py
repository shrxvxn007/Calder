"""Conservative single-gene inclusion rules; no heuristic combination splitting."""

from __future__ import annotations

import re
from collections import Counter


def eligible_perturbations(labels, gene_names, control_label: str, min_cells: int):
    genes = set(map(str, gene_names))
    counts = Counter(map(str, labels))
    retained: list[str] = []
    excluded: dict[str, str] = {}
    for label, count in sorted(counts.items()):
        if label == control_label:
            continue
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", label):
            excluded[label] = "ambiguous_or_combination_label"
        elif label not in genes:
            excluded[label] = "identifier_not_in_var_names"
        elif count < min_cells:
            excluded[label] = f"fewer_than_{min_cells}_retained_cells"
        else:
            retained.append(label)
    return retained, excluded
