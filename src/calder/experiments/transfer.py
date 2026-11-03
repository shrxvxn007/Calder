"""Context A: rerun a frozen procedure with new sequential RPE1 observations.

This is not direct transfer of K562 fitted coefficients (context B).
"""

from pathlib import Path

from calder.config import ExperimentConfig
from calder.experiments.runner import run_benchmark
from calder.utils.serialization import content_hash, write_json


def validate_contexts(
    k562: ExperimentConfig, rpe1: ExperimentConfig, output_dir: Path
) -> tuple[Path, Path]:
    for field in ("model", "acquisition", "protocol"):
        if getattr(k562, field) != getattr(rpe1, field):
            raise ValueError(f"Context validation requires frozen {field} settings")
    if k562.dataset.name != "replogle22k562" or rpe1.dataset.name != "replogle22rpe1":
        raise ValueError("Expected K562 development and RPE1 validation configs")
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(
        output_dir / "frozen_method.json",
        {
            "hash": content_hash(
                {k: getattr(k562, k).model_dump() for k in ("model", "acquisition", "protocol")}
            ),
            "mode": "new sequential measurements per context; no direct coefficient transfer",
        },
    )
    runs = (run_benchmark(k562, output_dir / "k562"), run_benchmark(rpe1, output_dir / "rpe1"))
    from calder.evaluation.reports import build_context_comparison

    build_context_comparison(*runs, output_dir / "comparison")
    return runs
