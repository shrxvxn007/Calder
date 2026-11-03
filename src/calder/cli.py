"""Command-line entry points with actionable research workflow errors."""

from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Annotated

import typer

from calder.config import ExperimentConfig, load_config
from calder.utils.serialization import json_default

app = typer.Typer(
    help="CALDER: reproducible, retrospective perturbation selection.", no_args_is_help=True
)
data_app = typer.Typer(help="Inspect, verify, download and prepare public or local AnnData.")
app.add_typer(data_app, name="data")


def _error(exc: Exception) -> None:
    typer.secho(f"Error: {exc}", fg=typer.colors.RED, err=True)
    raise typer.Exit(1)


@data_app.command("inspect")
def data_inspect(
    dataset: str = "replogle22k562",
    path: Path | None = None,
    perturbation_key: str = "perturbation",
    control_label: str = "control",
) -> None:
    """Inspect a local H5AD schema; without --path verify remote metadata only."""
    from calder.data.download import verify_remote
    from calder.data.inspect import inspect_dataset

    try:
        info = (
            inspect_dataset(path, perturbation_key=perturbation_key, control_label=control_label)
            if path
            else verify_remote(dataset)
        )
        typer.echo(json.dumps(info, indent=2, default=json_default))
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        _error(exc)


@data_app.command("download")
def data_download(
    dataset: str = "replogle22k562", destination: Path = Path("data/raw"), max_gb: float = 2.0
) -> None:
    """Verify metadata and stream a bounded download; never substitute synthetic data."""
    from calder.data.download import download_dataset

    try:
        typer.echo(str(download_dataset(dataset, destination, max_bytes=int(max_gb * 1e9))))
    except (OSError, ValueError, RuntimeError) as exc:
        _error(exc)


@data_app.command("prepare")
def data_prepare(
    dataset: str | None = None,
    path: Path | None = None,
    config: Path | None = None,
    output_dir: Path | None = None,
) -> None:
    """Prepare control-only features and separately held response tables."""
    from calder.data.download import dataset_url
    from calder.data.preprocess import prepare_dataset

    try:
        settings = load_config(config) if config else ExperimentConfig()
        ds = settings.dataset
        resolved_dataset = dataset or (ds.name if ds.name != "synthetic" else "replogle22k562")
        local = path or (
            Path(ds.local_path)
            if ds.local_path
            else Path(f"data/raw/{resolved_dataset}_processed_complete.h5ad")
        )
        target = output_dir or Path(f"data/processed/{resolved_dataset}")
        pool, _ = prepare_dataset(
            local,
            target,
            dataset=resolved_dataset,
            source_url=dataset_url(resolved_dataset),
            perturbation_key=ds.perturbation_key,
            control_label=ds.control_label,
            normalization_status=ds.normalization,
            min_cells=ds.min_cells,
            min_expressed_genes=ds.min_genes_per_cell,
            min_gene_cells=ds.min_cells_per_gene,
            hvg_n=ds.n_hvg,
            n_groups=ds.n_groups,
            n_components=ds.n_components,
            chunk_size=ds.chunk_size,
            max_embedding_controls=ds.max_embedding_controls,
            effect_weighting=ds.effect_weighting,
            seed=settings.synthetic.data_seed,
            bootstrap_replicates=ds.bootstrap_replicates,
            feature_mode=ds.feature_mode,
            network_path=ds.network_path,
            network_version=ds.network_version,
            network_mapping=ds.network_mapping,
        )
        typer.echo(f"Prepared {len(pool.gene_ids)} candidates in {target}")
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        _error(exc)


def _run(
    config_path: Path, run_dir: Path | None, synthetic_only: bool = False, ablate: bool = False
) -> None:
    try:
        config = load_config(config_path)
        if synthetic_only and config.dataset.name != "synthetic":
            raise ValueError("'synthetic' requires dataset.name: synthetic")
        if ablate:
            from calder.experiments.ablations import run_ablations

            result = run_ablations(config, run_dir)
        else:
            from calder.experiments.runner import run_benchmark

            result = run_benchmark(config, run_dir)
        typer.echo(f"Saved reproducible experiment: {result.resolve()}")
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        _error(exc)


@app.command()
def synthetic(
    config: Annotated[Path, typer.Option()] = Path("configs/synthetic.yaml"),
    run_dir: Path | None = None,
) -> None:
    """Run the known-truth finite-pool benchmark."""
    _run(config, run_dir, synthetic_only=True)


@app.command()
def benchmark(config: Annotated[Path, typer.Option()], run_dir: Path | None = None) -> None:
    """Compare policies with identical pilots, thresholds, surrogate and budgets."""
    _run(config, run_dir)


@app.command()
def ablate(
    config: Annotated[Path, typer.Option()] = Path("configs/ablations.yaml"),
    run_dir: Path | None = None,
    sensitivity: bool = False,
) -> None:
    """Run all declared components or the small one-factor sensitivity grid."""
    if not sensitivity:
        _run(config, run_dir, ablate=True)
        return
    from calder.experiments.ablations import sensitivity_grid
    from calder.experiments.runner import run_benchmark

    try:
        settings = load_config(config)
        for variant in sensitivity_grid(settings):
            typer.echo(
                str(
                    run_benchmark(variant, (run_dir / variant.experiment_name) if run_dir else None)
                )
            )
    except (OSError, ValueError, RuntimeError) as exc:
        _error(exc)


@app.command()
def report(run_dir: Annotated[Path, typer.Option()]) -> None:
    """Regenerate every available figure and summary from saved results."""
    from calder.evaluation.reports import build_report

    try:
        build_report(run_dir)
        typer.echo(f"Report: {(run_dir / 'summary.md').resolve()}")
    except (OSError, ValueError, RuntimeError) as exc:
        _error(exc)


@app.command("validate-contexts")
def validate_contexts(
    k562_config: Path = Path("configs/replogle_k562.yaml"),
    rpe1_config: Path = Path("configs/replogle_rpe1.yaml"),
    output_dir: Path = Path("results/context_validation"),
) -> None:
    """Run K562 then RPE1 with identical frozen method settings."""
    from calder.experiments.transfer import validate_contexts as validate

    try:
        typer.echo(str(validate(load_config(k562_config), load_config(rpe1_config), output_dir)))
    except (OSError, ValueError, RuntimeError) as exc:
        _error(exc)


@app.command()
def demo(run_dir: Annotated[Path, typer.Option()]) -> None:
    """Open optional saved-results explorer; never retrain automatically."""
    if not importlib.util.find_spec("streamlit"):
        _error(RuntimeError("Optional demo needs: uv sync --extra dev --extra demo"))
    script = Path(__file__).with_name("demo.py")
    raise typer.Exit(
        subprocess.call(
            [sys.executable, "-m", "streamlit", "run", str(script), "--", str(run_dir.resolve())]
        )
    )


@app.command("compare-contexts")
def compare_contexts(
    k562_run: Annotated[Path, typer.Option()],
    rpe1_run: Annotated[Path, typer.Option()],
    output_dir: Annotated[Path, typer.Option()],
) -> None:
    """Build context figures from completed runs with frozen hyperparameters."""
    from calder.evaluation.reports import build_context_comparison

    try:
        build_context_comparison(k562_run, rpe1_run, output_dir)
        typer.echo(str(output_dir.resolve()))
    except (OSError, ValueError, RuntimeError) as exc:
        _error(exc)


@app.command()
def doctor(config: Path = Path("configs/synthetic_smoke.yaml")) -> None:
    """Check dependencies, schema, storage and a minimal synthetic replay."""
    required = [
        "numpy",
        "scipy",
        "pandas",
        "sklearn",
        "anndata",
        "scanpy",
        "h5py",
        "pyarrow",
        "networkx",
        "matplotlib",
        "plotly",
        "pydantic",
        "typer",
        "yaml",
    ]
    info: dict = {
        "python": sys.version.split()[0],
        "python_compatible": (3, 11) <= sys.version_info[:2] < (3, 14),
        "dependencies": {name: importlib.util.find_spec(name) is not None for name in required},
        "disk_free_gb": shutil.disk_usage(Path.cwd()).free / 1e9,
        "gpu": "optional; not queried unless torch installed",
    }
    failed = not info["python_compatible"] or not all(info["dependencies"].values())
    try:
        settings = load_config(config)
        output = Path(settings.output_dir)
        output.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=output) as temporary:
            from calder.experiments.runner import run_benchmark

            smoke = ExperimentConfig()
            smoke.synthetic.n_candidates = 40
            smoke.synthetic.scenarios = ["A"]
            smoke.model.n_features = 16
            smoke.protocol.seeds = [0]
            smoke.protocol.threshold_quantiles = [0.8]
            smoke.protocol.policies = ["random", "probability_of_hit", "calder"]
            smoke.protocol.include_oracle = False
            smoke.protocol.budget_fractions = [0.1]
            run_benchmark(smoke, Path(temporary) / "smoke", make_report=False)
        info["smoke_test"] = "passed"
        if settings.dataset.name != "synthetic":
            from calder.data.adapters import load_prepared

            if not settings.dataset.prepared_path:
                raise ValueError("Configured real dataset lacks prepared_path")
            pool, _ = load_prepared(settings.dataset.prepared_path)
            info["prepared_candidates"] = len(pool.gene_ids)
        if importlib.util.find_spec("torch"):
            import torch

            info["gpu"] = bool(torch.cuda.is_available())
        info["output_writable"] = True
        info["disk_guidance"] = (
            "Keep raw H5AD outside Git; download verifies source size; allow room for raw plus scratch data."
        )
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        info["error"] = str(exc)
        failed = True
    typer.echo(json.dumps(info, indent=2))
    raise typer.Exit(1 if failed else 0)


if __name__ == "__main__":
    app()
