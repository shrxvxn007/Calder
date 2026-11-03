"""Stable JSON, hashes and environment provenance."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import subprocess
from pathlib import Path
from typing import Any

import numpy as np


def json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value)}")


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), default=json_default, allow_nan=False
    )


def content_hash(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, default=json_default, allow_nan=False) + "\n"
    )


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def environment_provenance() -> dict:
    def git(*args: str) -> str | None:
        result = subprocess.run(["git", *args], capture_output=True, text=True, check=False)
        return result.stdout.strip() if result.returncode == 0 else None

    dependencies = {}
    for package in (
        "numpy",
        "scipy",
        "pandas",
        "scikit-learn",
        "anndata",
        "scanpy",
        "h5py",
        "networkx",
        "matplotlib",
        "pydantic",
        "typer",
        "pyarrow",
    ):
        try:
            dependencies[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            dependencies[package] = "unavailable"
    project_root = Path(__file__).resolve().parents[3]
    source_paths = [
        *sorted((project_root / "src").glob("**/*.py")),
        project_root / "pyproject.toml",
        project_root / "uv.lock",
    ]
    source_files = {
        str(path.relative_to(project_root)): file_sha256(path)
        for path in source_paths
        if path.is_file()
    }
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "git_commit": git("rev-parse", "HEAD"),
        "git_dirty": bool(git("status", "--porcelain")),
        "git_tracked_dirty": bool(git("status", "--porcelain", "--untracked-files=no")),
        "source_tree_hash": content_hash(source_files),
        "source_files": source_files,
        "dependencies": dependencies,
    }
