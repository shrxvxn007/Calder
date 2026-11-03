"""Sparse, control-only representations and retrospective response preparation."""

from calder.data.adapters import load_prepared, prepare_dataset
from calder.data.download import download_dataset, verify_remote
from calder.data.inspect import inspect_dataset

__all__ = [
    "download_dataset",
    "inspect_dataset",
    "load_prepared",
    "prepare_dataset",
    "verify_remote",
]
