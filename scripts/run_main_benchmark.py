"""Run K562 by default; pass --config configs/replogle_rpe1.yaml for validation."""

import sys

from calder.cli import app

if __name__ == "__main__":
    args = sys.argv[1:]
    if "--config" not in args:
        args = ["--config", "configs/replogle_k562.yaml", *args]
    app(args=["benchmark", *args], prog_name="calder")
