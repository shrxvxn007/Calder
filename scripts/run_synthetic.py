"""Run the reproducible synthetic benchmark through the installed CLI."""

import sys

from calder.cli import app

if __name__ == "__main__":
    args = sys.argv[1:]
    if "--config" not in args:
        args = ["--config", "configs/synthetic.yaml", *args]
    app(args=["synthetic", *args], prog_name="calder")
