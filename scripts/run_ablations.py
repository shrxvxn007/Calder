"""Execute predeclared ablation variants via the CLI."""

import sys

from calder.cli import app

if __name__ == "__main__":
    args = sys.argv[1:]
    if "--config" not in args:
        args = ["--config", "configs/ablations.yaml", *args]
    app(args=["ablate", *args], prog_name="calder")
