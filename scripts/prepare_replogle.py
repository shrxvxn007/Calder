"""Prepare a verified local/downloaded Replogle AnnData file, K562 by default."""

import sys

from calder.cli import app

if __name__ == "__main__":
    args = sys.argv[1:]
    if "--dataset" not in args:
        args = ["--dataset", "replogle22k562", *args]
    app(args=["data", "prepare", *args], prog_name="calder")
