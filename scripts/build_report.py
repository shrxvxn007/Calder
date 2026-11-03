"""Build figures and a report from saved outputs: --run-dir results/<run_id>."""

import sys

from calder.cli import app

if __name__ == "__main__":
    app(args=["report", *sys.argv[1:]], prog_name="calder")
