"""Line-buffered append-only event audit, with pre-query records written first."""

from pathlib import Path

from calder.utils.serialization import canonical_json


class AuditLog:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def append(self, record: dict) -> None:
        with self.path.open("a") as handle:
            handle.write(canonical_json(record) + "\n")
            handle.flush()
