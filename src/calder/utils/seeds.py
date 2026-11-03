"""Independent deterministic substreams without Python's randomized hash()."""

import hashlib


def derived_seed(seed: int, name: str) -> int:
    return int.from_bytes(hashlib.sha256(f"{seed}:{name}".encode()).digest()[:4], "little")
