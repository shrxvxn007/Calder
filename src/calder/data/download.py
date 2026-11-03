"""Bounded, provenance-preserving public dataset downloads; never substitute data."""

from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

DATASETS = {
    name: f"https://storage.googleapis.com/scperteval/processed/{name}_processed_complete.h5ad"
    for name in ("replogle22k562", "replogle22rpe1")
}
SCPERTURB = {
    "replogle22k562": {
        "filename": "ReplogleWeissman2022_K562_essential.h5ad",
        "size_bytes": 1546729675,
        "md5_hex": "d8cba17576d1a8afc0f7d71b79cad0f7",
    },
    "replogle22rpe1": {
        "filename": "ReplogleWeissman2022_rpe1.h5ad",
        "size_bytes": 1236886900,
        "md5_hex": "cc7f1ec50aeb3a3e1b4a6cfa713d80fa",
    },
}


def dataset_url(dataset: str, source: str = "scperteval") -> str:
    """Return a declared source URL, rejecting unrecognized identifiers."""
    if dataset not in DATASETS:
        raise ValueError(f"Unknown dataset {dataset!r}; choose from {sorted(DATASETS)}")
    if source == "scperteval":
        return DATASETS[dataset]
    if source == "scperturb":
        return (
            f"https://zenodo.org/records/13350497/files/{SCPERTURB[dataset]['filename']}?download=1"
        )
    raise ValueError(
        "source must be 'scperteval' or 'scperturb'; alternatives are never substituted silently"
    )


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_remote(dataset: str, timeout: int = 30, *, source: str = "scperteval") -> dict:
    url = dataset_url(dataset, source)
    if source == "scperturb":
        expected = SCPERTURB[dataset]
        try:
            with urlopen(Request(url, method="HEAD"), timeout=timeout) as response:
                size = int(response.headers.get("Content-Length", "0"))
                if size != expected["size_bytes"]:
                    raise RuntimeError(
                        "scPerturb remote size differs from the versioned Zenodo release metadata"
                    )
                return {
                    "dataset": dataset,
                    "url": url,
                    "status": response.status,
                    **expected,
                    "source": source,
                    "source_release": "scPerturb 1.4 / Zenodo 13350497",
                    "source_record": "https://zenodo.org/records/13350497",
                    "content_type": response.headers.get("Content-Type"),
                    "normalization_note": "raw counts expected: inspect and validate before processing",
                    "verification_method": "HTTP HEAD against pinned Zenodo v1.4 size and published MD5",
                }
        except (HTTPError, URLError, TimeoutError, OSError, ValueError) as error:
            raise RuntimeError(
                f"Cannot verify explicit scPerturb source {url}: {error}; supply a verified local H5AD"
            ) from error
    try:
        with urlopen(Request(url, method="HEAD"), timeout=timeout) as response:
            length = response.headers.get("Content-Length")
            return {
                "dataset": dataset,
                "source": source,
                "url": url,
                "status": response.status,
                "size_bytes": int(length) if length else None,
                "etag": response.headers.get("ETag"),
                "google_hash": response.headers.get("x-goog-hash"),
                "content_type": response.headers.get("Content-Type"),
                "verification_method": "HTTP HEAD",
            }
    except (HTTPError, URLError, TimeoutError, OSError) as original:
        # GCS object metadata GET is useful when an HTTP proxy disallows HEAD.
        name = quote(f"processed/{dataset}_processed_complete.h5ad", safe="")
        metadata_url = f"https://storage.googleapis.com/storage/v1/b/scperteval/o/{name}"
        try:
            with urlopen(metadata_url, timeout=timeout) as response:
                metadata = json.load(response)
            return {
                "dataset": dataset,
                "source": source,
                "url": url,
                "status": 200,
                "size_bytes": int(metadata["size"]),
                "etag": metadata.get("etag"),
                "md5_base64": metadata.get("md5Hash"),
                "generation": metadata.get("generation"),
                "content_type": metadata.get("contentType"),
                "verification_method": "GCS public object metadata GET",
            }
        except (HTTPError, URLError, TimeoutError, OSError, KeyError, ValueError) as fallback:
            raise RuntimeError(
                f"Cannot verify {dataset} at {url}. HEAD: {original}; metadata: {fallback}. "
                "Check network access/source availability, or supply a local .h5ad with "
                "--path to calder data inspect/prepare. No substitute data were used."
            ) from fallback


def download_dataset(
    dataset: str,
    destination: str | Path,
    max_bytes: int = 2_000_000_000,
    timeout: int = 60,
    *,
    source: str = "scperteval",
) -> Path:
    metadata = verify_remote(dataset, timeout, source=source)
    size = metadata.get("size_bytes")
    if size is None:
        raise RuntimeError(
            "Remote size is unavailable; provide a local H5AD or verify source metadata first"
        )
    if size > max_bytes:
        raise RuntimeError(
            f"{dataset} is {size:,} bytes, exceeding the configured {max_bytes:,}-byte limit. "
            "Increase max_bytes explicitly if disk and bandwidth permit, or supply a local .h5ad."
        )
    import shutil

    destination = Path(destination)
    if destination.suffix != ".h5ad":
        destination = destination / metadata.get("filename", f"{dataset}_processed_complete.h5ad")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(destination.parent).free < size + 100_000_000:
        raise RuntimeError(f"Insufficient disk space for {size:,}-byte download at {destination}")
    if destination.exists():
        raise FileExistsError(
            f"{destination} already exists; inspect it or choose another destination"
        )
    temporary = destination.with_suffix(".h5ad.part")
    try:
        received = 0
        with urlopen(metadata["url"], timeout=timeout) as response, temporary.open("wb") as target:
            while block := response.read(8 * 1024 * 1024):
                received += len(block)
                if received > max_bytes or received > size:
                    raise RuntimeError(
                        "Downloaded body exceeds verified size or configured byte limit"
                    )
                target.write(block)
        if received != size:
            raise RuntimeError(f"Incomplete download: received {received:,} of {size:,} bytes")
        # HDF5 file signature check detects HTML/error payloads before publication.
        with temporary.open("rb") as handle:
            if handle.read(8) != b"\x89HDF\r\n\x1a\n":
                raise RuntimeError("Downloaded object is not an HDF5/H5AD file")
        metadata["sha256"] = file_sha256(temporary)
        remote_md5 = metadata.get("md5_base64")
        if remote_md5 is None and metadata.get("google_hash"):
            for item in metadata["google_hash"].split(","):
                key, _, value = item.strip().partition("=")
                if key == "md5":
                    remote_md5 = value
        if remote_md5 or metadata.get("md5_hex"):
            digest = hashlib.md5(usedforsecurity=False)
            with temporary.open("rb") as handle:
                for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                    digest.update(block)
            actual_md5 = base64.b64encode(digest.digest()).decode("ascii")
            if remote_md5 is not None and actual_md5 != remote_md5:
                raise RuntimeError("Downloaded object fails the source-provided MD5 checksum")
            if metadata.get("md5_hex") is not None and digest.hexdigest() != metadata["md5_hex"]:
                raise RuntimeError("Downloaded object fails the pinned Zenodo source MD5 checksum")
            metadata["source_md5_verified"] = True
        metadata["received_bytes"] = received
        temporary.replace(destination)
        destination.with_suffix(".download.json").write_text(json.dumps(metadata, indent=2) + "\n")
        return destination
    except (HTTPError, URLError, TimeoutError, OSError) as error:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(
            f"Download failed for {metadata['url']}: {error}; supply a local .h5ad"
        ) from error
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
