"""Download an immutable public input to a new file and admit it only after checks.

The destination stays suffixed ``.inprogress`` until byte count, SHA-256 and
optional gzip-stream validation all pass.  It never resumes or reads an older
partial file, because partial state is not provenance.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


def stream_gzip_test(path: Path) -> None:
    with gzip.open(path, "rb") as source:
        while source.read(1024 * 1024):
            pass


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--destination", required=True)
    parser.add_argument("--expected-bytes", required=True, type=int)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--version-or-build", required=True)
    parser.add_argument("--license-or-access", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--gzip-test", action="store_true")
    args = parser.parse_args()

    destination = Path(args.destination)
    temporary = destination.with_name(destination.name + ".inprogress")
    log = Path(args.log)
    destination.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or temporary.exists():
        raise SystemExit("refusing to overwrite an existing complete or in-progress destination")

    record = {
        "schema_version": "1.0", "url": args.url, "version_or_build": args.version_or_build,
        "license_or_access": args.license_or_access, "retrieval_started_utc": datetime.now(timezone.utc).isoformat(),
        "expected_byte_count": args.expected_bytes, "expected_sha256": args.expected_sha256.lower(),
        "destination": str(destination), "status": "in_progress",
    }
    log.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    digest = hashlib.sha256()
    byte_count = 0
    try:
        with urllib.request.urlopen(args.url, timeout=120) as response, temporary.open("xb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
                digest.update(chunk)
                byte_count += len(chunk)
                if byte_count % (100 * 1024 * 1024) < len(chunk):
                    print(f"downloaded_bytes={byte_count}", flush=True)
        observed_sha = digest.hexdigest()
        if byte_count != args.expected_bytes:
            raise ValueError(f"byte count mismatch: observed {byte_count}; expected {args.expected_bytes}")
        if observed_sha != args.expected_sha256.lower():
            raise ValueError(f"sha256 mismatch: observed {observed_sha}; expected {args.expected_sha256.lower()}")
        if args.gzip_test:
            stream_gzip_test(temporary)
        os.replace(temporary, destination)
        record.update({
            "retrieval_finished_utc": datetime.now(timezone.utc).isoformat(), "byte_count": byte_count,
            "sha256": observed_sha, "structural_check": "gzip_stream_passed" if args.gzip_test else "not_requested",
            "status": "verified",
        })
    except Exception as error:
        record.update({"retrieval_finished_utc": datetime.now(timezone.utc).isoformat(), "status": "failed", "error": str(error)})
        log.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        raise
    log.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": record["status"], "byte_count": byte_count, "sha256": observed_sha}, sort_keys=True))


if __name__ == "__main__":
    main()
