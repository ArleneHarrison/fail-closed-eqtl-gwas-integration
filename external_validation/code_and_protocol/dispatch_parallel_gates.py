#!/usr/bin/env python3
"""Run the unchanged external gate in independent shards and merge outputs.

This is an orchestration-only wrapper.  Every shard invokes
``run_external_validation.py gate``; the resulting per-unit ``gate.json`` is
copied byte-for-byte into one output tree.  BLAS thread counts are fixed at one
per worker to avoid nested oversubscription.  The wrapper refuses an existing
output directory and validates exact coverage of the supplied attempt registry.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_attempts(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    if not rows or "unit_id" not in fields:
        raise ValueError("attempt registry is empty or lacks unit_id")
    unit_ids = [row["unit_id"] for row in rows]
    if len(unit_ids) != len(set(unit_ids)):
        raise ValueError("attempt registry contains duplicate unit_id values")
    return fields, rows


def write_attempts(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_shard(
    shard_id: int,
    shard_attempts: Path,
    shard_out: Path,
    args: argparse.Namespace,
) -> dict[str, object]:
    command = [
        args.python,
        str(args.gate_script),
        "gate",
        "--protocol", str(args.protocol),
        "--metadata-lock", str(args.metadata_lock),
        "--attempts", str(shard_attempts),
        "--prepared", str(args.prepared),
        "--out", str(shard_out),
    ]
    environment = os.environ.copy()
    environment.update({
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    })
    completed = subprocess.run(command, text=True, capture_output=True, env=environment)
    return {
        "shard_id": shard_id,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "command": command,
    }


def main(args: argparse.Namespace) -> None:
    fields, attempts = read_attempts(args.attempts)
    jobs = min(max(1, args.jobs), len(attempts))
    if args.out.exists():
        raise FileExistsError(f"refusing existing output path: {args.out}")
    for required in (args.gate_script, args.protocol, args.metadata_lock, args.attempts):
        if not required.is_file():
            raise FileNotFoundError(required)
    if not args.prepared.is_dir():
        raise FileNotFoundError(args.prepared)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="parallel_gate_shards_", dir=args.out.parent) as temporary:
        temp_root = Path(temporary)
        shards: list[list[dict[str, str]]] = [[] for _ in range(jobs)]
        for index, row in enumerate(attempts):
            shards[index % jobs].append(row)

        futures = []
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            for shard_id, rows in enumerate(shards):
                shard_attempts = temp_root / f"attempts_{shard_id:03d}.tsv"
                shard_out = temp_root / f"gate_{shard_id:03d}"
                write_attempts(shard_attempts, fields, rows)
                futures.append(executor.submit(run_shard, shard_id, shard_attempts, shard_out, args))
            results = [future.result() for future in as_completed(futures)]

        failed = [result for result in results if result["returncode"] != 0]
        if failed:
            details = "\n".join(
                f"shard {item['shard_id']} rc={item['returncode']} stderr={item['stderr']}"
                for item in sorted(failed, key=lambda value: int(value["shard_id"]))
            )
            raise RuntimeError(f"one or more gate shards failed:\n{details}")

        args.out.mkdir(parents=False)
        units_out = args.out / "units"
        units_out.mkdir()
        trace: list[dict[str, str]] = []
        for index, attempt in enumerate(attempts):
            unit_id = attempt["unit_id"]
            shard_id = index % jobs
            source = temp_root / f"gate_{shard_id:03d}" / "units" / unit_id / "gate.json"
            if not source.is_file():
                raise RuntimeError(f"missing gate output for {unit_id}")
            destination = units_out / unit_id / "gate.json"
            destination.parent.mkdir()
            shutil.copy2(source, destination)
            gate = json.loads(destination.read_text(encoding="utf-8"))
            trace.append({
                "unit_id": unit_id,
                "stage": "gate",
                "status": str(gate.get("status", "")),
                "status_codes": ";".join(gate.get("status_codes", [])),
                "gate_sha256": sha256_file(destination),
            })

        observed = {path.parent.name for path in units_out.glob("*/gate.json")}
        expected = {row["unit_id"] for row in attempts}
        if observed != expected:
            raise RuntimeError("merged gate output does not exactly cover the attempt registry")

        with (args.out / "trace.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(trace[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(trace)
        first_shard = temp_root / "gate_000"
        shutil.copy2(first_shard / "environment.txt", args.out / "environment.txt")
        shutil.copy2(first_shard / "run_configuration.tsv", args.out / "run_configuration.tsv")
        audit = {
            "schema_version": "parallel_gate_dispatch_1.0",
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "attempts": len(attempts),
            "jobs": jobs,
            "gate_script_sha256": sha256_file(args.gate_script),
            "attempts_sha256": sha256_file(args.attempts),
            "protocol_sha256_file": sha256_file(args.protocol),
            "metadata_lock_sha256": sha256_file(args.metadata_lock),
            "per_worker_blas_threads": 1,
            "merge_policy": "byte_copy_per_unit_gate_json_in_original_attempt_order",
            "status": "PASS",
        }
        (args.out / "parallel_dispatch_audit.json").write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(json.dumps({"status": "PASS", "attempts": len(attempts), "jobs": jobs, "out": str(args.out)}))


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gate-script", type=Path, required=True)
    p.add_argument("--protocol", type=Path, required=True)
    p.add_argument("--metadata-lock", type=Path, required=True)
    p.add_argument("--attempts", type=Path, required=True)
    p.add_argument("--prepared", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--python", default=sys.executable)
    return p


if __name__ == "__main__":
    try:
        main(parser().parse_args())
    except Exception as error:
        print(f"PARALLEL_GATE_STOPPED: {error}", file=sys.stderr)
        raise SystemExit(2)
