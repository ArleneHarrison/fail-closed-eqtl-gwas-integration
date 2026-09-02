#!/usr/bin/env python3
"""Create or verify the immutable downstream model implementation lock."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPECTED = {
    "coloc_version": "5.2.3",
    "susieR_version": "0.14.2",
    "susie_L": 10,
    "seed": 20260901,
    "max_iter": 1000,
    "repeat_until_convergence": False,
    "estimate_residual_variance": False,
    "tol": 0.001,
    "coverage": 0.95,
    "min_abs_corr": 0.5,
    "scaled_prior_variance": 0.2,
    "check_prior": True,
    "z_ld_weight": 0,
    "coloc_p1": 1e-4,
    "coloc_p2": 1e-4,
    "coloc_p12": 5e-6,
    "coloc_overlap_min": 0.5,
    "coloc_trim_by_posterior": True,
    "coloc_back_calculate_lbf": False,
    "prior_variance_api_default": 50,
    "prior_variance_activity": "inactive_when_N_is_supplied_in_susieR_0.14.2",
    "back_calculate_lbf_activity": "formal_accepted_but_unreferenced_in_coloc_5.2.3_body",
    "prior_engineering_smoke_model_runs_seen": True,
    "formal_primary_model_results_seen_before_lock": False,
    "prior_smoke_outputs_excluded_from_formal_primary": True,
}


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_protocol_sha(path: Path) -> str:
    value = json.loads(path.read_text(encoding="utf-8"))
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    return sha256_bytes(encoded)


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def create(args: argparse.Namespace) -> None:
    if args.lock.exists() or args.freeze.exists():
        raise SystemExit("E_IMPLEMENTATION_LOCK_EXISTS: refusing to overwrite an existing lock")
    runner_sha = sha256_file(args.runner)
    protocol_sha = canonical_protocol_sha(args.protocol)
    lock = {
        "schema_version": "1.0",
        "state": "LOCKED_BEFORE_FORMAL_MODEL_RESULTS",
        "locked_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "implementation_correction": (
            "Final pre-primary-result conformance correction: make every result-affecting installed default "
            "reachable through the coloc 5.2.3/susieR 0.14.2 API explicit and prevent coloc::runsusie "
            "from increasing the prospectively declared 1000-iteration limit."
        ),
        "parameter_selection_basis": (
            "Conformance to protocol v2 max_iter and explicit locking of installed package defaults; "
            "not selected from posterior values or biological results."
        ),
        "scientific_design_changed": False,
        "protocol_canonical_json_sha256": protocol_sha,
        "runner_sha256": runner_sha,
        "runner_filename": args.runner.name,
        "model": "SuSiE-RSS_then_coloc.susie",
        **EXPECTED,
    }
    atomic_json(args.lock, lock)
    freeze = {
        "schema_version": "1.0",
        "state": "FROZEN",
        "lock_filename": args.lock.name,
        "lock_sha256": sha256_file(args.lock),
        "protocol_canonical_json_sha256": protocol_sha,
        "runner_sha256": runner_sha,
        "write_policy": "read_only_after_creation",
    }
    atomic_json(args.freeze, freeze)
    os.chmod(args.lock, 0o444)
    os.chmod(args.freeze, 0o444)
    print(json.dumps({"status": "CREATED", "lock_sha256": freeze["lock_sha256"], "runner_sha256": runner_sha}))


def check(args: argparse.Namespace) -> None:
    if not args.lock.is_file() or not args.freeze.is_file():
        raise SystemExit("E_IMPLEMENTATION_LOCK_MISSING")
    lock = json.loads(args.lock.read_text(encoding="utf-8"))
    freeze = json.loads(args.freeze.read_text(encoding="utf-8"))
    failures: list[str] = []
    protocol_sha = canonical_protocol_sha(args.protocol)
    runner_sha = sha256_file(args.runner)
    if lock.get("protocol_canonical_json_sha256") != protocol_sha:
        failures.append("E_IMPLEMENTATION_PROTOCOL_HASH_MISMATCH")
    if lock.get("runner_sha256") != runner_sha:
        failures.append("E_IMPLEMENTATION_RUNNER_HASH_MISMATCH")
    if freeze.get("lock_sha256") != sha256_file(args.lock):
        failures.append("E_IMPLEMENTATION_LOCK_HASH_MISMATCH")
    if freeze.get("protocol_canonical_json_sha256") != protocol_sha:
        failures.append("E_IMPLEMENTATION_FREEZE_PROTOCOL_HASH_MISMATCH")
    if freeze.get("runner_sha256") != runner_sha:
        failures.append("E_IMPLEMENTATION_FREEZE_RUNNER_HASH_MISMATCH")
    for key, expected in EXPECTED.items():
        if lock.get(key) != expected:
            failures.append(f"E_IMPLEMENTATION_PARAMETER_MISMATCH:{key}")
    if lock.get("state") != "LOCKED_BEFORE_FORMAL_MODEL_RESULTS" or freeze.get("state") != "FROZEN":
        failures.append("E_IMPLEMENTATION_STATE_INVALID")
    if failures:
        raise SystemExit(";".join(failures))
    print(json.dumps({
        "status": "PASS",
        "lock_sha256": sha256_file(args.lock),
        "runner_sha256": runner_sha,
        "protocol_canonical_json_sha256": protocol_sha,
    }, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name, function in (("create", create), ("check", check)):
        command = subparsers.add_parser(name)
        command.add_argument("--protocol", type=Path, required=True)
        command.add_argument("--runner", type=Path, required=True)
        command.add_argument("--lock", type=Path, required=True)
        command.add_argument("--freeze", type=Path, required=True)
        command.set_defaults(function=function)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
