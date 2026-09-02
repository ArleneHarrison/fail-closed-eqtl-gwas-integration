#!/usr/bin/env python3
"""Parallel, resumable, READY-only dispatcher for the frozen R model runner.

The dispatcher consumes attempt and gate registries, refuses every non-READY
unit, binds each model run to immutable inputs and the model implementation
lock, and publishes unit status atomically.  It performs no interpretation of
posterior results.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any

import numpy as np


REQUIRED_RESULTS = (
    "status.csv", "model_implementation_policy.csv", "ld_model_policy.csv",
    "summary.csv", "pairwise_results.csv", "coloc_susie_result.rds", "R_sessionInfo.txt",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(payload).hexdigest()


def read_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def protocol_sha(path: Path) -> str:
    return canonical_sha(read_json(path))


def verify_implementation_lock(args: argparse.Namespace) -> dict[str, Any]:
    lock, freeze = read_json(args.implementation_lock), read_json(args.implementation_freeze)
    expected_protocol_sha = protocol_sha(args.protocol)
    expected_runner_sha = sha256_file(args.runner)
    checks = {
        "lock_state": lock.get("state") == "LOCKED_BEFORE_FORMAL_MODEL_RESULTS",
        "freeze_state": freeze.get("state") == "FROZEN",
        "protocol_lock": lock.get("protocol_canonical_json_sha256") == expected_protocol_sha,
        "protocol_freeze": freeze.get("protocol_canonical_json_sha256") == expected_protocol_sha,
        "runner_lock": lock.get("runner_sha256") == expected_runner_sha,
        "runner_freeze": freeze.get("runner_sha256") == expected_runner_sha,
        "lock_hash": freeze.get("lock_sha256") == sha256_file(args.implementation_lock),
        "coloc_version": lock.get("coloc_version") == "5.2.3",
        "susieR_version": lock.get("susieR_version") == "0.14.2",
        "susie_L": lock.get("susie_L") == 10,
        "seed": lock.get("seed") == 20260901,
        "max_iter": lock.get("max_iter") == 1000,
        "repeat_until_convergence": lock.get("repeat_until_convergence") is False,
        "estimate_residual_variance": lock.get("estimate_residual_variance") is False,
        "tol": lock.get("tol") == 0.001,
        "coverage": lock.get("coverage") == 0.95,
        "min_abs_corr": lock.get("min_abs_corr") == 0.5,
        "scaled_prior_variance": lock.get("scaled_prior_variance") == 0.2,
        "check_prior": lock.get("check_prior") is True,
        "z_ld_weight": lock.get("z_ld_weight") == 0,
        "coloc_p1": lock.get("coloc_p1") == 1e-4,
        "coloc_p2": lock.get("coloc_p2") == 1e-4,
        "coloc_p12": lock.get("coloc_p12") == 5e-6,
        "coloc_overlap_min": lock.get("coloc_overlap_min") == 0.5,
        "coloc_trim_by_posterior": lock.get("coloc_trim_by_posterior") is True,
        "coloc_back_calculate_lbf": lock.get("coloc_back_calculate_lbf") is False,
        "prior_variance_api_default": lock.get("prior_variance_api_default") == 50,
        "prior_variance_activity": lock.get("prior_variance_activity") == "inactive_when_N_is_supplied_in_susieR_0.14.2",
        "back_calculate_lbf_activity": lock.get("back_calculate_lbf_activity") == "formal_accepted_but_unreferenced_in_coloc_5.2.3_body",
        "prior_engineering_smoke_model_runs_seen": lock.get("prior_engineering_smoke_model_runs_seen") is True,
        "formal_primary_model_results_seen_before_lock": lock.get("formal_primary_model_results_seen_before_lock") is False,
        "prior_smoke_outputs_excluded_from_formal_primary": lock.get("prior_smoke_outputs_excluded_from_formal_primary") is True,
    }
    failures = [key for key, passed in checks.items() if not passed]
    if failures:
        raise RuntimeError("E_IMPLEMENTATION_LOCK_INVALID:" + ",".join(failures))
    return {
        "protocol_canonical_json_sha256": expected_protocol_sha,
        "runner_sha256": expected_runner_sha,
        "implementation_lock_sha256": sha256_file(args.implementation_lock),
        "model_parameters": {
            "coloc_version": lock["coloc_version"], "susieR_version": lock["susieR_version"],
            "susie_L": lock["susie_L"], "seed": lock["seed"], "max_iter": lock["max_iter"],
            "repeat_until_convergence": lock["repeat_until_convergence"],
            "estimate_residual_variance": lock["estimate_residual_variance"],
            "tol": lock["tol"], "coverage": lock["coverage"],
            "min_abs_corr": lock["min_abs_corr"],
            "scaled_prior_variance": lock["scaled_prior_variance"],
            "check_prior": lock["check_prior"], "z_ld_weight": lock["z_ld_weight"],
            "coloc_p1": lock["coloc_p1"], "coloc_p2": lock["coloc_p2"], "coloc_p12": lock["coloc_p12"],
            "coloc_overlap_min": lock["coloc_overlap_min"],
            "coloc_trim_by_posterior": lock["coloc_trim_by_posterior"],
            "coloc_back_calculate_lbf": lock["coloc_back_calculate_lbf"],
            "prior_variance_api_default": lock["prior_variance_api_default"],
            "prior_variance_activity": lock["prior_variance_activity"],
            "back_calculate_lbf_activity": lock["back_calculate_lbf_activity"],
        },
    }


def parse_time(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    result: dict[str, Any] = {}
    for label, key, cast in (
        ("User time (seconds)", "user_seconds", str),
        ("System time (seconds)", "system_seconds", str),
        ("Elapsed (wall clock) time (h:mm:ss or m:ss)", "elapsed_clock", str),
        ("Maximum resident set size (kbytes)", "max_rss_kb", int),
        ("Exit status", "time_exit_status", int),
    ):
        match = re.search(rf"^\s*{re.escape(label)}:\s*(.+?)\s*$", text, re.MULTILINE)
        if match:
            result[key] = cast(match.group(1))
    return result


def materialise_locked_input(prepared: Path, destination: Path) -> dict[str, str]:
    destination.mkdir(parents=True, exist_ok=False)
    shutil.copy2(prepared / "summary.csv", destination / "summary.csv")
    archive = np.load(prepared / "ld.npz", allow_pickle=False)
    ld = archive["ld"]
    identifiers = [str(value) for value in archive["variant_ids"].tolist()]
    with (destination / "ld.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([""] + identifiers)
        for identifier, values in zip(identifiers, ld, strict=True):
            writer.writerow([identifier] + [format(float(value), ".17g") for value in values])
    return {
        "summary_csv_sha256": sha256_file(destination / "summary.csv"),
        "ld_csv_sha256": sha256_file(destination / "ld.csv"),
    }


def read_r_status(path: Path, returncode: int) -> tuple[str, str, str]:
    if returncode in {124, 137}:
        return "FAILED", "E_MODEL_TIMEOUT", "bounded model runtime exceeded"
    if not path.is_file():
        return "FAILED", "E_MODEL_STATUS_MISSING", "R runner did not create status.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle), {})
    r_status = row.get("status", "")
    code = row.get("stable_code", "") or "E_MODEL_RUNNER"
    reason = row.get("reason", "")
    if returncode == 0 and r_status in {"completed", "completed_no_comparable_credible_sets"}:
        return "COMPLETED", code, reason
    return "FAILED", code, reason or f"R runner return code {returncode}"


def valid_completed_resume(status: dict[str, Any], binding_sha: str, unit_dir: Path) -> bool:
    if status.get("status") != "COMPLETED" or status.get("binding_sha256") != binding_sha:
        return False
    run_dir_name = status.get("run_directory")
    if not run_dir_name:
        return False
    run_dir = unit_dir / run_dir_name
    hashes = status.get("result_sha256", {})
    return bool(hashes) and all((run_dir / name).is_file() and sha256_file(run_dir / name) == digest
                                for name, digest in hashes.items())


def rejection(unit_dir: Path, unit_id: str, code: str, reason: str, gate_sha: str | None = None) -> dict[str, Any]:
    status = {
        "unit_id": unit_id, "status": "REJECTED", "stable_code": code, "reason": reason,
        "gate_sha256": gate_sha, "model_invoked": False,
    }
    atomic_json(unit_dir / "status.json", status)
    return status


def dispatch_unit(
    attempt: dict[str, str], args: argparse.Namespace, lock_binding: dict[str, str],
    metadata: dict[str, Any],
) -> dict[str, Any]:
    unit_id = attempt["unit_id"]
    unit_dir = args.out / "units" / unit_id
    unit_dir.mkdir(parents=True, exist_ok=True)
    gate_path = args.gates / "units" / unit_id / "gate.json"
    if not gate_path.is_file():
        return rejection(unit_dir, unit_id, "E_GATE_MISSING", "gate JSON is absent")
    gate_sha, gate = sha256_file(gate_path), read_json(gate_path)
    if gate.get("status") != "READY" or gate.get("gate_ready_for_dispatch") is not True:
        return rejection(
            unit_dir, unit_id, "E_GATE_NOT_READY",
            "non-READY units are refused and cannot create a posterior artifact", gate_sha,
        )
    prepared = args.prepared / unit_id
    if not prepared.is_dir():
        return rejection(unit_dir, unit_id, "E_PREPARED_UNIT_MISSING", "prepared unit directory is absent", gate_sha)
    required_inputs = ("summary.csv", "ld.npz", "provenance.json", "overlap.json")
    if not all((prepared / name).is_file() for name in required_inputs):
        return rejection(unit_dir, unit_id, "E_PREPARED_INPUT_MISSING", "one or more prepared inputs are absent", gate_sha)
    source_hashes = {name: sha256_file(prepared / name) for name in required_inputs}
    if gate.get("prepared_input_sha256") and gate["prepared_input_sha256"] != source_hashes:
        return rejection(unit_dir, unit_id, "E_GATE_INPUT_BINDING_MISMATCH", "prepared input hashes differ from gate JSON", gate_sha)
    outcome = attempt.get("outcome_id", "")
    if outcome not in metadata.get("outcomes", {}):
        return rejection(unit_dir, unit_id, "E_OUTCOME_METADATA_MISSING", "outcome is absent from metadata lock", gate_sha)
    case_fraction = float(metadata["outcomes"][outcome]["case_fraction"])
    if not 0 < case_fraction < 1:
        return rejection(unit_dir, unit_id, "E_CASE_FRACTION_INVALID", "locked case fraction is invalid", gate_sha)
    binding = {
        "unit_id": unit_id, "outcome_id": outcome, "case_fraction": case_fraction,
        "gate_sha256": gate_sha, "gate_ready_for_dispatch": True,
        "prepared_input_sha256": source_hashes, **lock_binding,
    }
    binding_sha = canonical_sha(binding)
    status_path = unit_dir / "status.json"
    if status_path.is_file():
        previous = read_json(status_path)
        if valid_completed_resume(previous, binding_sha, unit_dir):
            return {**previous, "resume_action": "SKIPPED_VERIFIED_COMPLETED"}
        if (previous.get("binding_sha256") == binding_sha and previous.get("status") == "FAILED"
                and not args.retry_failed):
            return {**previous, "resume_action": "SKIPPED_PREVIOUS_FAILED_USE_RETRY_FAILED_TO_RERUN"}

    staging = Path(tempfile.mkdtemp(prefix=".staging_", dir=unit_dir))
    try:
        locked = staging / "locked_input"
        derived_hashes = materialise_locked_input(prepared, locked)
        atomic_json(locked / "binding.json", {**binding, **derived_hashes, "binding_sha256": binding_sha})
        stdout_path, stderr_path = staging / "stdout.log", staging / "stderr.log"
        time_path = staging / "gnu_time.txt"
        command = [
            "/usr/bin/time", "-v", "-o", str(time_path),
            "timeout", "--signal=TERM", "--kill-after=30s", str(args.timeout_seconds),
            args.rscript, str(args.runner), str(locked), str(case_fraction), str(staging),
        ]
        started = time.monotonic()
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        wall_seconds = round(time.monotonic() - started, 6)
        stdout_path.write_text(completed.stdout, encoding="utf-8")
        stderr_path.write_text(completed.stderr, encoding="utf-8")
        status, code, reason = read_r_status(staging / "status.csv", completed.returncode)
        missing_results = [name for name in REQUIRED_RESULTS if not (staging / name).is_file()]
        if status == "COMPLETED" and missing_results:
            status, code, reason = "FAILED", "E_MODEL_RESULT_INCOMPLETE", ",".join(missing_results)
        result_hashes = {
            name: sha256_file(staging / name) for name in REQUIRED_RESULTS if (staging / name).is_file()
        }
        run_name_base = f"run_{binding_sha[:16]}"
        run_name = run_name_base
        index = 1
        while (unit_dir / run_name).exists():
            index += 1
            run_name = f"{run_name_base}_retry{index}"
        final_run = unit_dir / run_name
        os.replace(staging, final_run)
        published = {
            "unit_id": unit_id, "status": status, "stable_code": code, "reason": reason,
            "model_invoked": True, "binding_sha256": binding_sha, "binding": binding,
            "run_directory": run_name, "returncode": completed.returncode,
            "elapsed_seconds": wall_seconds, "result_sha256": result_hashes,
            **parse_time(final_run / "gnu_time.txt"),
        }
        atomic_json(status_path, published)
        return published
    except Exception as error:
        if staging.exists():
            failed_name = f"failed_staging_{int(time.time())}_{os.getpid()}"
            os.replace(staging, unit_dir / failed_name)
        failed = {
            "unit_id": unit_id, "status": "FAILED", "stable_code": "E_DISPATCH_EXCEPTION",
            "reason": str(error), "model_invoked": False, "binding_sha256": binding_sha,
        }
        atomic_json(status_path, failed)
        return failed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--attempts", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--gates", type=Path, required=True)
    parser.add_argument("--metadata-lock", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--implementation-lock", type=Path, required=True)
    parser.add_argument("--implementation-freeze", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=12)
    parser.add_argument("--timeout-seconds", type=int, default=3600)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--retry-failed", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1 or args.timeout_seconds < 1:
        parser.error("jobs and timeout-seconds must be positive")
    attempts = list(csv.DictReader(args.attempts.open(newline="", encoding="utf-8"), delimiter="\t"))
    unit_ids = [attempt.get("unit_id", "") for attempt in attempts]
    if not unit_ids or any(not unit_id for unit_id in unit_ids) or len(unit_ids) != len(set(unit_ids)):
        raise SystemExit("E_ATTEMPT_REGISTRY_INVALID")
    metadata = read_json(args.metadata_lock)
    lock_binding = verify_implementation_lock(args)
    if metadata.get("protocol_canonical_json_sha256") != lock_binding["protocol_canonical_json_sha256"]:
        raise SystemExit("E_METADATA_PROTOCOL_MISMATCH")
    args.out.mkdir(parents=True, exist_ok=True)
    configuration = {
        "attempts_sha256": sha256_file(args.attempts), "metadata_lock_sha256": sha256_file(args.metadata_lock),
        "jobs": args.jobs, "timeout_seconds": args.timeout_seconds, "retry_failed": args.retry_failed,
        "attempt_count": len(attempts), **lock_binding,
    }
    atomic_json(args.out / "run_configuration.json", configuration)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda attempt: dispatch_unit(attempt, args, lock_binding, metadata), attempts))
    atomic_json(args.out / "dispatch_summary.json", {"configuration": configuration, "units": results})
    fields = (
        "unit_id", "status", "stable_code", "model_invoked", "elapsed_seconds", "max_rss_kb",
        "binding_sha256", "run_directory", "resume_action", "reason",
    )
    temporary = args.out / f"dispatch_summary.tsv.tmp.{os.getpid()}"
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            writer.writerow({field: result.get(field) for field in fields})
    os.replace(temporary, args.out / "dispatch_summary.tsv")
    if any(result["status"] == "FAILED" for result in results):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
