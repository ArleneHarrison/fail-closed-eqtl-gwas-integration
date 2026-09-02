#!/usr/bin/env python3
"""Audit and benchmark the pinned R downstream runner on prepared smoke units.

This is deliberately separate from raw preparation and the production dispatcher.
It never edits prepared inputs.  Each run receives a freshly materialised locked
input directory, a prospectively locked study-level case fraction, and a bounded
runtime.  Results include stable machine-readable status codes, GNU time resource
metrics, input audits, and byte-level repeat comparisons.
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
import time
from pathlib import Path
from typing import Any

import numpy as np


RESULT_FILES = (
    "status.csv",
    "ld_model_policy.csv",
    "summary.csv",
    "pairwise_results.csv",
    "coloc_susie_result.rds",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def parse_time_metrics(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    labels = {
        "User time (seconds)": "user_seconds",
        "System time (seconds)": "system_seconds",
        "Elapsed (wall clock) time (h:mm:ss or m:ss)": "elapsed_clock",
        "Maximum resident set size (kbytes)": "peak_rss_kb",
        "Exit status": "time_exit_status",
    }
    result: dict[str, Any] = {}
    for source, target in labels.items():
        match = re.search(rf"^\s*{re.escape(source)}:\s*(.+?)\s*$", text, re.MULTILINE)
        if match:
            raw = match.group(1)
            result[target] = int(raw) if target in {"peak_rss_kb", "time_exit_status"} else raw
    return result


def classify_run(returncode: int, status_file: Path, timed_out: bool) -> tuple[str, str, str]:
    if timed_out or returncode == 124:
        return "FAILED", "E_MODEL_TIMEOUT", "bounded benchmark runtime exceeded"
    if not status_file.is_file():
        return "FAILED", "E_MODEL_STATUS_MISSING", "runner did not create status.csv"
    with status_file.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle), {})
    status, reason = row.get("status", ""), row.get("reason", "")
    if returncode == 0 and status == "completed":
        return "COMPLETED", "OK", reason
    if returncode == 0 and status == "completed_no_comparable_credible_sets":
        return "COMPLETED", "OK_NO_COMPARABLE_CREDIBLE_SETS", reason
    if "version-pinned coloc and susieR packages are required" in reason:
        return "FAILED", "E_R_PACKAGES_MISSING", reason
    if "orders are not identical" in reason:
        return "FAILED", "E_LD_VARIANT_ORDER_MISMATCH", reason
    if "not numerically PSD" in reason:
        return "FAILED", "E_LD_NOT_PSD", reason
    if "LD matrix is asymmetric" in reason:
        return "FAILED", "E_LD_ASYMMETRIC", reason
    if "LD matrix diagonal is not one" in reason:
        return "FAILED", "E_LD_DIAGONAL_INVALID", reason
    return "FAILED", "E_MODEL_RUNNER", reason or f"runner return code {returncode}"


def audit_unit(unit: Path, gate_path: Path, outcome_lock: dict[str, Any]) -> dict[str, Any]:
    gate = read_json(gate_path)
    archive = np.load(unit / "ld.npz", allow_pickle=False)
    ld = archive["ld"]
    variant_ids = [str(value) for value in archive["variant_ids"].tolist()]
    with (unit / "summary.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0]) if rows else []
    required = {
        "target_variant", "position", "eqtl_beta", "eqtl_se", "gwas_beta", "gwas_se",
        "eqtl_maf", "gwas_maf", "eqtl_n", "gwas_n", "gwas_cases",
        "alignment_status", "vcf_status", "outcome_id",
    }
    missing = sorted(required.difference(columns))
    row_ids = [row.get("target_variant", "") for row in rows]
    outcome_ids = sorted({row.get("outcome_id", "") for row in rows})
    failures: list[str] = []
    if gate.get("status") != "READY":
        failures.append("E_GATE_NOT_READY")
    if missing:
        failures.append("E_REQUIRED_COLUMN_MISSING")
    if row_ids != variant_ids:
        failures.append("E_LD_VARIANT_ORDER_MISMATCH")
    if ld.shape != (len(rows), len(rows)):
        failures.append("E_LD_DIMENSION_MISMATCH")
    if len(outcome_ids) != 1:
        failures.append("E_OUTCOME_ID_HETEROGENEOUS")
    accepted = {"aligned", "flipped"}
    alignment_values = sorted({row.get("alignment_status", "").lower() for row in rows})
    if not set(alignment_values).issubset(accepted):
        failures.append("E_ALLELE_UNRESOLVED")

    def numeric(name: str) -> np.ndarray:
        try:
            return np.asarray([float(row[name]) for row in rows], dtype=float)
        except (KeyError, TypeError, ValueError):
            return np.asarray([], dtype=float)

    eqtl_maf, gwas_maf = numeric("eqtl_maf"), numeric("gwas_maf")
    if (eqtl_maf.size != len(rows) or gwas_maf.size != len(rows)
            or not np.isfinite(eqtl_maf).all() or not np.isfinite(gwas_maf).all()
            or np.any((eqtl_maf <= 0) | (eqtl_maf > 0.5))
            or np.any((gwas_maf <= 0) | (gwas_maf > 0.5))):
        failures.append("E_MAF_INVALID")
    maf_columns_separate = bool(
        eqtl_maf.size == len(rows) and gwas_maf.size == len(rows)
        and not np.array_equal(eqtl_maf, gwas_maf)
    )
    if not maf_columns_separate:
        failures.append("E_MAF_COLUMNS_NOT_DISTINCT")
    finite_effects = all(
        numeric(name).size == len(rows) and np.isfinite(numeric(name)).all()
        for name in ("eqtl_beta", "eqtl_se", "gwas_beta", "gwas_se")
    )
    positive_se = all(np.all(numeric(name) > 0) for name in ("eqtl_se", "gwas_se"))
    if not finite_effects or not positive_se:
        failures.append("E_EFFECT_OR_SE_INVALID")

    gwas_n, gwas_cases = numeric("gwas_n"), numeric("gwas_cases")
    row_fractions = gwas_cases / gwas_n if len(rows) else np.asarray([])
    locked_fraction = float(outcome_lock["case_fraction"])
    locked_n = int(outcome_lock["sample_size"])
    runner_gwas_n = float(np.median(gwas_n)) if gwas_n.size else None
    runner_eqtl_n = float(np.median(numeric("eqtl_n"))) if len(rows) else None
    return {
        "unit_id": unit.name,
        "gate_status": gate.get("status"),
        "gate_codes": gate.get("status_codes", []),
        "audit_status": "PASS" if not failures else "FAIL",
        "audit_codes": failures or ["OK"],
        "n_variants": len(rows),
        "ld_n": int(ld.shape[0]),
        "ld_rank_from_gate": gate.get("ld_diagnostics", {}).get("numerical_rank"),
        "summary_ld_order_identical": row_ids == variant_ids,
        "effect_orientation_contract": "both betas encoded to LD ALT; source orientation represented as aligned/flipped",
        "alignment_status_values": alignment_values,
        "finite_effects_and_positive_se": finite_effects and positive_se,
        "eqtl_maf_gwas_maf_separate": maf_columns_separate,
        "eqtl_maf_range": [float(eqtl_maf.min()), float(eqtl_maf.max())] if eqtl_maf.size else None,
        "gwas_maf_range": [float(gwas_maf.min()), float(gwas_maf.max())] if gwas_maf.size else None,
        "case_fraction_policy": "study_level",
        "locked_case_fraction": locked_fraction,
        "row_case_fraction_min": float(row_fractions.min()) if row_fractions.size else None,
        "row_case_fraction_max": float(row_fractions.max()) if row_fractions.size else None,
        "row_case_fraction_heterogeneous": bool(row_fractions.size and np.ptp(row_fractions) > 1e-12),
        "locked_study_n": locked_n,
        "runner_scalar_gwas_n_rule": "median(summary.gwas_n)",
        "runner_scalar_gwas_n": runner_gwas_n,
        "runner_scalar_eqtl_n_rule": "median(summary.eqtl_n)",
        "runner_scalar_eqtl_n": runner_eqtl_n,
        "input_sha256": {
            name: sha256_file(unit / name)
            for name in ("summary.csv", "ld.npz", "provenance.json", "overlap.json")
        },
    }


def materialise_locked_input(unit: Path, locked: Path) -> None:
    locked.mkdir(parents=True, exist_ok=True)
    shutil.copy2(unit / "summary.csv", locked / "summary.csv")
    archive = np.load(unit / "ld.npz", allow_pickle=False)
    ld = archive["ld"]
    ids = [str(value) for value in archive["variant_ids"].tolist()]
    with (locked / "ld.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([""] + ids)
        for identifier, values in zip(ids, ld, strict=True):
            writer.writerow([identifier] + [format(float(value), ".17g") for value in values])


def run_one(
    unit: Path,
    gate: Path,
    outcome_lock: dict[str, Any],
    runner: Path,
    out_root: Path,
    repeats: int,
    timeout_seconds: int,
) -> dict[str, Any]:
    audit = audit_unit(unit, gate, outcome_lock)
    unit_out = out_root / "units" / unit.name
    write_json(unit_out / "input_audit.json", audit)
    if audit["audit_status"] != "PASS":
        return {**audit, "run_status": "NOT_RUN", "run_code": "E_INPUT_AUDIT_FAILED"}
    locked = unit_out / "locked_input"
    materialise_locked_input(unit, locked)
    binding = {
        "summary_sha256": sha256_file(locked / "summary.csv"),
        "ld_csv_sha256": sha256_file(locked / "ld.csv"),
        "runner_sha256": sha256_file(runner),
        "case_fraction": float(outcome_lock["case_fraction"]),
    }
    write_json(locked / "binding.json", binding)
    repeat_rows: list[dict[str, Any]] = []
    for repeat in range(1, repeats + 1):
        run_out = unit_out / f"run_{repeat}"
        run_out.mkdir(parents=True, exist_ok=True)
        time_file = run_out / "gnu_time.txt"
        stdout_path, stderr_path = run_out / "stdout.log", run_out / "stderr.log"
        command = [
            "timeout", str(timeout_seconds), "/usr/bin/time", "-v", "-o", str(time_file),
            "Rscript", str(runner), str(locked), str(outcome_lock["case_fraction"]), str(run_out),
        ]
        started = time.time()
        timed_out = False
        try:
            completed = subprocess.run(command, capture_output=True, text=True, check=False)
            returncode = completed.returncode
            stdout_path.write_text(completed.stdout, encoding="utf-8")
            stderr_path.write_text(completed.stderr, encoding="utf-8")
        except subprocess.TimeoutExpired as error:
            timed_out, returncode = True, 124
            stdout_path.write_text(error.stdout or "", encoding="utf-8")
            stderr_path.write_text(error.stderr or "", encoding="utf-8")
        status, code, reason = classify_run(returncode, run_out / "status.csv", timed_out)
        hashes = {name: sha256_file(run_out / name) for name in RESULT_FILES if (run_out / name).is_file()}
        row = {
            "repeat": repeat,
            "status": status,
            "code": code,
            "reason": reason,
            "returncode": returncode,
            "wall_seconds_python": round(time.time() - started, 6),
            "result_sha256": hashes,
            **parse_time_metrics(time_file),
        }
        repeat_rows.append(row)
        write_json(run_out / "benchmark_status.json", row)
        if status != "COMPLETED":
            break
    comparable = len(repeat_rows) == repeats and all(row["status"] == "COMPLETED" for row in repeat_rows)
    byte_identical = bool(
        comparable and all(repeat_rows[index]["result_sha256"] == repeat_rows[0]["result_sha256"]
                           for index in range(1, len(repeat_rows)))
    )
    result = {
        **audit,
        "run_status": repeat_rows[-1]["status"],
        "run_code": repeat_rows[-1]["code"],
        "repeats_requested": repeats,
        "repeats_completed": sum(row["status"] == "COMPLETED" for row in repeat_rows),
        "byte_identical_repeats": byte_identical if comparable else None,
        "runs": repeat_rows,
    }
    write_json(unit_out / "unit_benchmark.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke-root", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--metadata-lock", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--timeout-seconds", type=int, default=3600)
    parser.add_argument("--unit", action="append", default=[])
    args = parser.parse_args()
    if args.jobs < 1 or args.repeats < 1 or args.timeout_seconds < 1:
        parser.error("jobs, repeats, and timeout-seconds must be positive")
    lock = read_json(args.metadata_lock)
    prepared = args.smoke_root / "prepared"
    units = [path for path in sorted(prepared.iterdir()) if path.is_dir()]
    if args.unit:
        selected = set(args.unit)
        units = [path for path in units if path.name in selected]
    if not units:
        raise SystemExit("E_NO_UNITS_SELECTED")
    args.out.mkdir(parents=True, exist_ok=True)
    configuration = {
        "smoke_root": str(args.smoke_root.resolve()),
        "runner": str(args.runner.resolve()),
        "runner_sha256": sha256_file(args.runner),
        "metadata_lock": str(args.metadata_lock.resolve()),
        "metadata_lock_sha256": sha256_file(args.metadata_lock),
        "jobs": args.jobs,
        "repeats": args.repeats,
        "timeout_seconds": args.timeout_seconds,
        "units": [path.name for path in units],
        "pid": os.getpid(),
    }
    write_json(args.out / "run_configuration.json", configuration)

    def invoke(unit: Path) -> dict[str, Any]:
        outcome = unit.name.rsplit("__", 1)[-1]
        return run_one(
            unit, args.smoke_root / "gates" / "units" / unit.name / "gate.json",
            lock["outcomes"][outcome], args.runner, args.out, args.repeats, args.timeout_seconds,
        )

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(invoke, units))
    write_json(args.out / "benchmark_summary.json", {"configuration": configuration, "units": results})
    fields = (
        "unit_id", "gate_status", "audit_status", "n_variants", "ld_rank_from_gate",
        "runner_scalar_eqtl_n", "runner_scalar_gwas_n", "locked_case_fraction",
        "row_case_fraction_heterogeneous", "run_status", "run_code", "repeats_completed",
        "byte_identical_repeats",
    )
    with (args.out / "benchmark_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            writer.writerow({field: result.get(field) for field in fields})
    if any(result["audit_status"] != "PASS" or result["run_status"] != "COMPLETED" for result in results):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
