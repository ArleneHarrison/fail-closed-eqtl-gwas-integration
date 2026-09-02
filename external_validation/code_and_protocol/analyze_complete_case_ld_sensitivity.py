#!/usr/bin/env python3
"""Describe complete-case LD sensitivity without changing the primary gate.

This is a secondary, read-only diagnostic.  It compares the frozen primary
mean-filled LD matrix with the separately materialized complete-case matrix.
It neither chooses a difference threshold nor writes beneath the prepared or
primary-gate roots, and it cannot alter gate eligibility or model dispatch.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from external_gate import validate_ld  # noqa: E402


OUTPUT_COLUMNS = [
    "unit_id", "region_id", "anchor_outcome", "tissue_id", "outcome_id",
    "gene_id", "n_variants", "complete_case_samples", "variant_order_identical",
    "analysis_status", "primary_gate_status", "primary_gate_status_codes",
    "primary_ld_integrity_status", "primary_ld_failure_codes",
    "complete_case_ld_integrity_status", "complete_case_ld_failure_codes",
    "primary_vs_complete_case_ld_status_concordant", "primary_ld_rank",
    "complete_case_ld_rank", "primary_vs_complete_case_rank_concordant",
    "primary_ld_condition_number_2", "primary_ld_condition_number_is_infinite",
    "complete_case_ld_condition_number_2", "complete_case_ld_condition_number_is_infinite",
    "complete_case_finite", "complete_case_symmetric_max_abs_error",
    "complete_case_diagonal_max_abs_error_from_one", "complete_case_min_eigenvalue",
    "complete_case_max_eigenvalue", "complete_case_rank_tolerance",
    "complete_case_positive_semidefinite", "relative_frobenius_distance",
    "max_absolute_off_diagonal_delta_r", "median_absolute_off_diagonal_delta_r",
    "p95_absolute_off_diagonal_delta_r",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def ordered_hash(identifiers: list[str]) -> str:
    return hashlib.sha256(("\n".join(identifiers) + "\n").encode("utf-8")).hexdigest()


def require(condition: bool, code: str, detail: str) -> None:
    if not condition:
        raise RuntimeError(f"{code}:{detail}")


def load_archive(path: Path) -> tuple[np.ndarray, list[str]]:
    require(path.is_file(), "E_SENSITIVITY_INPUT_MISSING", str(path))
    with np.load(path, allow_pickle=False) as archive:
        require(set(archive.files) == {"ld", "variant_ids"}, "E_LD_ARCHIVE_FIELDS", str(path))
        matrix = np.asarray(archive["ld"], dtype=np.float64)
        identifiers = [str(value) for value in archive["variant_ids"].tolist()]
    require(matrix.ndim == 2 and matrix.shape[0] == matrix.shape[1], "E_LD_SHAPE", str(path))
    require(matrix.shape[0] == len(identifiers), "E_LD_IDENTIFIER_DIMENSION", str(path))
    require(len(identifiers) == len(set(identifiers)), "E_LD_DUPLICATE_IDENTIFIER", str(path))
    return matrix, identifiers


def gate_ld_status(gate: dict[str, Any], n_variants: int) -> tuple[str, list[str]]:
    if n_variants == 0:
        return "N_A", []
    diagnostics = gate.get("ld_diagnostics")
    require(isinstance(diagnostics, dict), "E_PRIMARY_GATE_LD_DIAGNOSTICS_MISSING", str(gate.get("status")))
    require(int(diagnostics.get("n_variants", -1)) == n_variants,
            "E_PRIMARY_GATE_LD_DIMENSION_MISMATCH", f"{diagnostics.get('n_variants')}!={n_variants}")
    for key in ("finite", "symmetric_max_abs_error", "diagonal_max_abs_error_from_one",
                "numerical_rank", "rank_tolerance", "positive_semidefinite",
                "condition_number_is_infinite"):
        require(key in diagnostics, "E_PRIMARY_GATE_LD_DIAGNOSTIC_FIELD_MISSING", key)
    codes = sorted({str(code) for code in gate.get("status_codes", []) if str(code).startswith("E_LD_")})
    return ("FAIL" if codes else "PASS"), codes


def bool_text(value: bool | None) -> str:
    return "NA" if value is None else ("TRUE" if value else "FALSE")


def number_text(value: Any) -> str:
    if value is None:
        return "NA"
    number = float(value)
    if math.isnan(number):
        return "NA"
    if math.isinf(number):
        return "Inf" if number > 0 else "-Inf"
    return format(number, ".17g")


def percentile_linear(values: np.ndarray, percentile: float) -> float:
    try:
        return float(np.percentile(values, percentile, method="linear"))
    except TypeError:  # NumPy < 1.22
        return float(np.percentile(values, percentile, interpolation="linear"))


def difference_metrics(primary: np.ndarray, complete: np.ndarray) -> dict[str, str]:
    na = {
        "relative_frobenius_distance": "NA",
        "max_absolute_off_diagonal_delta_r": "NA",
        "median_absolute_off_diagonal_delta_r": "NA",
        "p95_absolute_off_diagonal_delta_r": "NA",
    }
    if primary.size == 0 or not np.isfinite(primary).all() or not np.isfinite(complete).all():
        return na
    denominator = float(np.linalg.norm(primary, ord="fro"))
    require(denominator > 0, "E_PRIMARY_LD_FROBENIUS_ZERO", "non-empty primary LD")
    delta = complete - primary
    relative = float(np.linalg.norm(delta, ord="fro") / denominator)
    if primary.shape[0] < 2:
        off_values = np.asarray([], dtype=float)
    else:
        off_values = np.abs(delta[np.triu_indices(primary.shape[0], k=1)])
    return {
        "relative_frobenius_distance": number_text(relative),
        "max_absolute_off_diagonal_delta_r": number_text(float(np.max(off_values))) if off_values.size else "NA",
        "median_absolute_off_diagonal_delta_r": number_text(float(np.median(off_values))) if off_values.size else "NA",
        "p95_absolute_off_diagonal_delta_r": number_text(percentile_linear(off_values, 95.0)) if off_values.size else "NA",
    }


def analyze_unit(
    attempt: dict[str, str], prepared_root: Path, gates_root: Path,
    max_condition_number: float,
) -> dict[str, str]:
    unit_id = str(attempt["unit_id"])
    unit = prepared_root / unit_id
    primary_path = unit / "ld.npz"
    complete_path = unit / "ld_complete_case.npz"
    trace_path = unit / "ld_missingness_trace.json"
    gate_path = gates_root / "units" / unit_id / "gate.json"
    for required_path in (primary_path, complete_path, trace_path, gate_path):
        require(required_path.is_file(), "E_SENSITIVITY_INPUT_MISSING", str(required_path))

    primary, primary_ids = load_archive(primary_path)
    complete, complete_ids = load_archive(complete_path)
    require(primary_ids == complete_ids, "E_COMPLETE_CASE_VARIANT_ORDER_MISMATCH", unit_id)
    require(primary.shape == complete.shape, "E_COMPLETE_CASE_LD_DIMENSION_MISMATCH", unit_id)

    trace = read_json(trace_path)
    require(isinstance(trace, dict), "E_MISSINGNESS_TRACE_INVALID", unit_id)
    require(trace.get("ordered_variant_sha256") == ordered_hash(primary_ids),
            "E_MISSINGNESS_TRACE_ORDER_HASH_MISMATCH", unit_id)
    complete_samples = trace.get("complete_case_samples")
    primary_samples = trace.get("primary_reference_samples")
    require(isinstance(complete_samples, int) and complete_samples >= 0,
            "E_COMPLETE_CASE_SAMPLE_COUNT_INVALID", unit_id)
    require(isinstance(primary_samples, int) and complete_samples <= primary_samples,
            "E_COMPLETE_CASE_SAMPLE_COUNT_INVALID", unit_id)

    gate = read_json(gate_path)
    require(isinstance(gate, dict), "E_PRIMARY_GATE_INVALID", unit_id)
    gate_input_hashes = gate.get("prepared_input_sha256")
    if primary_ids or isinstance(gate_input_hashes, dict):
        require(isinstance(gate_input_hashes, dict) and "ld.npz" in gate_input_hashes,
                "E_PRIMARY_GATE_INPUT_BINDING_MISSING", unit_id)
        require(gate_input_hashes["ld.npz"] == sha256_file(primary_path),
                "E_PRIMARY_GATE_INPUT_BINDING_MISMATCH", unit_id)
    primary_status, primary_failures = gate_ld_status(gate, len(primary_ids))
    primary_diag = gate.get("ld_diagnostics", {})

    base = {column: "NA" for column in OUTPUT_COLUMNS}
    for key in ("unit_id", "region_id", "anchor_outcome", "tissue_id", "outcome_id", "gene_id"):
        base[key] = str(attempt.get(key, "")) or "NA"
    base.update({
        "n_variants": str(len(primary_ids)),
        "complete_case_samples": str(complete_samples),
        "variant_order_identical": "TRUE",
        "primary_gate_status": str(gate.get("status", "NA")),
        "primary_gate_status_codes": ";".join(map(str, gate.get("status_codes", []))) or "NA",
        "primary_ld_integrity_status": primary_status,
        "primary_ld_failure_codes": ";".join(primary_failures) or "NA",
    })

    if not primary_ids:
        base.update({
            "analysis_status": "N_A_ZERO_OVERLAP",
            "complete_case_ld_integrity_status": "N_A",
            "complete_case_ld_failure_codes": "NA",
            "primary_vs_complete_case_ld_status_concordant": "NA",
            "primary_vs_complete_case_rank_concordant": "NA",
        })
        return base

    failures, complete_diag = validate_ld(
        complete, max_condition_number=max_condition_number, rank_policy="diagnostic",
    )
    cc_codes = sorted({failure.code for failure in failures})
    cc_status = "FAIL" if cc_codes else "PASS"
    primary_rank = int(primary_diag["numerical_rank"])
    complete_rank = complete_diag.get("numerical_rank")
    base.update({
        "analysis_status": "DESCRIPTIVE_COMPLETE_CASE_SENSITIVITY",
        "complete_case_ld_integrity_status": cc_status,
        "complete_case_ld_failure_codes": ";".join(cc_codes) or "NA",
        "primary_vs_complete_case_ld_status_concordant": bool_text(primary_status == cc_status),
        "primary_ld_rank": str(primary_rank),
        "complete_case_ld_rank": "NA" if complete_rank is None else str(int(complete_rank)),
        "primary_vs_complete_case_rank_concordant": (
            "NA" if complete_rank is None else bool_text(primary_rank == int(complete_rank))
        ),
        "primary_ld_condition_number_2": number_text(primary_diag.get("condition_number_2")),
        "primary_ld_condition_number_is_infinite": bool_text(primary_diag.get("condition_number_is_infinite")),
        "complete_case_ld_condition_number_2": number_text(complete_diag.get("condition_number_2")),
        "complete_case_ld_condition_number_is_infinite": bool_text(complete_diag.get("condition_number_is_infinite")),
        "complete_case_finite": bool_text(complete_diag.get("finite")),
        "complete_case_symmetric_max_abs_error": number_text(complete_diag.get("symmetric_max_abs_error")),
        "complete_case_diagonal_max_abs_error_from_one": number_text(complete_diag.get("diagonal_max_abs_error_from_one")),
        "complete_case_min_eigenvalue": number_text(complete_diag.get("min_eigenvalue")),
        "complete_case_max_eigenvalue": number_text(complete_diag.get("max_eigenvalue")),
        "complete_case_rank_tolerance": number_text(complete_diag.get("rank_tolerance")),
        "complete_case_positive_semidefinite": bool_text(complete_diag.get("positive_semidefinite")),
    })
    base.update(difference_metrics(primary, complete))
    return base


def _analyze_unit_task(payload: tuple[dict[str, str], Path, Path, float]) -> dict[str, str]:
    """Pickle-safe worker wrapper; executor.map preserves attempt-registry order."""
    return analyze_unit(*payload)


def write_tsv(path: Path, rows: Iterable[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def analyze(
    attempts_path: Path, prepared_root: Path, gates_root: Path, protocol_path: Path,
    output_dir: Path, expected_attempts: int = 192, jobs: int = 1,
) -> None:
    require(attempts_path.is_file(), "E_ATTEMPTS_MISSING", str(attempts_path))
    require(prepared_root.is_dir(), "E_PREPARED_ROOT_MISSING", str(prepared_root))
    require(gates_root.is_dir(), "E_GATES_ROOT_MISSING", str(gates_root))
    protocol = read_json(protocol_path)
    require(protocol.get("gate", {}).get("ld_rank_policy") == "diagnostic",
            "E_PROTOCOL_RANK_POLICY", str(protocol.get("gate", {}).get("ld_rank_policy")))
    max_condition_number = float(protocol["gate"]["max_condition_number"])

    attempts_frame = pd.read_csv(attempts_path, sep="\t", dtype=str, keep_default_na=False)
    require("unit_id" in attempts_frame.columns, "E_ATTEMPTS_COLUMNS", "unit_id")
    require(len(attempts_frame) == expected_attempts,
            "E_ATTEMPT_COUNT", f"observed={len(attempts_frame)} expected={expected_attempts}")
    require(attempts_frame["unit_id"].is_unique, "E_ATTEMPT_UNIT_DUPLICATE", "unit_id")
    require(jobs >= 1, "E_JOBS_INVALID", str(jobs))
    attempts = attempts_frame.to_dict(orient="records")

    parent = output_dir.resolve().parent
    parent.mkdir(parents=True, exist_ok=True)
    require(not output_dir.exists() or not any(output_dir.iterdir()),
            "E_OUTPUT_NOT_EMPTY", str(output_dir))
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.staging.", dir=parent))
    try:
        payloads = ((attempt, prepared_root, gates_root, max_condition_number) for attempt in attempts)
        if jobs == 1:
            rows = [_analyze_unit_task(payload) for payload in payloads]
        else:
            with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as executor:
                rows = list(executor.map(_analyze_unit_task, payloads))
        result_path = staging / "complete_case_ld_sensitivity.tsv"
        write_tsv(result_path, rows)
        counts = Counter(row["analysis_status"] for row in rows)
        status_counts = Counter(row["complete_case_ld_integrity_status"] for row in rows)
        summary = {
            "schema_version": "complete_case_ld_sensitivity_1.0",
            "status": "PASS",
            "analysis_role": "descriptive_secondary_sensitivity_only",
            "primary_gate_or_dispatch_mutation": False,
            "difference_threshold_selected": False,
            "rank_policy": "diagnostic",
            "max_condition_number_inherited_from_frozen_primary_protocol": max_condition_number,
            "normalized_frobenius_definition": "frobenius(complete_case-primary)/frobenius(primary)",
            "expected_attempts": expected_attempts,
            "observed_attempts": len(rows),
            "parallel_workers": jobs,
            "analysis_status_counts": dict(sorted(counts.items())),
            "complete_case_ld_integrity_status_counts": dict(sorted(status_counts.items())),
            "inputs": {
                "attempts_sha256": sha256_file(attempts_path),
                "protocol_sha256": sha256_file(protocol_path),
            },
            "implementation_sha256": sha256_file(Path(__file__).resolve()),
            "result_sha256": sha256_file(result_path),
        }
        summary_path = staging / "run_summary.json"
        summary_path.write_text(canonical_json(summary) + "\n", encoding="utf-8", newline="\n")
        manifest_path = staging / "manifest.sha256"
        manifest_path.write_text(
            f"{sha256_file(result_path)}  {result_path.name}\n"
            f"{sha256_file(summary_path)}  {summary_path.name}\n",
            encoding="utf-8", newline="\n",
        )
        if output_dir.exists():
            output_dir.rmdir()
        os.replace(staging, output_dir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attempts", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--gates", type=Path, required=True, help="Frozen primary gates_v1 directory")
    parser.add_argument("--protocol", type=Path, default=SCRIPT_DIR / "config/protocol.v2.json")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--expected-attempts", type=int, default=192)
    parser.add_argument("--jobs", type=int, default=1,
                        help="Independent unit workers; output order remains the frozen attempt order")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        analyze(args.attempts, args.prepared, args.gates, args.protocol, args.out_dir,
                args.expected_attempts, args.jobs)
    except Exception as error:
        print(str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
