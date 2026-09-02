"""Deterministic fault-injection and scalability benchmark for the preflight gate."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import platform
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
from hcsmr.reanalysis.preflight_gate import run_preflight


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def toeplitz_ld(size: int, rho: float = 0.2) -> np.ndarray:
    indices = np.arange(size)
    return rho ** np.abs(indices[:, None] - indices[None, :])


def base_contract(size: int) -> dict[str, Any]:
    variants = [f"chr11_{5_000_000 + index}_A_G" for index in range(size)]
    rows = [
        {
            "target_variant": variant,
            "eqtl_se": 0.02,
            "gwas_se": 0.03,
            "maf": 0.2,
            "eqtl_an": 838,
            "gwas_n": 1000,
            "gwas_cases": 100,
            "alignment_status": "aligned" if index % 2 == 0 else "flipped",
            "vcf_status": "exact_ref_alt_match",
        }
        for index, variant in enumerate(variants)
    ]
    provenance = [
        {
            "label": "synthetic fixture",
            "url": "fixture://complete-contract",
            "version_or_build": "synthetic-v1",
            "license_or_access": "CC0 synthetic fixture",
            "retrieval_date": "2026-08-28",
            "sha256": "0" * 64,
            "byte_count": 1,
            "structural_check": "generated in memory",
        }
    ]
    return {
        "rows": rows,
        "ld": toeplitz_ld(size),
        "ld_variant_ids": variants,
        "provenance": provenance,
        "eqtl_trait_type": "quant",
        "gwas_trait_type": "cc",
        "case_fraction": 0.1,
        "max_condition_number": 1e12,
    }


def change_row(field: str, value: Any) -> Callable[[dict[str, Any], int], None]:
    def mutation(contract: dict[str, Any], seed: int) -> None:
        contract["rows"][seed % len(contract["rows"])][field] = value
    return mutation


def mutation_registry() -> list[tuple[str, str, Callable[[dict[str, Any], int], None]]]:
    def provenance(contract: dict[str, Any], _seed: int) -> None:
        contract["provenance"][0]["sha256"] = ""

    def duplicate(contract: dict[str, Any], seed: int) -> None:
        index = seed % (len(contract["rows"]) - 1)
        contract["rows"][index + 1]["target_variant"] = contract["rows"][index]["target_variant"]

    def heterogeneous_cases(contract: dict[str, Any], seed: int) -> None:
        contract["rows"][seed % len(contract["rows"])]["gwas_cases"] = 101

    def missing_ld_ids(contract: dict[str, Any], _seed: int) -> None:
        contract["ld_variant_ids"] = None

    def mismatched_ld_ids(contract: dict[str, Any], seed: int) -> None:
        contract["ld_variant_ids"][seed % len(contract["ld_variant_ids"])] = "chr11_999999999_C_T"

    def duplicated_ld_ids(contract: dict[str, Any], seed: int) -> None:
        index = seed % (len(contract["ld_variant_ids"]) - 1)
        contract["ld_variant_ids"][index + 1] = contract["ld_variant_ids"][index]

    def shortened_ld_ids(contract: dict[str, Any], _seed: int) -> None:
        contract["ld_variant_ids"] = contract["ld_variant_ids"][:-1]

    def reordered_ld_ids(contract: dict[str, Any], _seed: int) -> None:
        contract["ld_variant_ids"][0], contract["ld_variant_ids"][1] = (
            contract["ld_variant_ids"][1], contract["ld_variant_ids"][0]
        )

    def asymmetric(contract: dict[str, Any], seed: int) -> None:
        index = seed % (len(contract["rows"]) - 1)
        contract["ld"][index, index + 1] += 1e-4

    def diagonal(contract: dict[str, Any], seed: int) -> None:
        index = seed % len(contract["rows"])
        contract["ld"][index, index] += 1e-4

    def nonfinite(contract: dict[str, Any], seed: int) -> None:
        index = seed % len(contract["rows"])
        contract["ld"][index, index] = np.nan

    def ill_conditioned(contract: dict[str, Any], _seed: int) -> None:
        contract["ld"] = np.eye(len(contract["rows"]))
        # Keep the block full-rank at the scale-aware tolerance while making
        # its standard 2-norm condition number exceed the declared threshold.
        contract["ld"][0, 1] = contract["ld"][1, 0] = 1.0 - 1e-13

    def not_psd(contract: dict[str, Any], _seed: int) -> None:
        contract["ld"] = np.eye(len(contract["rows"]))
        contract["ld"][:3, :3] = np.array([
            [1.0, -0.6, -0.6],
            [-0.6, 1.0, -0.6],
            [-0.6, -0.6, 1.0],
        ])

    return [
        ("provenance_incomplete", "E_PROVENANCE_INCOMPLETE", provenance),
        ("variant_id_missing", "E_VARIANT_ID_INVALID", change_row("target_variant", "")),
        ("variant_duplicate", "E_VARIANT_DUPLICATE", duplicate),
        ("eqtl_se_invalid", "E_EQTL_SE_INVALID", change_row("eqtl_se", 0)),
        ("gwas_se_invalid", "E_GWAS_SE_INVALID", change_row("gwas_se", -1)),
        ("maf_invalid", "E_MAF_INVALID", change_row("maf", 0.5001)),
        ("eqtl_n_invalid", "E_EQTL_N_INVALID", change_row("eqtl_an", "")),
        ("gwas_n_invalid", "E_GWAS_N_INVALID", change_row("gwas_n", 0)),
        ("allele_unresolved", "E_ALLELE_UNRESOLVED", change_row("alignment_status", "strand_ambiguous")),
        ("vcf_absent", "E_VCF_ABSENT", change_row("vcf_status", "position_absent")),
        ("case_fraction_heterogeneous", "E_CASE_FRACTION_HETEROGENEOUS", heterogeneous_cases),
        ("ld_variant_ids_missing", "E_LD_VARIANT_IDS_MISSING", missing_ld_ids),
        ("ld_variant_set_mismatch", "E_LD_VARIANT_SET_MISMATCH", mismatched_ld_ids),
        ("ld_variant_ids_duplicate", "E_LD_VARIANT_DUPLICATE", duplicated_ld_ids),
        ("ld_variant_count_mismatch", "E_LD_VARIANT_COUNT_MISMATCH", shortened_ld_ids),
        ("ld_variant_order_mismatch", "E_LD_VARIANT_ORDER_MISMATCH", reordered_ld_ids),
        ("eqtl_trait_invalid", "E_EQTL_TRAIT_TYPE_INVALID", lambda c, _s: c.update(eqtl_trait_type="other")),
        ("gwas_trait_invalid", "E_GWAS_TRAIT_TYPE_INVALID", lambda c, _s: c.update(gwas_trait_type="other")),
        ("case_fraction_invalid", "E_CASE_FRACTION_INVALID", lambda c, _s: c.update(case_fraction=1.0)),
        ("ld_asymmetric", "E_LD_ASYMMETRIC", asymmetric),
        ("ld_diagonal_invalid", "E_LD_DIAGONAL_INVALID", diagonal),
        ("ld_nonfinite", "E_LD_NONFINITE", nonfinite),
        ("ld_ill_conditioned", "E_LD_ILL_CONDITIONED", ill_conditioned),
        ("ld_not_psd", "E_LD_NOT_PSD", not_psd),
        ("ld_threshold_unspecified", "E_LD_THRESHOLD_UNSPECIFIED", lambda c, _s: c.update(max_condition_number=None)),
    ]


def execute(contract: dict[str, Any]) -> tuple[dict[str, Any], float]:
    started = time.perf_counter_ns()
    result = run_preflight(**contract)
    elapsed_ms = (time.perf_counter_ns() - started) / 1_000_000
    return result, elapsed_ms


def run_correctness(replicates: int, size: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    run_rows: list[dict[str, Any]] = []
    cases = [("valid_complete_contract", "OK", lambda _c, _s: None), *mutation_registry()]
    for case_name, expected_code, mutation in cases:
        for seed in range(replicates):
            contract = base_contract(size)
            mutation(contract, seed)
            result, runtime_ms = execute(contract)
            observed = list(result["status_codes"])
            expected_present = expected_code in observed
            run_rows.append({
                "case": case_name,
                "seed": seed,
                "n_variants": size,
                "expected_code": expected_code,
                "gate_status": result["status"],
                "observed_codes": ";".join(observed),
                "expected_code_present": expected_present,
                "runtime_ms": runtime_ms,
            })
    runs = pd.DataFrame(run_rows)
    summaries = []
    for case_name, group in runs.groupby("case", sort=False):
        successes = int(group["expected_code_present"].sum())
        trials = len(group)
        summaries.append({
            "case": case_name,
            "expected_code": group["expected_code"].iloc[0],
            "replicates": trials,
            "expected_code_detected": successes,
            "detection_rate": successes / trials,
            "median_runtime_ms": float(group["runtime_ms"].median()),
            "p95_runtime_ms": float(group["runtime_ms"].quantile(0.95)),
        })
    return runs, pd.DataFrame(summaries)


def run_boundaries(size: int) -> pd.DataFrame:
    cases: list[tuple[str, str, Callable[[dict[str, Any]], None]]] = []
    for delta in (5e-9, 2e-8):
        cases.append((
            f"asymmetry_delta_{delta:g}", "READY" if delta <= 1e-8 else "INELIGIBLE",
            lambda c, value=delta: c["ld"].__setitem__((0, 1), c["ld"][0, 1] + value),
        ))
        cases.append((
            f"diagonal_delta_{delta:g}", "READY" if delta <= 1e-8 else "INELIGIBLE",
            lambda c, value=delta: c["ld"].__setitem__((0, 0), 1.0 + value),
        ))
    for epsilon in (5e-12, 5e-13):
        def condition_mutation(contract: dict[str, Any], value: float = epsilon) -> None:
            contract["ld"] = np.eye(size)
            contract["ld"][0, 1] = contract["ld"][1, 0] = 1.0 - value
        cases.append((
            f"condition_block_epsilon_{epsilon:g}", "READY" if 2.0 / epsilon <= 1e12 else "INELIGIBLE",
            condition_mutation,
        ))
    rows = []
    for case_name, expected_status, mutation in cases:
        contract = base_contract(size)
        mutation(contract)
        result, runtime_ms = execute(contract)
        rows.append({
            "case": case_name,
            "expected_status": expected_status,
            "observed_status": result["status"],
            "status_codes": ";".join(result["status_codes"]),
            "condition_number": result["ld_diagnostics"].get("condition_number_2"),
            "symmetry_max_error": result["ld_diagnostics"].get("symmetric_max_abs_error"),
            "diagonal_max_error": result["ld_diagnostics"].get("diagonal_max_abs_error_from_one"),
            "expected_status_matched": result["status"] == expected_status,
            "runtime_ms": runtime_ms,
        })
    return pd.DataFrame(rows)


def run_scaling(sizes: list[int], repeats: int) -> pd.DataFrame:
    rows = []
    for size in sizes:
        for repeat in range(repeats):
            result, runtime_ms = execute(base_contract(size))
            rows.append({
                "n_variants": size,
                "repeat": repeat,
                "status": result["status"],
                "status_codes": ";".join(result["status_codes"]),
                "runtime_ms": runtime_ms,
            })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--correctness-replicates", type=int, default=100)
    parser.add_argument("--correctness-size", type=int, default=50)
    parser.add_argument("--scaling-sizes", default="10,50,100,250,500,1000")
    parser.add_argument("--scaling-repeats", type=int, default=5)
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    correctness_runs, correctness_summary = run_correctness(args.correctness_replicates, args.correctness_size)
    boundaries = run_boundaries(args.correctness_size)
    scaling = run_scaling([int(value) for value in args.scaling_sizes.split(",")], args.scaling_repeats)
    correctness_runs.to_csv(out_dir / "mutation_runs.tsv", sep="\t", index=False)
    correctness_summary.to_csv(out_dir / "mutation_summary.tsv", sep="\t", index=False)
    boundaries.to_csv(out_dir / "tolerance_boundary.tsv", sep="\t", index=False)
    scaling.to_csv(out_dir / "scaling_runs.tsv", sep="\t", index=False)
    valid = correctness_runs.loc[correctness_runs["case"] == "valid_complete_contract"]
    mutations = correctness_runs.loc[correctness_runs["case"] != "valid_complete_contract"]
    scaling_summary = (
        scaling.groupby("n_variants", as_index=False)
        .agg(median_runtime_ms=("runtime_ms", "median"), min_runtime_ms=("runtime_ms", "min"), max_runtime_ms=("runtime_ms", "max"))
    )
    scaling_summary.to_csv(out_dir / "scaling_summary.tsv", sep="\t", index=False)
    summary = {
        "schema_version": "1.0",
        "scope": "synthetic contract-mutation and computational scaling benchmark; no biological or model result",
        "correctness_size": args.correctness_size,
        "correctness_replicates_per_case": args.correctness_replicates,
        "valid_cases_ready": int((valid["gate_status"] == "READY").sum()),
        "valid_cases_total": len(valid),
        "fault_cases_expected_code_detected": int(mutations["expected_code_present"].sum()),
        "fault_cases_total": len(mutations),
        "all_boundary_expectations_matched": bool(boundaries["expected_status_matched"].all()),
        "scaling_all_ready": bool((scaling["status"] == "READY").all()),
        "total_runtime_seconds": time.perf_counter() - started,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "interpretation_boundary": "This validates declared software contracts and timing only; it does not estimate biological validity, statistical power, or superiority.",
    }
    write_json(out_dir / "benchmark_summary.json", summary)
    manifest = {
        "command": [sys.executable, *sys.argv],
        "files": [
            {"path": path.name, "byte_count": path.stat().st_size, "sha256": sha256_file(path)}
            for path in sorted(out_dir.iterdir()) if path.is_file() and path.name != "artifact_manifest.json"
        ],
    }
    write_json(out_dir / "artifact_manifest.json", manifest)
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
