"""Fail-closed, machine-readable preflight checks for regional analyses.

This module deliberately does not condition, prune, impute, or otherwise
"repair" inputs.  It records every failed contract and returns an
``INELIGIBLE`` status when any contract fails.  A caller must provide the
numerical stopping threshold before a model is invoked.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Iterable, Mapping

import numpy as np


REQUIRED_PROVENANCE_FIELDS = (
    "url", "version_or_build", "license_or_access", "retrieval_date",
    "sha256", "byte_count", "structural_check",
)
ACCEPTED_ALIGNMENT_STATUSES = {"aligned", "flipped"}


@dataclass(frozen=True)
class GateFailure:
    code: str
    message: str
    affected_rows: int | None = None


def _finite_positive(value: Any) -> bool:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return False
    return bool(np.isfinite(numeric) and numeric > 0)


def _first_present(row: Mapping[str, Any], names: Iterable[str]) -> Any:
    for name in names:
        if name in row:
            return row[name]
    return None


def _failure(code: str, message: str, affected_rows: int | None = None) -> GateFailure:
    return GateFailure(code, message, affected_rows)


def prepare_analysis_rows(
    rows: list[Mapping[str, Any]],
    *,
    ld_variant_ids: list[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Create a deterministic analysis view while retaining a row-level audit.

    Rows absent from the locked LD identifier set are excluded explicitly. Exact
    duplicate mappings are collapsed by stable input order. Discordant duplicate
    mappings remain in the analysis view so that ``E_VARIANT_DUPLICATE`` blocks
    them rather than choosing among conflicting values.
    """
    ld_set = {str(value) for value in ld_variant_ids}
    selected: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    first_by_identifier: dict[str, dict[str, Any]] = {}
    for index, source in enumerate(rows):
        row = dict(source)
        identifier = str(_first_present(row, ("target_variant", "variant", "snp")) or "")
        if identifier not in ld_set:
            decision = "excluded_not_in_locked_ld"
        elif identifier not in first_by_identifier:
            first_by_identifier[identifier] = row
            selected.append(row)
            decision = "retained_first_occurrence"
        elif row == first_by_identifier[identifier]:
            decision = "excluded_exact_duplicate"
        else:
            selected.append(row)
            decision = "retained_discordant_duplicate_for_gate_failure"
        audit.append({"source_row_index": index, "target_variant": identifier, "decision": decision})
    return selected, audit


def validate_provenance(records: Iterable[Mapping[str, Any]]) -> list[GateFailure]:
    """Require a complete immutable-source record for each declared input."""
    failures: list[GateFailure] = []
    for index, record in enumerate(records):
        missing = [field for field in REQUIRED_PROVENANCE_FIELDS if not record.get(field)]
        if missing:
            label = str(record.get("label", f"input_{index}"))
            failures.append(_failure(
                "E_PROVENANCE_INCOMPLETE",
                f"{label} is missing provenance fields: {','.join(missing)}",
            ))
    return failures


def validate_summary_rows(
    rows: list[Mapping[str, Any]],
    *,
    enforce_constant_per_row_case_fraction: bool = True,
) -> list[GateFailure]:
    """Validate per-variant fields without discarding any input row."""
    failures: list[GateFailure] = []
    identifiers = [str(_first_present(row, ("target_variant", "variant", "snp")) or "") for row in rows]
    duplicate_count = len(identifiers) - len(set(identifiers))
    if not identifiers or any(not value for value in identifiers):
        failures.append(_failure("E_VARIANT_ID_INVALID", "summary rows contain a missing variant identifier"))
    if duplicate_count:
        failures.append(_failure("E_VARIANT_DUPLICATE", "summary rows contain non-unique variant identifiers", duplicate_count))

    checks = (
        ("E_EQTL_SE_INVALID", "eQTL standard error", ("eqtl_se", "se_eqtl", "se")),
        ("E_GWAS_SE_INVALID", "GWAS standard error", ("gwas_se", "se_gwas", "standard_error")),
        ("E_MAF_INVALID", "minor-allele frequency", ("maf", "MAF")),
        ("E_EQTL_N_INVALID", "eQTL sample size", ("eqtl_n", "eqtl_an", "an", "N_eqtl")),
        ("E_GWAS_N_INVALID", "GWAS sample size", ("gwas_n", "N", "n")),
    )
    for code, label, fields in checks:
        invalid = sum(not _finite_positive(_first_present(row, fields)) for row in rows)
        if invalid:
            failures.append(_failure(code, f"{label} is absent, non-finite, or non-positive", invalid))
    maf_invalid = sum(
        not (0 < float(_first_present(row, ("maf", "MAF"))) <= 0.5)
        if _finite_positive(_first_present(row, ("maf", "MAF"))) else True
        for row in rows
    )
    if maf_invalid:
        failures = [failure for failure in failures if failure.code != "E_MAF_INVALID"]
        failures.append(_failure("E_MAF_INVALID", "minor-allele frequency is absent or outside (0, 0.5]", maf_invalid))

    alignment_bad = sum(
        str(row.get("alignment_status", "")).lower() not in ACCEPTED_ALIGNMENT_STATUSES
        for row in rows
    )
    if alignment_bad:
        failures.append(_failure("E_ALLELE_UNRESOLVED", "effect-allele orientation is unresolved", alignment_bad))
    vcf_absent = sum(str(row.get("vcf_status", "")).lower() in {"position_absent", "vcf_absent"} for row in rows)
    if vcf_absent:
        failures.append(_failure("E_VCF_ABSENT", "a requested variant is absent from the LD VCF", vcf_absent))
    case_values = [_first_present(row, ("gwas_cases", "cases")) for row in rows]
    if any(value is not None for value in case_values):
        n_values = [_first_present(row, ("gwas_n", "N", "n")) for row in rows]
        invalid_cases = sum(
            not _finite_positive(cases) or not _finite_positive(sample_size) or float(cases) > float(sample_size)
            for cases, sample_size in zip(case_values, n_values, strict=True)
        )
        if invalid_cases:
            failures.append(_failure("E_GWAS_CASES_INVALID", "GWAS cases are absent, non-finite, or exceed the row sample size", invalid_cases))
        else:
            per_row_fractions = {round(float(cases) / float(sample_size), 12) for cases, sample_size in zip(case_values, n_values, strict=True)}
            if enforce_constant_per_row_case_fraction and len(per_row_fractions) != 1:
                failures.append(_failure(
                    "E_CASE_FRACTION_HETEROGENEOUS",
                    "per-row GWAS case fractions are heterogeneous; a single case fraction is not valid", len(rows),
                ))
    return failures


def validate_ld(
    ld: np.ndarray,
    *,
    max_condition_number: float | None,
    rank_policy: str = "require_full_rank",
    symmetry_tolerance: float = 1e-8,
    diagonal_tolerance: float = 1e-8,
) -> tuple[list[GateFailure], dict[str, float | int | bool | None]]:
    """Validate a supplied LD matrix and report, never modify, its diagnostics."""
    if rank_policy not in {"require_full_rank", "diagnostic"}:
        raise ValueError("rank_policy must be 'require_full_rank' or 'diagnostic'")
    matrix = np.asarray(ld, dtype=float)
    diagnostics: dict[str, float | int | bool | None] = {
        "n_variants": int(matrix.shape[0]) if matrix.ndim == 2 else 0,
        "finite": bool(np.isfinite(matrix).all()),
        "symmetric_max_abs_error": None,
        "diagonal_max_abs_error_from_one": None,
        "min_eigenvalue": None,
        "max_eigenvalue": None,
        "smallest_absolute_eigenvalue": None,
        "numerical_rank": None,
        "rank_tolerance": None,
        "rank_deficient": None,
        "positive_semidefinite": None,
        "psd_tolerance": None,
        "condition_number_2": None,
        "condition_number_is_infinite": None,
        "max_condition_number": max_condition_number,
        "rank_policy": rank_policy,
        "condition_number_threshold_applied": None,
    }
    failures: list[GateFailure] = []
    if matrix.ndim != 2 or matrix.shape[0] == 0 or matrix.shape[0] != matrix.shape[1]:
        return [_failure("E_LD_SHAPE", "LD matrix is not a non-empty square matrix")], diagnostics
    if not diagnostics["finite"]:
        return [_failure("E_LD_NONFINITE", "LD matrix contains non-finite values")], diagnostics
    symmetric_error = float(np.max(np.abs(matrix - matrix.T)))
    diagonal_error = float(np.max(np.abs(np.diag(matrix) - 1.0)))
    diagnostics["symmetric_max_abs_error"] = symmetric_error
    diagnostics["diagonal_max_abs_error_from_one"] = diagonal_error
    if symmetric_error > symmetry_tolerance:
        failures.append(_failure("E_LD_ASYMMETRIC", f"LD asymmetry {symmetric_error:.3g} exceeds tolerance"))
    if diagonal_error > diagonal_tolerance:
        failures.append(_failure("E_LD_DIAGONAL_INVALID", f"LD diagonal error {diagonal_error:.3g} exceeds tolerance"))
    eigenvalues = np.linalg.eigvalsh(matrix)
    minimum, maximum = float(eigenvalues[0]), float(eigenvalues[-1])
    absolute_eigenvalues = np.abs(eigenvalues)
    largest_absolute = float(np.max(absolute_eigenvalues))
    smallest_absolute = float(np.min(absolute_eigenvalues))
    rank_tolerance = float(matrix.shape[0] * np.finfo(float).eps * largest_absolute)
    numerical_rank = int(np.count_nonzero(absolute_eigenvalues > rank_tolerance))
    rank_deficient = numerical_rank < matrix.shape[0]
    positive_semidefinite = minimum >= -rank_tolerance
    condition_number = None if rank_deficient else float(largest_absolute / max(smallest_absolute, np.finfo(float).tiny))
    diagnostics.update(
        min_eigenvalue=minimum,
        max_eigenvalue=maximum,
        smallest_absolute_eigenvalue=smallest_absolute,
        numerical_rank=numerical_rank,
        rank_tolerance=rank_tolerance,
        rank_deficient=rank_deficient,
        positive_semidefinite=positive_semidefinite,
        psd_tolerance=rank_tolerance,
        condition_number_2=condition_number,
        condition_number_is_infinite=rank_deficient,
    )
    if not positive_semidefinite:
        failures.append(_failure(
            "E_LD_NOT_PSD",
            f"LD minimum eigenvalue {minimum:.6g} is below the numerical PSD tolerance {-rank_tolerance:.6g}",
        ))
    if max_condition_number is None:
        failures.append(_failure("E_LD_THRESHOLD_UNSPECIFIED", "no prospective LD condition-number stop threshold was supplied"))
        diagnostics["condition_number_threshold_applied"] = False
    elif rank_deficient and rank_policy == "require_full_rank":
        failures.append(_failure(
            "E_LD_RANK_DEFICIENT",
            f"LD matrix has numerical rank {numerical_rank} of {matrix.shape[0]} at tolerance {rank_tolerance:.6g}",
        ))
        diagnostics["condition_number_threshold_applied"] = False
    elif rank_deficient:
        diagnostics["condition_number_threshold_applied"] = False
    elif condition_number is not None and condition_number > max_condition_number:
        diagnostics["condition_number_threshold_applied"] = True
        failures.append(_failure(
            "E_LD_ILL_CONDITIONED",
            f"LD condition number {condition_number:.6g} exceeds the supplied stop threshold {max_condition_number:.6g}",
        ))
    else:
        diagnostics["condition_number_threshold_applied"] = True
    return failures, diagnostics


def run_preflight(
    *,
    rows: list[Mapping[str, Any]],
    ld: np.ndarray,
    provenance: list[Mapping[str, Any]],
    eqtl_trait_type: str,
    gwas_trait_type: str,
    case_fraction: float | None,
    max_condition_number: float | None,
    ld_variant_ids: list[str] | None = None,
    case_fraction_policy: str = "require_constant_per_row",
    ld_rank_policy: str = "require_full_rank",
) -> dict[str, Any]:
    """Run all contracts and emit one stable, machine-readable status object."""
    if case_fraction_policy not in {"require_constant_per_row", "study_level"}:
        raise ValueError("case_fraction_policy must be 'require_constant_per_row' or 'study_level'")
    failures = validate_provenance(provenance) + validate_summary_rows(
        rows,
        enforce_constant_per_row_case_fraction=case_fraction_policy == "require_constant_per_row",
    )
    row_variant_order = [str(_first_present(row, ("target_variant", "variant", "snp")) or "") for row in rows]
    row_variant_ids = set(row_variant_order)
    if ld_variant_ids is None:
        failures.append(_failure("E_LD_VARIANT_IDS_MISSING", "LD archive has no auditable variant identifier order"))
    else:
        ld_variant_order = [str(value) for value in ld_variant_ids]
        ld_variant_set = set(ld_variant_order)
        duplicate_ld_ids = len(ld_variant_order) - len(ld_variant_set)
        if duplicate_ld_ids:
            failures.append(_failure(
                "E_LD_VARIANT_DUPLICATE",
                "LD archive contains non-unique variant identifiers",
                duplicate_ld_ids,
            ))
        matrix = np.asarray(ld)
        matrix_dimension = int(matrix.shape[0]) if matrix.ndim == 2 else None
        if matrix_dimension is None or len(ld_variant_order) != matrix_dimension:
            failures.append(_failure(
                "E_LD_VARIANT_COUNT_MISMATCH",
                "LD identifier count does not equal the LD matrix dimension",
                abs(len(ld_variant_order) - matrix_dimension) if matrix_dimension is not None else None,
            ))
        mismatch = row_variant_ids.symmetric_difference(ld_variant_set)
        if mismatch:
            failures.append(_failure(
                "E_LD_VARIANT_SET_MISMATCH",
                "summary and LD archive do not contain the same unique variant set", len(mismatch),
            ))
        elif row_variant_order != ld_variant_order:
            failures.append(_failure(
                "E_LD_VARIANT_ORDER_MISMATCH",
                "summary-row order does not exactly match the locked LD identifier order",
                sum(left != right for left, right in zip(row_variant_order, ld_variant_order, strict=False))
                + abs(len(row_variant_order) - len(ld_variant_order)),
            ))
    if eqtl_trait_type != "quant":
        failures.append(_failure("E_EQTL_TRAIT_TYPE_INVALID", "eQTL trait type must be quant"))
    if gwas_trait_type != "cc":
        failures.append(_failure("E_GWAS_TRAIT_TYPE_INVALID", "GWAS trait type must be cc"))
    if case_fraction is None or not (0 < float(case_fraction) < 1):
        failures.append(_failure("E_CASE_FRACTION_INVALID", "case fraction must be strictly between zero and one"))
    ld_failures, diagnostics = validate_ld(
        ld, max_condition_number=max_condition_number, rank_policy=ld_rank_policy,
    )
    failures.extend(ld_failures)
    return {
        "schema_version": "1.2",
        "status": "READY" if not failures else "INELIGIBLE",
        "status_codes": [failure.code for failure in failures] or ["OK"],
        "failures": [asdict(failure) for failure in failures],
        "n_summary_rows": len(rows), "n_ld_variant_ids": len(ld_variant_ids) if ld_variant_ids is not None else None,
        "ld_diagnostics": diagnostics,
        "case_fraction_policy": case_fraction_policy,
        "decision_scope": "declared_gate_policy_conformance_only",
        "model_eligibility": None,
        "actions": [] if not failures else ["Do not run coloc or apply an unrecorded input-repair policy."],
    }
