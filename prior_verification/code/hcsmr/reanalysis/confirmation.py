"""Candidate confirmation gates for locked reanalysis outputs."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping


@dataclass(frozen=True)
class CandidateDecision:
    primary_candidate: bool
    failed_criteria: tuple[str, ...]


@dataclass(frozen=True)
class ValidationResult:
    valid: bool
    errors: tuple[str, ...] = ()


CRITERION_LABELS = {
    "evidence_class_pass": "evidence_class",
    "coordinate_pass": "coordinate_harmonization",
    "instrument_pass": "instrument_strength",
    "convergence_pass": "convergence",
    "sensitivity_pass": "sensitivity",
    "multiplicity_pass": "multiplicity_control",
    "coloc_pass": "formal_colocalization",
}


def evaluate_candidate(record: Mapping[str, object]) -> CandidateDecision:
    """Return an all-or-nothing primary-candidate decision."""
    failed = tuple(
        label for key, label in CRITERION_LABELS.items() if not bool(record.get(key))
    )
    return CandidateDecision(primary_candidate=not failed, failed_criteria=failed)


def validate_manuscript_values(values: Mapping[str, object]) -> ValidationResult:
    """Require a locked run identifier and manifest hash for numerical claims."""
    errors = tuple(
        field
        for field in ("locked_run_id", "input_manifest_hash")
        if not values.get(field)
    )
    return ValidationResult(valid=not errors, errors=errors)
