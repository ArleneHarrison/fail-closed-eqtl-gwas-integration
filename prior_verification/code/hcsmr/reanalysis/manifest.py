"""Validation for locked analysis manifests."""
from __future__ import annotations

from typing import Mapping

from .confirmation import ValidationResult


def validate_manifest(manifest: Mapping[str, object]) -> ValidationResult:
    """Check that a manifest identifies the run, version, and hashed inputs."""
    errors: list[str] = []
    if not manifest.get("locked_run_id"):
        errors.append("locked_run_id")
    if not manifest.get("analysis_version"):
        errors.append("analysis_version")
    input_hashes = manifest.get("input_hashes")
    if not isinstance(input_hashes, Mapping) or not input_hashes:
        errors.append("input_hashes")
    return ValidationResult(valid=not errors, errors=tuple(errors))
