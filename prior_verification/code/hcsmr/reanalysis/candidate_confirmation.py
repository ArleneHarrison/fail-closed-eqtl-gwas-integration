"""Assemble conservative, auditable primary-candidate decisions."""
from __future__ import annotations

from typing import Sequence

import pandas as pd

from .confirmation import evaluate_candidate
from .context_registry import ContextRecord


def build_confirmation_table(
    diagnostics: pd.DataFrame,
    contexts: Sequence[ContextRecord],
    coloc: pd.DataFrame,
    coordinates: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Return one all-or-nothing confirmation row per diagnostic candidate.

    The current audit intentionally keeps convergence, sensitivity, and
    multiplicity gates closed until their locked reanalysis outputs exist.
    """
    required = {"candidate_id", "instrument_pass"}
    missing = required.difference(diagnostics.columns)
    if missing:
        raise ValueError(f"diagnostics missing required columns: {sorted(missing)}")
    context_by_name = {record.name: record for record in contexts}
    coloc_by_region = (
        coloc.set_index("region_id").to_dict("index") if not coloc.empty else {}
    )
    coordinate_by_candidate = (
        coordinates.set_index("candidate_id").to_dict("index")
        if coordinates is not None and not coordinates.empty
        else {}
    )
    rows = []
    for diagnostic in diagnostics.to_dict("records"):
        candidate_id = str(diagnostic["candidate_id"])
        parts = candidate_id.split("|", maxsplit=2)
        context_name = parts[2] if len(parts) == 3 else ""
        context = context_by_name.get(context_name)
        coloc_record = coloc_by_region.get(candidate_id, {})
        coordinate_record = coordinate_by_candidate.get(candidate_id, {})
        coloc_status = str(coloc_record.get("status", "missing_coloc_record"))
        record = {
            "evidence_class_pass": context is not None,
            "coordinate_pass": bool(coordinate_record.get("coordinate_pass", False)),
            "instrument_pass": bool(diagnostic["instrument_pass"]),
            "convergence_pass": False,
            "sensitivity_pass": False,
            "multiplicity_pass": False,
            "coloc_pass": coloc_status == "confirmed",
        }
        decision = evaluate_candidate(record)
        rows.append(
            {
                "candidate_id": candidate_id,
                "context_evidence_class": context.evidence_class if context else "unknown",
                "context_source": context.source_study if context else "unknown",
                "coordinate_status": coordinate_record.get("status", "missing_coordinate_audit"),
                "coordinate_required_action": coordinate_record.get(
                    "required_action", "run_coordinate_audit"
                ),
                "coloc_status": coloc_status,
                "coloc_ineligibility_reason": coloc_record.get(
                    "ineligibility_reason", "missing_coloc_record"
                ),
                **record,
                "primary_candidate": decision.primary_candidate,
                "failed_criteria": ";".join(decision.failed_criteria),
            }
        )
    return pd.DataFrame(rows)
