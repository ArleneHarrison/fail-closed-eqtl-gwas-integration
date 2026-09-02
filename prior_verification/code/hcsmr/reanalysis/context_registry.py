"""Evidence-class registry and conservative context-comparison guardrails."""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Sequence

EvidenceClass = Literal["immune_population", "immune_stimulation", "bulk_tissue", "proxy"]


@dataclass(frozen=True)
class ContextRecord:
    name: str
    evidence_class: EvidenceClass
    biological_source: str
    covariance_available: bool
    source_study: str


@dataclass(frozen=True)
class ContextDecision:
    analysis_mode: Literal["joint", "separate"]
    can_claim_cardiac_cell_state_localization: bool


def load_context_registry(path: str | Path) -> list[ContextRecord]:
    """Load context metadata from a CSV with the ContextRecord field names."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        return [
            ContextRecord(
                name=row["name"],
                evidence_class=row["evidence_class"],  # type: ignore[arg-type]
                biological_source=row["biological_source"],
                covariance_available=row["covariance_available"].strip().lower() == "true",
                source_study=row["source_study"],
            )
            for row in rows
        ]


def validate_context_pair(
    records: Sequence[ContextRecord], left: str, right: str
) -> ContextDecision:
    """Allow joint modelling only for comparable, covariance-aware contexts.

    This function never treats PBMC or bulk tissue contexts as evidence of a
    cardiac cell state. That claim requires a dedicated cardiac cell-state
    eQTL resource and is intentionally outside this registry's scope.
    """
    by_name = {record.name: record for record in records}
    left_record = by_name[left]
    right_record = by_name[right]
    same_class = left_record.evidence_class == right_record.evidence_class
    same_source = left_record.biological_source == right_record.biological_source
    joint = (
        same_class
        and same_source
        and left_record.covariance_available
        and right_record.covariance_available
    )
    return ContextDecision(
        analysis_mode="joint" if joint else "separate",
        can_claim_cardiac_cell_state_localization=False,
    )
