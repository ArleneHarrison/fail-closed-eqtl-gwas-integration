"""Genome-build contracts that prevent unverified cross-resource inference."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import pandas as pd


@dataclass(frozen=True)
class SourceCoordinateRecord:
    source: str
    genome_build: str
    role: str


@dataclass(frozen=True)
class CoordinateDecision:
    status: str
    coordinate_pass: bool
    required_action: str


EQTL_SOURCE_BY_CONTEXT_PREFIX = {
    "GTEx_": "GTEx v8",
    "Fairfax_": "Fairfax eQTL Catalogue",
    "OneK1K_": "OneK1K eQTL Catalogue",
    "Alasoo_": "Alasoo eQTL Catalogue",
    "BLUEPRINT_": "BLUEPRINT eQTL Catalogue",
}
GWAS_SOURCE_BY_OUTCOME = {
    "CAD": "CAD Aragam 2022",
    "HF": "HERMES HF",
    "AF": "Nielsen AF 2018",
    "IS": "GIGASTROKE IS",
}


def validate_source_pair(
    records: Sequence[SourceCoordinateRecord], eqtl_source: str, gwas_source: str
) -> CoordinateDecision:
    by_source = {record.source: record for record in records}
    eqtl = by_source.get(eqtl_source)
    gwas = by_source.get(gwas_source)
    if eqtl is None or gwas is None or "unknown" in {eqtl.genome_build, gwas.genome_build}:
        return CoordinateDecision("unknown_build", False, "document source build and allele convention")
    if eqtl.genome_build != gwas.genome_build:
        return CoordinateDecision(
            "needs_liftover",
            False,
            f"liftover one source to a single build and re-harmonize ref/alt alleles ({eqtl.genome_build} vs {gwas.genome_build})",
        )
    return CoordinateDecision(
        "build_matched_pending_allele_audit",
        False,
        "complete allele, palindromic-SNP, and reference-panel match audit",
    )


def _source_for_context(context: str) -> str:
    for prefix, source in EQTL_SOURCE_BY_CONTEXT_PREFIX.items():
        if context.startswith(prefix):
            return source
    return "unknown_eQTL_source"


def build_coordinate_audit(
    candidates: pd.DataFrame, records: Sequence[SourceCoordinateRecord]
) -> pd.DataFrame:
    """Assign a coordinate decision to every candidate without silent defaults."""
    rows = []
    for candidate_id in candidates["candidate_id"].astype(str):
        outcome, _, context = candidate_id.split("|", maxsplit=2)
        eqtl_source = _source_for_context(context)
        gwas_source = GWAS_SOURCE_BY_OUTCOME.get(outcome, "unknown_GWAS_source")
        decision = validate_source_pair(records, eqtl_source, gwas_source)
        rows.append(
            {
                "candidate_id": candidate_id,
                "eqtl_source": eqtl_source,
                "gwas_source": gwas_source,
                "status": decision.status,
                "coordinate_pass": decision.coordinate_pass,
                "required_action": decision.required_action,
            }
        )
    return pd.DataFrame(rows)
