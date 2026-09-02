"""Strict effect-allele alignment for regional eQTL/GWAS summary data."""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class AlignmentResult:
    status: str
    beta: float | None


def align_outcome_to_exposure(
    exposure_ref: str,
    exposure_alt: str,
    outcome_effect: str,
    outcome_other: str,
    outcome_beta: float,
) -> AlignmentResult:
    ref, alt = exposure_ref.upper(), exposure_alt.upper()
    effect, other = outcome_effect.upper(), outcome_other.upper()
    if len(ref) != 1 or len(alt) != 1 or len(effect) != 1 or len(other) != 1:
        return AlignmentResult("indel_unresolved", None)
    if {ref, alt} in ({"A", "T"}, {"C", "G"}):
        return AlignmentResult("palindromic_ambiguous", None)
    if effect == alt and other == ref:
        return AlignmentResult("aligned", float(outcome_beta))
    if effect == ref and other == alt:
        return AlignmentResult("flipped", -float(outcome_beta))
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if effect == complement.get(alt) and other == complement.get(ref):
        return AlignmentResult("strand_unresolved", None)
    if effect == complement.get(ref) and other == complement.get(alt):
        return AlignmentResult("strand_unresolved", None)
    return AlignmentResult("allele_mismatch", None)


def join_by_lifted_position(
    eqtl: pd.DataFrame,
    lifted: pd.DataFrame,
    gwas: pd.DataFrame,
    gwas_chrom_column: str,
    gwas_position_column: str,
) -> pd.DataFrame:
    """Join eQTL rows to coordinate-based GWAS data after liftover.

    The function preserves alleles for a subsequent orientation decision. It
    does not assume that position alone determines a unique biallelic variant.
    """
    mapped = lifted.loc[lifted["status"] == "mapped"].copy()
    components = mapped["target_variant"].str.split("_", n=3, expand=True)
    mapped["lifted_chrom"] = components[0].str.removeprefix("chr")
    mapped["lifted_position"] = pd.to_numeric(components[1], errors="raise")
    merged = eqtl.merge(mapped, left_on="variant", right_on="source_variant", how="inner")
    gwas_copy = gwas.copy()
    gwas_copy["lifted_chrom"] = gwas_copy[gwas_chrom_column].astype(str).str.removeprefix("chr")
    gwas_copy["lifted_position"] = pd.to_numeric(gwas_copy[gwas_position_column], errors="raise")
    return merged.merge(
        gwas_copy,
        on=["lifted_chrom", "lifted_position"],
        how="inner",
        suffixes=("", "_gwas"),
    )


def build_coordinate_coloc_table(
    eqtl: pd.DataFrame,
    lifted: pd.DataFrame,
    gwas: pd.DataFrame,
    *,
    gwas_chrom_column: str,
    gwas_position_column: str,
    gwas_effect_column: str,
    gwas_other_column: str,
    gwas_beta_column: str,
    gwas_se_column: str,
    gwas_n_column: str,
) -> pd.DataFrame:
    """Build an allele-aligned regional table with the coloc input contract."""
    merged = join_by_lifted_position(
        eqtl, lifted, gwas, gwas_chrom_column, gwas_position_column
    )
    def gwas_column(column: str) -> str:
        return f"{column}_gwas" if column in eqtl.columns else column

    effect_column = gwas_column(gwas_effect_column)
    other_column = gwas_column(gwas_other_column)
    beta_column = gwas_column(gwas_beta_column)
    se_column = gwas_column(gwas_se_column)
    n_column = gwas_column(gwas_n_column)
    alignments = [
        align_outcome_to_exposure(
            row.ref,
            row.alt,
            getattr(row, effect_column),
            getattr(row, other_column),
            getattr(row, beta_column),
        )
        for row in merged.itertuples(index=False)
    ]
    merged["alignment_status"] = [result.status for result in alignments]
    merged["gwas_beta_aligned"] = [result.beta for result in alignments]
    merged["se_eqtl"] = pd.to_numeric(merged["se"], errors="coerce")
    merged["se_gwas"] = pd.to_numeric(merged[se_column], errors="coerce")
    merged["N"] = pd.to_numeric(merged[n_column], errors="coerce")
    merged["BP"] = pd.to_numeric(merged["lifted_position"], errors="raise")
    return merged


def build_coordinate_coloc_audit(
    eqtl: pd.DataFrame,
    lifted: pd.DataFrame,
    gwas: pd.DataFrame,
    *,
    gwas_chrom_column: str,
    gwas_position_column: str,
    gwas_effect_column: str,
    gwas_other_column: str,
    gwas_beta_column: str,
    gwas_se_column: str,
    gwas_n_column: str,
) -> pd.DataFrame:
    """Return one auditable record for every eQTL row, including failures.

    This intentionally retains missing and non-mapped liftover rows, positions
    absent from the outcome source, allele mismatches and palindromic pairs.
    It does not infer a strand flip or normalize an indel.
    """
    audit = eqtl.merge(lifted, left_on="variant", right_on="source_variant", how="left")
    audit["alignment_status"] = "liftover_record_missing"
    has_liftover = audit["status"].notna()
    audit.loc[has_liftover, "alignment_status"] = (
        "liftover_" + audit.loc[has_liftover, "status"].astype(str)
    )
    mapped = audit["status"].eq("mapped") & audit["target_variant"].notna()
    components = audit.loc[mapped, "target_variant"].str.split("_", n=3, expand=True)
    audit.loc[mapped, "lifted_chrom"] = components[0].str.removeprefix("chr").to_numpy()
    audit.loc[mapped, "lifted_position"] = pd.to_numeric(components[1], errors="raise").to_numpy()

    gwas_copy = gwas.copy()
    gwas_copy["lifted_chrom"] = gwas_copy[gwas_chrom_column].astype(str).str.removeprefix("chr")
    gwas_copy["lifted_position"] = pd.to_numeric(gwas_copy[gwas_position_column], errors="raise")
    mapped_rows = audit.loc[mapped].merge(
        gwas_copy,
        on=["lifted_chrom", "lifted_position"],
        how="left",
        suffixes=("", "_gwas"),
    )
    if mapped_rows.empty:
        return audit

    def gwas_column(column: str) -> str:
        return f"{column}_gwas" if column in eqtl.columns else column

    effect_column = gwas_column(gwas_effect_column)
    other_column = gwas_column(gwas_other_column)
    beta_column = gwas_column(gwas_beta_column)
    se_column = gwas_column(gwas_se_column)
    n_column = gwas_column(gwas_n_column)
    has_outcome = mapped_rows[gwas_chrom_column].notna()
    mapped_rows["alignment_status"] = "no_outcome_at_lifted_coordinate"
    alignments = [
        align_outcome_to_exposure(
            row.ref, row.alt, getattr(row, effect_column), getattr(row, other_column),
            getattr(row, beta_column),
        )
        if has else AlignmentResult("no_outcome_at_lifted_coordinate", None)
        for row, has in zip(mapped_rows.itertuples(index=False), has_outcome, strict=True)
    ]
    mapped_rows["alignment_status"] = [result.status for result in alignments]
    mapped_rows["gwas_beta_aligned"] = [result.beta for result in alignments]
    mapped_rows["se_eqtl"] = pd.to_numeric(mapped_rows["se"], errors="coerce")
    mapped_rows["se_gwas"] = pd.to_numeric(mapped_rows[se_column], errors="coerce")
    mapped_rows["N"] = pd.to_numeric(mapped_rows[n_column], errors="coerce")
    mapped_rows["BP"] = pd.to_numeric(mapped_rows["lifted_position"], errors="coerce")
    return pd.concat([audit.loc[~mapped], mapped_rows], ignore_index=True, sort=False)
