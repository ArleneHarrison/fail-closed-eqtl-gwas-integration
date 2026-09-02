"""Streaming-safe filters for complete regional eQTL and GWAS summaries."""
from __future__ import annotations

import pandas as pd


def filter_eqtl_region(
    frame: pd.DataFrame, gene_id: str, chrom: str, start: int, end: int
) -> pd.DataFrame:
    required = {"gene_id", "chromosome", "position"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"eQTL frame missing columns: {sorted(missing)}")
    chromosome = frame["chromosome"].astype(str).str.removeprefix("chr")
    return frame.loc[
        (frame["gene_id"] == gene_id)
        & (chromosome == str(chrom).removeprefix("chr"))
        & frame["position"].between(start, end)
    ].copy()


def filter_gwas_region(
    frame: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    chrom_column: str,
    position_column: str,
) -> pd.DataFrame:
    required = {chrom_column, position_column}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"GWAS frame missing columns: {sorted(missing)}")
    chromosome = frame[chrom_column].astype(str).str.removeprefix("chr")
    return frame.loc[
        (chromosome == str(chrom).removeprefix("chr"))
        & pd.to_numeric(frame[position_column], errors="coerce").between(start, end)
    ].copy()
