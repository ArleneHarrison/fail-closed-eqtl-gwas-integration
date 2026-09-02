"""Prepare identically ordered regional summaries and LD for coloc-SuSiE."""
from __future__ import annotations

import numpy as np
import pandas as pd


def align_harmonized_to_ld(
    harmonized: pd.DataFrame,
    ld_variants: list[str],
    ld: np.ndarray,
    molecular_trait_id: str | None = None,
) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    if ld.shape != (len(ld_variants), len(ld_variants)):
        raise ValueError("LD dimensions do not match its variant identifier list")
    if molecular_trait_id is not None:
        if "molecular_trait_id" not in harmonized.columns:
            raise ValueError("molecular_trait_id selection requested but column is absent")
        harmonized = harmonized.loc[
            harmonized["molecular_trait_id"].astype(str) == str(molecular_trait_id)
        ].copy()
        if harmonized.empty:
            raise ValueError(f"no rows found for molecular_trait_id={molecular_trait_id}")
    unique = harmonized.drop_duplicates("target_variant", keep=False)
    by_variant = unique.set_index("target_variant", drop=False)
    retained = [variant for variant in ld_variants if variant in by_variant.index]
    if len(retained) < 2:
        raise ValueError("fewer than two harmonized variants overlap the LD matrix")
    index = {variant: position for position, variant in enumerate(ld_variants)}
    positions = [index[variant] for variant in retained]
    dropped = sorted(set(harmonized["target_variant"]).difference(retained))
    return by_variant.loc[retained].reset_index(drop=True), ld[np.ix_(positions, positions)], dropped
