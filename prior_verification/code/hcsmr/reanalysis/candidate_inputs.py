"""Build traceable diagnostic inputs from exploratory candidate tables."""
from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from .diagnostics import with_f_statistics


def build_candidate_variants(
    hits: pd.DataFrame, outcome_tables: Mapping[str, pd.DataFrame]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return variant rows and one explicit status record per exploratory candidate."""
    required = {
        "outcome",
        "gene_id",
        "cell_state",
        "is_robust_corrected",
        "passes_bonferroni_corrected",
    }
    missing = required.difference(hits.columns)
    if missing:
        raise ValueError(f"hits missing required columns: {sorted(missing)}")
    selected = hits.loc[
        hits["is_robust_corrected"].eq(1)
        & hits["passes_bonferroni_corrected"].eq(1),
        ["outcome", "gene_id", "cell_state"],
    ].drop_duplicates()
    variant_frames: list[pd.DataFrame] = []
    status_rows: list[dict[str, object]] = []
    for row in selected.itertuples(index=False):
        candidate_id = f"{row.outcome}|{row.gene_id}|{row.cell_state}"
        table = outcome_tables.get(row.outcome)
        if table is None:
            status_rows.append({"candidate_id": candidate_id, "status": "missing_outcome_table"})
            continue
        subset = table.loc[
            table["gene_id"].eq(row.gene_id) & table["state"].eq(row.cell_state)
        ].copy()
        if subset.empty:
            status_rows.append({"candidate_id": candidate_id, "status": "missing_variant_rows"})
            continue
        subset["candidate_id"] = candidate_id
        variant_frames.append(with_f_statistics(subset))
        status_rows.append({"candidate_id": candidate_id, "status": "ready_for_diagnostics"})
    variants = pd.concat(variant_frames, ignore_index=True) if variant_frames else pd.DataFrame()
    return variants, pd.DataFrame(status_rows, columns=["candidate_id", "status"])
