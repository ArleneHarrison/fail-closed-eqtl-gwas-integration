"""Eligibility rules for joint context-level regional analyses."""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class RegionalEligibility:
    joint_fit_allowed: bool
    n_shared_variants: int
    reason: str


def assess_regional_eligibility(
    variants: pd.DataFrame,
    *,
    min_shared_variants: int = 3,
    min_f_stat: float = 10.0,
) -> RegionalEligibility:
    """Assess whether contexts share enough strong variants for joint fitting.

    The function records the reason a region cannot be jointly analysed. It
    intentionally does not turn separate context-wise estimates into a
    localization claim.
    """
    required = {"context", "variant_id", "f_stat"}
    missing = required.difference(variants.columns)
    if missing:
        raise ValueError(f"variants missing required columns: {sorted(missing)}")
    if variants.empty:
        return RegionalEligibility(False, 0, "no_variants")

    by_context = variants.groupby("context", sort=False)["variant_id"].agg(set)
    variant_sets = list(by_context)
    shared = set.intersection(*variant_sets) if len(variant_sets) > 1 else variant_sets[0]
    weak_contexts = variants.groupby("context", sort=False)["f_stat"].min().lt(min_f_stat)

    reasons: list[str] = []
    if len(shared) < min_shared_variants:
        reasons.append("insufficient_shared_variants")
    if weak_contexts.any():
        reasons.append("weak_instruments")
    return RegionalEligibility(
        joint_fit_allowed=not reasons,
        n_shared_variants=len(shared),
        reason=";".join(reasons) if reasons else "eligible",
    )
