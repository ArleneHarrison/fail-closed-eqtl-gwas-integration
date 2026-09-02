"""Schema-stable diagnostics for real-data candidate associations."""
from __future__ import annotations

from typing import Any

import pandas as pd


def with_f_statistics(variants: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with conventional single-variant F statistics."""
    required = {"beta", "se"}
    missing = required.difference(variants.columns)
    if missing:
        raise ValueError(f"variants missing required columns: {sorted(missing)}")
    result = variants.copy()
    beta = pd.to_numeric(result["beta"], errors="coerce")
    standard_error = pd.to_numeric(result["se"], errors="coerce")
    result["f_stat"] = (beta / standard_error) ** 2
    return result


def diagnostic_row(candidate_id: str, variants: pd.DataFrame) -> dict[str, Any]:
    """Build a diagnostic record without converting unavailable checks to passes."""
    if "f_stat" not in variants.columns:
        raise ValueError("variants must contain f_stat")
    f_stats = pd.to_numeric(variants["f_stat"], errors="coerce").dropna()
    if f_stats.empty:
        raise ValueError("variants contain no valid F statistics")
    return {
        "candidate_id": candidate_id,
        "n_iv": int(len(variants)),
        "min_f_stat": float(f_stats.min()),
        "mean_f_stat": float(f_stats.mean()),
        "instrument_pass": bool((f_stats >= 10.0).all()),
        "leave_one_out_pass": "not_assessed",
        "strongest_iv_pass": "not_assessed",
        "heterogeneity_status": "not_assessed",
        "ld_threshold_sensitivity": "not_assessed",
        "posterior_predictive_status": "not_assessed",
        "convergence_pass": "not_assessed",
        "overlap_status": "not_assessed",
    }
