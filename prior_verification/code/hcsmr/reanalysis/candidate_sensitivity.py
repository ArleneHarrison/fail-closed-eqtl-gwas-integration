"""Transparent single-context sensitivity diagnostics for candidate records.

These functions deliberately do not turn a small set of fine-mapped signals
into a confirmed MR finding.  They report what can be calculated from the
available independent-signal table and make unavailable diagnostics explicit.
"""
from __future__ import annotations

from math import erfc, sqrt

import numpy as np
import pandas as pd
from scipy.stats import chi2


REQUIRED_COLUMNS = {
    "candidate_id",
    "beta",
    "se",
    "gwas_beta",
    "gwas_se",
    "f_stat",
}


def _two_sided_p(z_score: float) -> float:
    return erfc(abs(z_score) / sqrt(2.0))


def _valid_instruments(variants: pd.DataFrame) -> pd.DataFrame:
    missing = REQUIRED_COLUMNS.difference(variants.columns)
    if missing:
        raise ValueError(f"variants missing required columns: {sorted(missing)}")
    frame = variants.copy()
    numeric = ["beta", "se", "gwas_beta", "gwas_se", "f_stat"]
    frame[numeric] = frame[numeric].apply(pd.to_numeric, errors="coerce")
    valid = (
        np.isfinite(frame["beta"])
        & np.isfinite(frame["se"])
        & np.isfinite(frame["gwas_beta"])
        & np.isfinite(frame["gwas_se"])
        & np.isfinite(frame["f_stat"])
        & (frame["beta"] != 0)
        & (frame["se"] > 0)
        & (frame["gwas_se"] > 0)
    )
    return frame.loc[valid].copy()


def _variant_label(row: pd.Series, fallback: object) -> str:
    """Prefer rsID, then the source variant key, without emitting literal NaN."""
    for field in ("rsid", "variant", "chrpos"):
        value = row.get(field)
        if pd.notna(value) and str(value).strip():
            return str(value)
    return str(fallback)


def ivw_summary(variants: pd.DataFrame) -> dict[str, float]:
    """Return fixed-effect IVW and Cochran's Q for a single candidate context."""
    frame = _valid_instruments(variants)
    if len(frame) < 2:
        raise ValueError("at least two valid instruments are required for IVW")
    beta = frame["beta"].to_numpy(float)
    gamma = frame["gwas_beta"].to_numpy(float)
    se_y = frame["gwas_se"].to_numpy(float)
    weights = 1.0 / np.square(se_y)
    denominator = float(np.sum(weights * np.square(beta)))
    estimate = float(np.sum(weights * beta * gamma) / denominator)
    standard_error = float(sqrt(1.0 / denominator))
    z_score = estimate / standard_error
    residual = gamma - estimate * beta
    q_statistic = float(np.sum(weights * np.square(residual)))
    q_df = int(len(frame) - 1)
    return {
        "n_iv": int(len(frame)),
        "estimate": estimate,
        "standard_error": standard_error,
        "p_value": _two_sided_p(z_score),
        "q_statistic": q_statistic,
        "q_df": q_df,
        "q_p_value": float(chi2.sf(q_statistic, q_df)),
    }


def candidate_sensitivity(variants: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Produce candidate-level and leave-one-out diagnostic tables.

    A leave-one-out stability decision requires at least four original signals,
    leaving three after omission.  Pairwise LD, overlap, and joint-posterior
    diagnostics are intentionally not inferred from a lead-signal table.
    """
    frame = _valid_instruments(variants)
    summaries: list[dict[str, object]] = []
    leave_one_out_rows: list[dict[str, object]] = []
    for candidate_id, group in frame.groupby("candidate_id", sort=False):
        primary = ivw_summary(group)
        n_iv = int(primary["n_iv"])
        strongest_index = group["f_stat"].idxmax()
        strongest_label = _variant_label(group.loc[strongest_index], strongest_index)
        if n_iv >= 4:
            for omitted_index, omitted in group.iterrows():
                reduced = ivw_summary(group.drop(index=omitted_index))
                leave_one_out_rows.append(
                    {
                        "candidate_id": candidate_id,
                        "omitted_variant": _variant_label(omitted, omitted_index),
                        **reduced,
                    }
                )
            leave_one_out = pd.DataFrame(
                [row for row in leave_one_out_rows if row["candidate_id"] == candidate_id]
            )
            stable_sign = bool(np.all(np.sign(leave_one_out["estimate"]) == np.sign(primary["estimate"])))
            within_interval = np.abs(leave_one_out["estimate"] - primary["estimate"]) <= (
                1.959963984540054
                * np.sqrt(np.square(leave_one_out["standard_error"]) + primary["standard_error"] ** 2)
            )
            leave_one_out_status = "stable" if stable_sign and bool(np.all(within_interval)) else "unstable"
            strongest = leave_one_out.loc[
                leave_one_out["omitted_variant"].eq(strongest_label)
            ].iloc[0]
            strongest_status = "stable" if (
                np.sign(strongest["estimate"]) == np.sign(primary["estimate"])
                and abs(strongest["estimate"] - primary["estimate"])
                <= 1.959963984540054
                * sqrt(strongest["standard_error"] ** 2 + primary["standard_error"] ** 2)
            ) else "unstable"
        else:
            leave_one_out_status = "not_assessed_fewer_than_four_signals"
            strongest_status = "not_assessed_fewer_than_four_signals"
        summaries.append(
            {
                "candidate_id": candidate_id,
                **primary,
                "min_f_stat": float(group["f_stat"].min()),
                "mean_f_stat": float(group["f_stat"].mean()),
                "instrument_strength_status": "pass" if bool((group["f_stat"] >= 10.0).all()) else "fail",
                "leave_one_out_status": leave_one_out_status,
                "strongest_iv_status": strongest_status,
                "ld_threshold_sensitivity": "not_identifiable_without_pairwise_ld",
                "sample_overlap_status": "unknown_without_participant_linkage",
                "joint_model_status": "not_fit_noncomparable_or_covariance_unavailable",
            }
        )
    return pd.DataFrame(summaries), pd.DataFrame(leave_one_out_rows)
