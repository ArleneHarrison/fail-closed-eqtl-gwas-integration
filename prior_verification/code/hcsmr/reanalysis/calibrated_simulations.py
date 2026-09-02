"""Calibrated stress-test simulations for transparent comparator reporting."""
from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

import numpy as np
import pandas as pd

from .benchmarks import run_comparator_suite
from .simulation_protocol import SimulationScenarioV2, binomial_interval


@dataclass(frozen=True)
class ComparatorData:
    beta_hat: np.ndarray
    gamma_hat: np.ndarray
    se_x: np.ndarray
    se_y: np.ndarray

    @property
    def n_instruments(self) -> int:
        return len(self.beta_hat)


def _draw_summary_data(
    scenario: SimulationScenarioV2, true_theta: float, seed: int
) -> ComparatorData:
    rng = np.random.default_rng(seed)
    n = scenario.n_instruments
    # Sparse settings are deliberately weaker, rather than silently borrowing
    # information from a dense-instrument regime.
    median_f = 8.0 if n <= 5 else 25.0
    f_stat = np.clip(rng.lognormal(np.log(median_f), 0.45, n), 2.0, 100.0)
    if scenario.ld == "residual_ld":
        innovations = rng.normal(0.0, 0.025, n)
        beta_true = np.empty(n)
        beta_true[0] = 0.12 + innovations[0]
        for index in range(1, n):
            beta_true[index] = 0.8 * beta_true[index - 1] + 0.024 + innovations[index]
    else:
        beta_true = rng.normal(0.12, 0.035, n)
    beta_true = np.clip(beta_true, 0.02, None)
    se_x = beta_true / np.sqrt(f_stat)
    se_y = np.full(n, 0.025)
    correlation = scenario.sample_overlap
    errors_x = rng.normal(0.0, se_x)
    errors_y = correlation * (errors_x / se_x) * se_y + sqrt(1 - correlation**2) * rng.normal(0.0, se_y)
    if scenario.pleiotropy == "balanced":
        alpha = rng.normal(0.0, 0.025, n)
    elif scenario.pleiotropy == "directional":
        alpha = np.abs(rng.normal(0.02, 0.02, n))
    elif scenario.pleiotropy == "correlated":
        alpha = 0.18 * beta_true + rng.normal(0.0, 0.012, n)
    elif scenario.pleiotropy == "inside_violation":
        alpha = 0.05 * (beta_true - beta_true.mean()) / beta_true.std() + rng.normal(0.0, 0.012, n)
    else:
        raise ValueError(f"unknown pleiotropy setting: {scenario.pleiotropy}")
    return ComparatorData(
        beta_hat=beta_true + errors_x,
        gamma_hat=true_theta * beta_true + alpha + errors_y,
        se_x=se_x,
        se_y=se_y,
    )


def run_calibrated_replicates(
    scenario: SimulationScenarioV2,
    *,
    replicates: int,
    true_theta: float,
    seed_offset: int = 0,
) -> pd.DataFrame:
    """Run a reproducible comparator stress test without discarding failures."""
    if replicates <= 0:
        raise ValueError("replicates must be positive")
    rows = []
    for replicate in range(replicates):
        data = _draw_summary_data(scenario, true_theta, seed_offset + replicate)
        methods = run_comparator_suite(data, ("ivw", "egger", "weighted_median"))
        for result in methods.to_dict("records"):
            eligible = bool(result["eligible"])
            p_value = float(result["p_value"])
            estimate = float(result["estimate"])
            se = float(result["standard_error"])
            rows.append(
                {
                    "replicate": replicate,
                    "n_instruments": scenario.n_instruments,
                    "pleiotropy": scenario.pleiotropy,
                    "ld": scenario.ld,
                    "sample_overlap": scenario.sample_overlap,
                    "true_theta": true_theta,
                    **result,
                    "positive": eligible and np.isfinite(p_value) and p_value < 0.05,
                    "covered": eligible and np.isfinite(se) and abs(estimate - true_theta) <= 1.959963984540054 * se,
                }
            )
    return pd.DataFrame(rows)


def summarize_method_performance(results: pd.DataFrame, *, true_theta: float) -> pd.DataFrame:
    """Summarize calibration/power with Wilson intervals and all failures shown."""
    rows = []
    for method, group in results.groupby("method", sort=False):
        eligible = group.loc[group["eligible"]]
        n = len(group)
        eligible_n = len(eligible)
        positives = int(eligible["positive"].sum())
        covered = int(eligible["covered"].sum())
        positive_low, positive_high = binomial_interval(positives, eligible_n) if eligible_n else (np.nan, np.nan)
        coverage_low, coverage_high = binomial_interval(covered, eligible_n) if eligible_n else (np.nan, np.nan)
        rows.append(
            {
                "method": method,
                "n_replicates": n,
                "eligible_replicates": eligible_n,
                "ineligible_replicates": n - eligible_n,
                "mean_estimate": float(eligible["estimate"].mean()) if eligible_n else np.nan,
                "bias": float(eligible["estimate"].mean() - true_theta) if eligible_n else np.nan,
                "positive_rate": positives / eligible_n if eligible_n else np.nan,
                "type1_rate" if true_theta == 0 else "power": positives / eligible_n if eligible_n else np.nan,
                "type1_ci_low" if true_theta == 0 else "power_ci_low": positive_low,
                "type1_ci_high" if true_theta == 0 else "power_ci_high": positive_high,
                "coverage": covered / eligible_n if eligible_n else np.nan,
                "coverage_ci_low": coverage_low,
                "coverage_ci_high": coverage_high,
            }
        )
    return pd.DataFrame(rows)
