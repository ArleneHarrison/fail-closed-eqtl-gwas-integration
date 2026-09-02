"""Scenario-level validation for no, complete, and partial pooling.

The validation intentionally uses a small, biologically interpretable
two-state system.  It distinguishes the state-specific no-pooling estimand
from the all-signal complete-pooling estimand, so a pooled result is never
misreported as a cell-state effect.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import erfc, sqrt
from typing import Iterable

import numpy as np
import pandas as pd

from hcsmr.baselines import ivw_estimator
from hcsmr.simulate import SimulationScenario, simulate_scenario


@dataclass(frozen=True)
class HierarchicalValidationScenario:
    name: str
    ivs_per_state: int
    causal_theta: float
    pleiotropy_fraction: float
    cell_state_freq: float = 0.1


def required_hierarchical_validation_scenarios() -> tuple[HierarchicalValidationScenario, ...]:
    """Return the locked validation grid used by the executable report."""
    return (
        HierarchicalValidationScenario("null_sparse_clean", 3, 0.0, 0.0),
        HierarchicalValidationScenario("signal_sparse_clean", 3, 0.15, 0.0),
        HierarchicalValidationScenario("signal_sparse_pleiotropy", 3, 0.15, 0.15),
        HierarchicalValidationScenario("signal_dense_clean", 10, 0.15, 0.0),
        HierarchicalValidationScenario("signal_dense_pleiotropy", 10, 0.15, 0.15),
    )


def _p_value(estimate: float, standard_error: float) -> float:
    return erfc(abs(estimate / standard_error) / sqrt(2.0)) if standard_error > 0 else float("nan")


def _state_slice(simulated: dict, state: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mask = simulated["iv_in_state"][:, state]
    return (
        simulated["beta_hat"][mask, state],
        simulated["gamma_hat"][mask],
        simulated["se_X"][mask, state],
        simulated["se_Y"][mask],
    )


def _pooled_slice(simulated: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mask = simulated["iv_in_state"]
    row, state = np.where(mask)
    return (
        simulated["beta_hat"][row, state],
        simulated["gamma_hat"][row],
        simulated["se_X"][row, state],
        simulated["se_Y"][row],
    )


def _classical_row(method: str, arrays: tuple[np.ndarray, ...], true_theta: float) -> dict[str, object]:
    estimate = ivw_estimator(*arrays)
    theta = float(estimate["theta_hat"])
    standard_error = float(estimate["se"])
    eligible = np.isfinite(theta) and np.isfinite(standard_error)
    return {
        "method": method,
        "eligible": eligible,
        "ineligibility_reason": "" if eligible else "fewer_than_three_instruments",
        "estimate": theta,
        "standard_error": standard_error,
        "ci_lower": theta - 1.959963984540054 * standard_error if eligible else float("nan"),
        "ci_upper": theta + 1.959963984540054 * standard_error if eligible else float("nan"),
        "p_value": _p_value(theta, standard_error) if eligible else float("nan"),
        "n_iv": int(estimate["n_iv"]),
        "n_divergences": float("nan"),
        "max_rhat": float("nan"),
        "posterior_predictive_tail_probability": float("nan"),
        "true_theta": true_theta,
    }


def run_validation_replicate(
    scenario: HierarchicalValidationScenario,
    *,
    replicate: int,
    warmup: int,
    samples: int,
    chains: int,
    target_accept_prob: float = 0.98,
) -> pd.DataFrame:
    """Fit all three pooling strategies to one simulated data set."""
    simulated = simulate_scenario(
        SimulationScenario(
            name=scenario.name,
            n_types=1,
            states_per_type=2,
            causal_states=(0,),
            causal_theta=scenario.causal_theta,
            ivs_per_state=scenario.ivs_per_state,
            F_stat_mean=30.0,
            pleiotropy_fraction=scenario.pleiotropy_fraction,
            cell_state_freq=scenario.cell_state_freq,
            rng_seed=202608150 + replicate,
        )
    )
    rows = [
        _classical_row("no_pooling_state_ivw", _state_slice(simulated, 0), scenario.causal_theta),
        _classical_row("complete_pooling_all_signal_ivw", _pooled_slice(simulated), scenario.causal_theta),
    ]
    try:
        from hcsmr.model import fit_hcsmr

        fitted = fit_hcsmr(
            beta_hat=simulated["beta_hat"],
            se_X=simulated["se_X"],
            gamma_hat=simulated["gamma_hat"],
            se_Y=simulated["se_Y"],
            state_to_type=simulated["state_to_type"],
            iv_in_state=simulated["iv_in_state"],
            n_warmup=warmup,
            n_samples=samples,
            n_chains=chains,
            target_accept_prob=target_accept_prob,
            seed=202608150 + replicate,
            progress_bar=False,
            chain_method="sequential",
        )
        theta = fitted["summary"]["theta"]
        diagnostics = fitted["diagnostics"]
        rows.append(
            {
                "method": "partial_pooling_hcsmr",
                "eligible": True,
                "ineligibility_reason": "",
                "estimate": float(theta["mean"][0]),
                "standard_error": float(theta["sd"][0]),
                "ci_lower": float(theta["q025"][0]),
                "ci_upper": float(theta["q975"][0]),
                "p_value": float("nan"),
                "n_iv": int(simulated["iv_in_state"][:, 0].sum()),
                "n_divergences": diagnostics["n_divergences"],
                "max_rhat": float(max(theta["rhat"])),
                "posterior_predictive_tail_probability": diagnostics["posterior_predictive"].get(
                    "tail_probability", float("nan")
                ),
                "true_theta": scenario.causal_theta,
            }
        )
    except Exception as exc:  # Preserve failed fits in the publication output.
        rows.append(
            {
                "method": "partial_pooling_hcsmr",
                "eligible": False,
                "ineligibility_reason": f"fit_failed:{type(exc).__name__}",
                "estimate": float("nan"),
                "standard_error": float("nan"),
                "ci_lower": float("nan"),
                "ci_upper": float("nan"),
                "p_value": float("nan"),
                "n_iv": int(simulated["iv_in_state"][:, 0].sum()),
                "n_divergences": float("nan"),
                "max_rhat": float("nan"),
                "posterior_predictive_tail_probability": float("nan"),
                "true_theta": scenario.causal_theta,
            }
        )
    for row in rows:
        row.update(
            {
                "scenario": scenario.name,
                "replicate": replicate,
                "ivs_per_state": scenario.ivs_per_state,
                "pleiotropy_fraction": scenario.pleiotropy_fraction,
                "covered": bool(
                    row["eligible"]
                    and row["ci_lower"] <= scenario.causal_theta <= row["ci_upper"]
                ),
            }
        )
    return pd.DataFrame(rows)


def summarize_hierarchical_validation(rows: pd.DataFrame) -> pd.DataFrame:
    """Summarize precision, calibration, and fitting failures without dropping them."""
    output = []
    for (scenario, method), group in rows.groupby(["scenario", "method"], sort=False):
        eligible = group.loc[group["eligible"]]
        true_theta = float(group["true_theta"].iloc[0])
        finite_divergences = eligible["n_divergences"].dropna()
        finite_rhat = eligible["max_rhat"].dropna()
        finite_ppc = eligible["posterior_predictive_tail_probability"].dropna()
        output.append(
            {
                "scenario": scenario,
                "method": method,
                "n_replicates": int(len(group)),
                "eligible_replicates": int(len(eligible)),
                "failed_or_ineligible_replicates": int(len(group) - len(eligible)),
                "mean_estimate": float(eligible["estimate"].mean()) if len(eligible) else float("nan"),
                "bias_vs_state0": float(eligible["estimate"].mean() - true_theta) if len(eligible) else float("nan"),
                "coverage": float(eligible["covered"].mean()) if len(eligible) else float("nan"),
                "mean_n_divergences": float(finite_divergences.mean()) if len(finite_divergences) else float("nan"),
                "max_rhat": float(finite_rhat.max()) if len(finite_rhat) else float("nan"),
                "median_posterior_predictive_tail_probability": float(finite_ppc.median()) if len(finite_ppc) else float("nan"),
            }
        )
    return pd.DataFrame(output)
