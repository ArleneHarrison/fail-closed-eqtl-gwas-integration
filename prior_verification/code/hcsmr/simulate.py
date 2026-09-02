"""Generate synthetic (sc-eQTL, GWAS) summary statistics for hcsMR validation.

We simulate a single-gene scenario in which:
  - There are K cell types, each with several cell states (M total states).
  - Each cell state c has a true cell-state causal effect theta_c on the
    outcome; many of these are 0, a few are non-zero (typically restricted to
    one or two states per gene).
  - SNPs serve as IVs for a single cell state each (extension: shared IVs).
  - Per-IV sc-eQTL effect beta_{j,c} has its own mean per state, drawn from a
    Normal centred on a within-type mu_beta_type.
  - Outcome effect gamma_j is generated as theta_{c(j)} * beta_{j,c(j)} +
    alpha_j  where alpha_j is non-zero with probability pi_alpha (horizontal
    pleiotropy).
  - Both beta_hat and gamma_hat have Gaussian sampling noise calibrated to a
    target F-statistic.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass
class SimulationScenario:
    name: str
    n_types: int = 5
    states_per_type: int = 3
    causal_states: Sequence[int] = (0,)      # which state indices have non-zero theta
    causal_theta: float = 0.15                # OR-per-SD style effect
    ivs_per_state: int = 20
    F_stat_mean: float = 30.0                 # target F-stat for typical IV
    pleiotropy_fraction: float = 0.0          # fraction of IVs with non-zero alpha
    alpha_sd: float = 0.05
    cell_state_freq: float = 0.1              # only affects sample-size proxy (we
                                              # inflate se for rare states)
    rng_seed: int = 0

    @property
    def n_states(self) -> int:
        return self.n_types * self.states_per_type

    @property
    def n_ivs(self) -> int:
        return self.n_states * self.ivs_per_state


def simulate_scenario(s: SimulationScenario):
    rng = np.random.default_rng(s.rng_seed)

    # state -> type map
    state_to_type = np.repeat(np.arange(s.n_types), s.states_per_type)

    # true theta_c: zero except at causal states
    theta_true = np.zeros(s.n_states, dtype=float)
    for c in s.causal_states:
        theta_true[c] = s.causal_theta

    # per-state beta means
    mu_beta_state = rng.normal(0.0, 0.5, size=s.n_states)

    # build IVs: each IV is assigned to one state
    iv_states = np.repeat(np.arange(s.n_states), s.ivs_per_state)  # (J,)
    J = iv_states.size

    # True per-(j, c) beta: noise around its assigned state's mean, near zero in
    # other states (we encode IV specificity by giving non-assigned states a tiny
    # beta plus se = NaN for the mask).
    beta_true = np.full((J, s.n_states), np.nan, dtype=float)
    se_X = np.full((J, s.n_states), np.nan, dtype=float)
    for j in range(J):
        c = iv_states[j]
        beta_true[j, c] = rng.normal(mu_beta_state[c], 0.2)
        # se calibrated to target F-stat:  F = (beta / se)^2  -> se = |beta| / sqrt(F)
        target_F = s.F_stat_mean * (1.0 + 0.5 * rng.standard_normal())
        target_F = max(5.0, float(target_F))
        se = max(abs(beta_true[j, c]), 0.02) / np.sqrt(target_F)
        # rare states get noisier se
        if s.cell_state_freq < 0.05:
            se *= 1.5
        if s.cell_state_freq < 0.01:
            se *= 2.0
        se_X[j, c] = se

    iv_in_state = ~np.isnan(beta_true)

    # observed beta_hat: true + Gaussian noise (filled only where iv_in_state)
    beta_hat = np.where(
        iv_in_state,
        beta_true + rng.normal(0.0, np.where(iv_in_state, np.nan_to_num(se_X, nan=0.1), 0.1), size=beta_true.shape),
        np.nan,
    )

    # outcome side
    gamma_true = np.zeros(J, dtype=float)
    for j in range(J):
        c = iv_states[j]
        gamma_true[j] = theta_true[c] * beta_true[j, c]

    # pleiotropy
    alpha_true = np.zeros(J, dtype=float)
    n_pleio = int(round(s.pleiotropy_fraction * J))
    if n_pleio > 0:
        pleio_idx = rng.choice(J, size=n_pleio, replace=False)
        alpha_true[pleio_idx] = rng.normal(0.0, s.alpha_sd, size=n_pleio)

    se_Y = np.full(J, 0.01, dtype=float)  # GWAS se ~ matched to large GWAS
    gamma_hat = gamma_true + alpha_true + rng.normal(0.0, se_Y)

    return {
        "scenario": s.name,
        "state_to_type": state_to_type,
        "theta_true": theta_true,
        "beta_true": beta_true,
        "alpha_true": alpha_true,
        "iv_states": iv_states,
        "beta_hat": beta_hat,
        "se_X": se_X,
        "gamma_hat": gamma_hat,
        "se_Y": se_Y,
        "iv_in_state": iv_in_state,
    }
