"""Classical MR baselines used to benchmark hcsMR.

For each (gene, state) we extract the IVs assigned to that state and apply:
  - IVW (Burgess 2013)
  - MR-Egger (Bowden 2015): intercept + slope
  - Weighted median (Bowden 2016)
  - Weighted mode (Hartwig 2017) — simplified
  - "Bulk-IVW": pool all states' IVs as if bulk, ignoring cell-state labels

All baselines operate on (beta_hat[j], gamma_hat[j], se_X[j], se_Y[j]) at the
state level — i.e. the per-state slice of the simulated arrays.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm


def _weights_ivw(se_X: np.ndarray, se_Y: np.ndarray, ratio: np.ndarray) -> np.ndarray:
    """Delta-method weights for IVW: w_j ∝ 1 / (se_Y^2 / beta_hat^2)."""
    se_ratio = np.sqrt(se_Y ** 2 / np.maximum(np.abs(se_X) ** 2, 1e-12))
    # safer: use Wald ratio se from delta-method approx using se_Y / beta_hat
    return 1.0 / np.maximum(se_ratio ** 2, 1e-12)


def ivw_estimator(beta_hat: np.ndarray, gamma_hat: np.ndarray, se_X: np.ndarray, se_Y: np.ndarray):
    mask = np.isfinite(beta_hat) & np.isfinite(gamma_hat) & (se_X > 0) & (se_Y > 0)
    bh, gh, sx, sy = beta_hat[mask], gamma_hat[mask], se_X[mask], se_Y[mask]
    if bh.size < 3:
        return {"theta_hat": np.nan, "se": np.nan, "p": np.nan, "n_iv": int(bh.size)}
    ratio = gh / bh
    # delta-method se per IV
    se_ratio = np.sqrt(sy ** 2 / np.maximum(bh ** 2, 1e-12))
    w = 1.0 / np.maximum(se_ratio ** 2, 1e-12)
    theta = np.sum(w * ratio) / np.sum(w)
    se = 1.0 / np.sqrt(np.sum(w))
    z = theta / se if se > 0 else 0.0
    p = 2 * (1 - norm.cdf(abs(z)))
    return {"theta_hat": float(theta), "se": float(se), "p": float(p), "n_iv": int(bh.size)}


def mr_egger(beta_hat: np.ndarray, gamma_hat: np.ndarray, se_X: np.ndarray, se_Y: np.ndarray):
    mask = np.isfinite(beta_hat) & np.isfinite(gamma_hat) & (se_X > 0) & (se_Y > 0)
    bh, gh, sx, sy = beta_hat[mask], gamma_hat[mask], se_X[mask], se_Y[mask]
    if bh.size < 4:
        return {"theta_hat": np.nan, "se": np.nan, "p": np.nan, "intercept": np.nan, "intercept_p": np.nan, "n_iv": int(bh.size)}
    # Egger regression with weights w = 1 / sy^2
    w = 1.0 / np.maximum(sy ** 2, 1e-12)
    sw = w.sum()
    x_bar = np.sum(w * bh) / sw
    y_bar = np.sum(w * gh) / sw
    Sxx = np.sum(w * (bh - x_bar) ** 2)
    Sxy = np.sum(w * (bh - x_bar) * (gh - y_bar))
    beta_slope = Sxy / Sxx if Sxx > 0 else np.nan
    alpha_int = y_bar - beta_slope * x_bar
    # residual variance
    resid = gh - alpha_int - beta_slope * bh
    rss = np.sum(w * resid ** 2)
    dof = max(bh.size - 2, 1)
    sigma2 = rss / dof
    se_slope = np.sqrt(sigma2 / Sxx) if Sxx > 0 else np.nan
    se_alpha = np.sqrt(sigma2 * (1.0 / sw + x_bar ** 2 / Sxx)) if Sxx > 0 else np.nan
    z_slope = beta_slope / se_slope if se_slope > 0 else 0.0
    z_alpha = alpha_int / se_alpha if se_alpha > 0 else 0.0
    return {
        "theta_hat": float(beta_slope),
        "se": float(se_slope),
        "p": float(2 * (1 - norm.cdf(abs(z_slope)))),
        "intercept": float(alpha_int),
        "intercept_p": float(2 * (1 - norm.cdf(abs(z_alpha)))),
        "n_iv": int(bh.size),
    }


def weighted_median(beta_hat: np.ndarray, gamma_hat: np.ndarray, se_X: np.ndarray, se_Y: np.ndarray):
    mask = np.isfinite(beta_hat) & np.isfinite(gamma_hat) & (se_X > 0) & (se_Y > 0)
    bh, gh, sx, sy = beta_hat[mask], gamma_hat[mask], se_X[mask], se_Y[mask]
    if bh.size < 3:
        return {"theta_hat": np.nan, "p": np.nan, "n_iv": int(bh.size)}
    ratio = gh / bh
    se_ratio = np.sqrt(sy ** 2 / np.maximum(bh ** 2, 1e-12))
    w = 1.0 / np.maximum(se_ratio ** 2, 1e-12)
    w = w / w.sum()
    order = np.argsort(ratio)
    cw = np.cumsum(w[order])
    # weighted median: smallest i with cw[i] >= 0.5
    idx = np.searchsorted(cw, 0.5)
    theta = float(ratio[order][min(idx, len(ratio) - 1)])
    # bootstrap se
    rng = np.random.default_rng(0)
    boots = []
    for _ in range(500):
        b_bh = bh + rng.normal(0.0, sx)
        b_gh = gh + rng.normal(0.0, sy)
        b_ratio = b_gh / b_bh
        order_b = np.argsort(b_ratio)
        cw_b = np.cumsum(w[order_b])
        idx_b = np.searchsorted(cw_b, 0.5)
        boots.append(b_ratio[order_b][min(idx_b, len(b_ratio) - 1)])
    se = float(np.std(boots))
    z = theta / se if se > 0 else 0.0
    p = 2 * (1 - norm.cdf(abs(z)))
    return {"theta_hat": theta, "se": se, "p": float(p), "n_iv": int(bh.size)}


def run_state_baselines(beta_hat_per_state, se_X_per_state, gamma_hat, se_Y, iv_states):
    """For each cell state, gather its IVs and run the three classical baselines.

    beta_hat_per_state, se_X_per_state are (J, C) arrays; iv_states is (J,) telling
    which state each IV belongs to (we use this to pick which column to read).
    """
    J, C = beta_hat_per_state.shape
    results = {"ivw": [None] * C, "egger": [None] * C, "wmedian": [None] * C}
    for c in range(C):
        idx = np.where(iv_states == c)[0]
        if idx.size < 3:
            for k in results:
                results[k][c] = {"theta_hat": np.nan, "se": np.nan, "p": np.nan, "n_iv": int(idx.size)}
            continue
        bh = beta_hat_per_state[idx, c]
        sx = se_X_per_state[idx, c]
        gh = gamma_hat[idx]
        sy = se_Y[idx]
        results["ivw"][c] = ivw_estimator(bh, gh, sx, sy)
        results["egger"][c] = mr_egger(bh, gh, sx, sy)
        results["wmedian"][c] = weighted_median(bh, gh, sx, sy)
    return results


def bulk_ivw(beta_hat_per_state, se_X_per_state, gamma_hat, se_Y, iv_states, state_to_type, target_type=None):
    """`Bulk` IVW: pretend bulk eQTL exists.  For each cell type k, average the
    per-state betas across the states of that type (weighted by 1/se^2) to
    derive a synthetic bulk beta_hat, then run IVW on that pooled signal.
    Effectively this simulates what a researcher with only GTEx-tissue eQTL
    would obtain.
    """
    J, C = beta_hat_per_state.shape
    out = {}
    types = np.unique(state_to_type) if target_type is None else [target_type]
    for k in types:
        # gather all IVs whose assigned state belongs to type k
        idx = np.where(np.isin(iv_states, np.where(state_to_type == k)[0]))[0]
        if idx.size < 3:
            out[int(k)] = {"theta_hat": np.nan, "n_iv": int(idx.size)}
            continue
        # synthetic bulk beta_hat: weighted average across states of type k
        bulk_bh = []
        bulk_sx = []
        for j in idx:
            c = iv_states[j]
            # the bulk effect is the within-type weighted average — emulate dilution by
            # averaging beta_hat[j, c] (the per-state effect) against beta_hat[j, c'] for
            # other states c' of the same type which equal zero in our simulation.
            states_k = np.where(state_to_type == k)[0]
            vals = beta_hat_per_state[j, states_k]
            ses = se_X_per_state[j, states_k]
            mask = np.isfinite(vals) & np.isfinite(ses)
            if mask.sum() == 0:
                bulk_bh.append(np.nan)
                bulk_sx.append(np.nan)
                continue
            w = 1.0 / np.maximum(ses[mask] ** 2, 1e-12)
            bulk_bh.append(np.sum(w * vals[mask]) / np.sum(w))
            bulk_sx.append(np.sqrt(1.0 / np.sum(w)))
        bulk_bh = np.asarray(bulk_bh)
        bulk_sx = np.asarray(bulk_sx)
        gh = gamma_hat[idx]
        sy = se_Y[idx]
        out[int(k)] = ivw_estimator(bulk_bh, gh, bulk_sx, sy)
    return out
