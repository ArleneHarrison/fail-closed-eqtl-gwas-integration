"""Comparator adapters that preserve ineligible-method records."""
from __future__ import annotations

from math import erfc, sqrt
from typing import Iterable

import numpy as np
import pandas as pd

def _ineligible(method: str, reason: str) -> dict[str, object]:
    return {
        "method": method,
        "estimate": float("nan"),
        "standard_error": float("nan"),
        "p_value": float("nan"),
        "eligible": False,
        "ineligibility_reason": reason,
    }


def _p_value(z_score: float) -> float:
    return erfc(abs(z_score) / sqrt(2.0))


def _valid_arrays(data):
    mask = (
        np.isfinite(data.beta_hat)
        & np.isfinite(data.gamma_hat)
        & np.isfinite(data.se_x)
        & np.isfinite(data.se_y)
        & (data.se_x > 0)
        & (data.se_y > 0)
    )
    return data.beta_hat[mask], data.gamma_hat[mask], data.se_x[mask], data.se_y[mask]


def _ivw(data) -> dict[str, float]:
    beta, gamma, _, se_y = _valid_arrays(data)
    weights = 1.0 / np.maximum(se_y**2, 1e-12)
    denominator = float(np.sum(weights * beta**2))
    estimate = float(np.sum(weights * beta * gamma) / denominator)
    standard_error = float(sqrt(1.0 / denominator))
    return {"theta_hat": estimate, "se": standard_error, "p": _p_value(estimate / standard_error)}


def _egger(data) -> dict[str, float]:
    beta, gamma, _, se_y = _valid_arrays(data)
    weights = 1.0 / np.maximum(se_y**2, 1e-12)
    total_weight = float(weights.sum())
    beta_mean = float(np.sum(weights * beta) / total_weight)
    gamma_mean = float(np.sum(weights * gamma) / total_weight)
    sxx = float(np.sum(weights * (beta - beta_mean) ** 2))
    slope = float(np.sum(weights * (beta - beta_mean) * (gamma - gamma_mean)) / sxx)
    intercept = gamma_mean - slope * beta_mean
    residual = gamma - intercept - slope * beta
    scale = float(np.sum(weights * residual**2) / max(len(beta) - 2, 1))
    standard_error = float(sqrt(scale / sxx))
    return {"theta_hat": slope, "se": standard_error, "p": _p_value(slope / standard_error)}


def _weighted_median(data) -> dict[str, float]:
    beta, gamma, _, se_y = _valid_arrays(data)
    ratios = gamma / beta
    weights = (beta**2) / np.maximum(se_y**2, 1e-12)
    order = np.argsort(ratios)
    ratios, weights = ratios[order], weights[order]
    estimate = float(ratios[np.searchsorted(np.cumsum(weights) / weights.sum(), 0.5)])
    # The bootstrap is deterministic for reproducible reporting.
    generator = np.random.default_rng(0)
    bootstrap_indices = generator.integers(0, len(ratios), size=(500, len(ratios)))
    bootstrap = np.median(ratios[bootstrap_indices], axis=1)
    standard_error = float(np.std(bootstrap, ddof=1))
    p_value = _p_value(estimate / standard_error) if standard_error > 0 else float("nan")
    return {"theta_hat": estimate, "se": standard_error, "p": p_value}


def _supported(method: str, data) -> dict[str, object]:
    if method == "ivw":
        result = _ivw(data)
    elif method == "egger":
        result = _egger(data)
    elif method == "weighted_median":
        result = _weighted_median(data)
    else:
        return _ineligible(method, "method_not_implemented")
    return {
        "method": method,
        "estimate": result.get("theta_hat", float("nan")),
        "standard_error": result.get("se", float("nan")),
        "p_value": result.get("p", float("nan")),
        "eligible": True,
        "ineligibility_reason": "",
    }


def run_comparator_suite(data, methods: Iterable[str]) -> pd.DataFrame:
    """Run supported methods and retain explicit records for ineligible methods."""
    rows: list[dict[str, object]] = []
    for method in methods:
        if method == "egger" and data.n_instruments < 4:
            rows.append(_ineligible(method, "fewer_than_four_instruments"))
        elif method in {"ivw", "weighted_median"} and data.n_instruments < 3:
            rows.append(_ineligible(method, "fewer_than_three_instruments"))
        else:
            rows.append(_supported(method, data))
    return pd.DataFrame(
        rows,
        columns=[
            "method",
            "estimate",
            "standard_error",
            "p_value",
            "eligible",
            "ineligibility_reason",
        ],
    )
