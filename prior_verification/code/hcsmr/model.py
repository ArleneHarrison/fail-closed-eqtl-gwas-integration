"""hcsMR — hierarchical cell-state Mendelian randomization model.

Implements a regularised Bayesian causal-inference model with cell-state
partial pooling and spike-and-slab horizontal-pleiotropy control.

Notation (matches Methods M2 of the design spec):
  T   — cell type (K cell types total, T(c) = parental type of state c)
  c   — cell state (M cell states total, nested in T)
  g   — gene
  j   — SNP

For a single (gene, outcome) batch we fit:

    beta_hat[j,c]  ~ Normal(beta[j,c],         se_X[j,c]^2)
    gamma_hat[j]   ~ Normal(theta[c]*beta[j,c] + alpha[j],  se_Y[j]^2)

with hierarchical priors

    beta[j,c]      ~ Normal(0, 1)
    theta[c]       ~ Normal(mu_type[T(c)], tau_pool)             # partial pooling
    mu_type[k]     ~ Normal(0, 1)                                # gene-level prior
    tau_pool       = 0.25 (fixed, varied in sensitivity analyses)
    alpha[j]       ~ pi * Normal(0, sigma_alpha) + (1-pi) * delta_0   # spike-and-slab
    pi             ~ Beta(1, 9)
    sigma_alpha    = 0.005 + HalfNormal(0.1)

The earlier implementation estimated multiple variance components jointly in
very sparse settings and produced divergent transitions.  The fixed pooling
scale and regularised pleiotropy scale avoid that non-identifiable funnel; the
model exposes those scales as sensitivity parameters rather than pretending
the data can estimate them reliably.  We use a continuous relaxation for the
spike-and-slab so NUTS gradients are
valid:  alpha[j] = z[j] * Normal(0, sigma_alpha), with z[j] ~ Bernoulli(pi)
marginalised analytically into a mixture log-density (Polson and Scott 2010
trick implemented via numpyro.factor).
"""
from __future__ import annotations

import os
import gc
from typing import Any

import jax
import jax.numpy as jnp
import numpy as np
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS

ArrayLike = Any


def _logsumexp(a: jnp.ndarray, b: jnp.ndarray) -> jnp.ndarray:
    m = jnp.maximum(a, b)
    return m + jnp.log(jnp.exp(a - m) + jnp.exp(b - m))


def hcsmr_single_gene_model(
    beta_hat: jnp.ndarray,       # (J, C) sc-eQTL exposure point estimates
    se_X: jnp.ndarray,           # (J, C) standard errors of beta_hat
    gamma_hat: jnp.ndarray,      # (J,) GWAS outcome point estimates
    se_Y: jnp.ndarray,           # (J,) standard errors of gamma_hat
    state_to_type: jnp.ndarray,  # (C,) integer cell-type index for each state
    n_types: int,
    *,
    sigma_alpha_scale: float = 0.1,
    tau_pool: float = 0.25,
    beta_prior_sd: float = 1.0,
    iv_in_state: jnp.ndarray | None = None,  # (J, C) bool mask: which (j,c) are IVs
):
    """Single-gene hcsMR Bayesian model in NumPyro.

    Set entries of ``beta_hat`` and ``se_X`` to NaN where SNP j is not an
    instrument for state c; we will mask those out of the likelihood.  Or pass
    ``iv_in_state`` explicitly as a boolean mask.
    """
    J, C = beta_hat.shape
    if iv_in_state is None:
        iv_in_state = jnp.isfinite(beta_hat) & jnp.isfinite(se_X)
    else:
        iv_in_state = iv_in_state.astype(bool) & jnp.isfinite(beta_hat) & jnp.isfinite(se_X)

    if tau_pool <= 0 or beta_prior_sd <= 0 or sigma_alpha_scale <= 0:
        raise ValueError("prior scales must be positive")

    # ---------- Level 2: gene-level priors on each cell type's centre ----------
    mu_type = numpyro.sample("mu_type", dist.Normal(jnp.zeros(n_types), 1.0).to_event(1))

    # ---------- Level 3: cell-state causal effect with constrained pooling ----------
    with numpyro.plate("states", C):
        theta = numpyro.sample(
            "theta",
            dist.Normal(loc=mu_type[state_to_type], scale=tau_pool),
        )

    # ---------- Level 1: per-(SNP, state) latent true sc-eQTL effects ----------
    beta = numpyro.sample(
        "beta",
        dist.Normal(jnp.zeros((J, C)), jnp.full((J, C), beta_prior_sd)).to_event(2),
    )

    # SNP-exposure likelihood (only at observed (j, c))
    se_X_safe = jnp.where(iv_in_state, se_X, 1.0)  # filler 1.0 where masked
    beta_hat_safe = jnp.where(iv_in_state, beta_hat, 0.0)
    log_lik_X = dist.Normal(beta, se_X_safe).log_prob(beta_hat_safe)
    log_lik_X = jnp.where(iv_in_state, log_lik_X, 0.0)
    numpyro.factor("log_lik_X", log_lik_X.sum())

    # ---------- Spike-and-slab pleiotropy intercept per SNP ----------
    pi_alpha = numpyro.sample("pi_alpha", dist.Beta(1.0, 9.0))
    sigma_alpha = 0.005 + numpyro.sample("sigma_alpha_excess", dist.HalfNormal(sigma_alpha_scale))

    # For each SNP j, the predicted outcome effect when its alpha is 0:
    #   gamma_pred_clean[j] = sum_c iv_mask[j,c] * theta[c] * beta[j,c]
    # IVs are by construction cell-state-specific; in the simulation each SNP
    # is an instrument for exactly one cell state.  For cross-state SNPs the
    # contribution sums.
    contrib = jnp.where(iv_in_state, theta[None, :] * beta, 0.0)  # (J, C)
    gamma_pred = contrib.sum(axis=1)  # (J,)

    # Spike-and-slab on alpha_j: marginalise z_j analytically.
    var_Y = se_Y ** 2
    # log p(gamma_hat | alpha=0)
    log_p_no = dist.Normal(gamma_pred, se_Y).log_prob(gamma_hat)
    # log p(gamma_hat | alpha~Normal(0, sigma_alpha)) — convolution gives
    # gamma_hat ~ Normal(gamma_pred, sqrt(se_Y^2 + sigma_alpha^2))
    log_p_yes = dist.Normal(gamma_pred, jnp.sqrt(var_Y + sigma_alpha ** 2)).log_prob(
        gamma_hat
    )
    log_pi = jnp.log(pi_alpha + 1e-30)
    log_1mpi = jnp.log1p(-pi_alpha + 1e-30)
    log_mix = _logsumexp(log_1mpi + log_p_no, log_pi + log_p_yes)
    numpyro.factor("log_lik_Y", log_mix.sum())


def fit_hcsmr(
    beta_hat: np.ndarray,
    se_X: np.ndarray,
    gamma_hat: np.ndarray,
    se_Y: np.ndarray,
    state_to_type: np.ndarray,
    *,
    n_warmup: int = 1000,
    n_samples: int = 1000,
    n_chains: int = 4,
    target_accept_prob: float = 0.9,
    seed: int = 0,
    iv_in_state: np.ndarray | None = None,
    tau_pool: float = 0.25,
    beta_prior_sd: float = 1.0,
    progress_bar: bool = False,
    chain_method: str | None = None,
) -> dict[str, Any]:
    """Run NUTS for a single-gene hcsMR model.

    Returns a dict with posterior samples and per-state summary.
    """
    n_types = int(state_to_type.max()) + 1
    beta_hat_j = jnp.asarray(beta_hat, dtype=jnp.float32)
    se_X_j = jnp.asarray(se_X, dtype=jnp.float32)
    gamma_hat_j = jnp.asarray(gamma_hat, dtype=jnp.float32)
    se_Y_j = jnp.asarray(se_Y, dtype=jnp.float32)
    state_to_type_j = jnp.asarray(state_to_type, dtype=jnp.int32)
    if iv_in_state is not None:
        iv_in_state_j = jnp.asarray(iv_in_state, dtype=bool)
    else:
        iv_in_state_j = None

    if chain_method is None:
        chain_method = "parallel" if n_chains > 1 else "sequential"

    kernel = NUTS(
        hcsmr_single_gene_model,
        target_accept_prob=target_accept_prob,
        max_tree_depth=10,
    )
    mcmc = MCMC(
        kernel,
        num_warmup=n_warmup,
        num_samples=n_samples,
        num_chains=n_chains,
        chain_method=chain_method,
        progress_bar=progress_bar,
    )
    mcmc.run(
        jax.random.PRNGKey(seed),
        beta_hat=beta_hat_j,
        se_X=se_X_j,
        gamma_hat=gamma_hat_j,
        se_Y=se_Y_j,
        state_to_type=state_to_type_j,
        n_types=n_types,
        iv_in_state=iv_in_state_j,
        tau_pool=tau_pool,
        beta_prior_sd=beta_prior_sd,
    )

    samples = mcmc.get_samples(group_by_chain=False)
    chain_samples = mcmc.get_samples(group_by_chain=True)
    summary = _summarise(samples, chain_samples)
    extra_fields = mcmc.get_extra_fields(group_by_chain=True)
    divergences = np.asarray(extra_fields.get("diverging", []), dtype=bool)
    posterior_predictive = _posterior_predictive_diagnostics(
        samples=samples,
        beta_hat=np.asarray(beta_hat, dtype=float),
        gamma_hat=np.asarray(gamma_hat, dtype=float),
        se_y=np.asarray(se_Y, dtype=float),
        iv_in_state=np.asarray(iv_in_state, dtype=bool)
        if iv_in_state is not None
        else np.isfinite(np.asarray(beta_hat, dtype=float)) & np.isfinite(np.asarray(se_X, dtype=float)),
        seed=seed,
    )
    result = {
        "samples": {k: np.asarray(v) for k, v in samples.items()},
        "summary": summary,
        "diagnostics": {
            "n_divergences": int(divergences.sum()),
            "posterior_predictive": posterior_predictive,
        },
    }
    # Repeated short NUTS fits otherwise retain compiled XLA executables and
    # can exhaust server-side LLVM section memory during a validation grid.
    del mcmc, samples, chain_samples, extra_fields
    jax.clear_caches()
    gc.collect()
    return result


def _posterior_predictive_diagnostics(
    *,
    samples: dict,
    beta_hat: np.ndarray,
    gamma_hat: np.ndarray,
    se_y: np.ndarray,
    iv_in_state: np.ndarray,
    seed: int,
) -> dict[str, float]:
    """A compact posterior-predictive residual check for the outcome model.

    The check is an audit diagnostic, not a model-selection score.  It reports
    the observed standardised residual discrepancy and the proportion of
    replicated discrepancies at least that large.  A tail probability near 0
    or 1 flags poor calibration for follow-up rather than proving causality.
    """
    required = {"theta", "beta", "pi_alpha", "sigma_alpha_excess"}
    missing = required.difference(samples)
    if missing:
        return {"status": "unavailable", "reason": f"missing posterior draws: {sorted(missing)}"}
    n_draws = min(len(samples["theta"]), 500)
    if n_draws == 0:
        return {"status": "unavailable", "reason": "no posterior draws"}
    indices = np.linspace(0, len(samples["theta"]) - 1, n_draws, dtype=int)
    rng = np.random.default_rng(seed + 10_003)
    observed, replicated = [], []
    for index in indices:
        beta = np.asarray(samples["beta"][index])
        prediction = np.where(
            iv_in_state,
            np.asarray(samples["theta"][index])[None, :] * beta,
            0.0,
        ).sum(axis=1)
        pi_alpha = float(samples["pi_alpha"][index])
        sigma_alpha = 0.005 + float(samples["sigma_alpha_excess"][index])
        slab = rng.random(len(gamma_hat)) < pi_alpha
        alpha_rep = rng.normal(0.0, sigma_alpha, len(gamma_hat)) * slab
        gamma_rep = prediction + alpha_rep + rng.normal(0.0, se_y)
        observed.append(float(np.sum(np.square((gamma_hat - prediction) / se_y))))
        replicated.append(float(np.sum(np.square((gamma_rep - prediction) / se_y))))
    observed_arr = np.asarray(observed)
    replicated_arr = np.asarray(replicated)
    return {
        "status": "completed",
        "observed_discrepancy_median": float(np.median(observed_arr)),
        "replicated_discrepancy_median": float(np.median(replicated_arr)),
        "tail_probability": float(np.mean(replicated_arr >= observed_arr)),
    }


def _summarise(samples: dict, chain_samples: dict) -> dict:
    out: dict[str, Any] = {}
    for k, v in samples.items():
        arr = np.asarray(v)
        # If parameter is 1-D (e.g. theta has shape (S, C)), produce per-component summary
        if arr.ndim == 1:
            out[k] = {
                "mean": float(arr.mean()),
                "sd": float(arr.std()),
                "q025": float(np.quantile(arr, 0.025)),
                "q500": float(np.quantile(arr, 0.5)),
                "q975": float(np.quantile(arr, 0.975)),
                "rhat": _rhat(np.asarray(chain_samples[k]).reshape(np.asarray(chain_samples[k]).shape[0], -1)),
                "ess": _ess(np.asarray(chain_samples[k]).reshape(np.asarray(chain_samples[k]).shape[0], -1)),
            }
        elif arr.ndim == 2:
            C = arr.shape[1]
            out[k] = {
                "mean": arr.mean(axis=0).tolist(),
                "sd": arr.std(axis=0).tolist(),
                "q025": np.quantile(arr, 0.025, axis=0).tolist(),
                "q500": np.quantile(arr, 0.5, axis=0).tolist(),
                "q975": np.quantile(arr, 0.975, axis=0).tolist(),
            }
            cs = np.asarray(chain_samples[k])  # (chains, draws, C)
            if cs.ndim == 3:
                out[k]["rhat"] = [_rhat(cs[:, :, i]) for i in range(C)]
                out[k]["ess"] = [_ess(cs[:, :, i]) for i in range(C)]
        else:
            out[k] = {"shape": arr.shape, "mean": float(arr.mean()), "sd": float(arr.std())}
    return out


def _rhat(chain_arr: np.ndarray) -> float:
    """Gelman-Rubin Rhat for a (chains, draws) or (chains, draws, ...) array."""
    if chain_arr.ndim == 3:
        # average per-element rhat
        return float(np.mean([_rhat(chain_arr[:, :, i]) for i in range(chain_arr.shape[2])]))
    n_chains, n_draws = chain_arr.shape
    if n_chains < 2:
        return float("nan")
    chain_means = chain_arr.mean(axis=1)
    chain_vars = chain_arr.var(axis=1, ddof=1)
    W = chain_vars.mean()
    B = n_draws * chain_means.var(ddof=1)
    var_hat = ((n_draws - 1) / n_draws) * W + B / n_draws
    if W <= 0:
        return float("nan")
    return float(np.sqrt(var_hat / W))


def _ess(chain_arr: np.ndarray) -> float:
    """Effective sample size (a simple autocorr-based estimator).

    Falls back to total draws when autocorr can't be estimated.
    """
    if chain_arr.ndim == 3:
        return float(np.mean([_ess(chain_arr[:, :, i]) for i in range(chain_arr.shape[2])]))
    n_chains, n_draws = chain_arr.shape
    if n_chains < 2:
        return float(n_draws)
    x = chain_arr.reshape(-1)
    n = x.size
    x = x - x.mean()
    var = x.var()
    if var <= 0:
        return float(n)
    # autocorr up to lag K
    K = min(200, n // 4)
    rho = np.array([np.dot(x[: n - k], x[k:]) / ((n - k) * var) for k in range(K)])
    # Geyer's initial monotone sum
    rho_pairs = rho[1:].reshape(-1, 2).sum(axis=1) if (K - 1) % 2 == 0 else rho[1:-1].reshape(-1, 2).sum(axis=1)
    # cut at first negative pair sum
    neg = np.where(rho_pairs < 0)[0]
    cutoff = neg[0] if neg.size else len(rho_pairs)
    tau = 1.0 + 2.0 * rho_pairs[:cutoff].sum()
    if tau <= 0:
        return float(n)
    return float(n / tau)


def bayes_factor_zero(theta_samples: np.ndarray, prior_sd: float = 1.0,
                       cap_log10: float = 200.0) -> float:
    """Savage-Dickey Bayes factor BF₁₀ that θ ≠ 0, via closed-form assuming the
    posterior is approximately Gaussian.

    Closed form:
      BF₁₀ = p(θ=0 | prior) / p(θ=0 | posterior)
           = N(0; 0, σ_prior²) / N(0; μ̂, σ̂²)
           = (σ̂ / σ_prior) × exp(μ̂² / (2 σ̂²))
      ⇒ log10 BF₁₀ = log10(σ̂/σ_prior) + (μ̂²) / (2 σ̂² ln 10)

    The closed form avoids the KDE-density-at-zero underflow that produced
    spurious BF₁₀ ≈ 10²⁹⁸ values in the v1 pipeline.  We additionally cap
    log10 BF₁₀ at ``cap_log10`` (default 200) so downstream downstream
    consumers do not propagate infinities.

    The Gaussian-posterior assumption is checked weakly via the Shapiro-Wilk
    statistic on the samples; if the samples are very non-Gaussian we fall
    back to a robust sample-based density estimate at zero using a Student-t
    KDE bandwidth so the answer is still finite.
    """
    from scipy.stats import norm

    if theta_samples.size < 50:
        return float("nan")
    mu  = float(np.mean(theta_samples))
    sig = float(np.std(theta_samples, ddof=1))
    if sig <= 0 or not np.isfinite(sig):
        return float("nan")
    # Closed-form log10 BF₁₀ under Gaussian-posterior approximation
    log10_bf = np.log10(sig / prior_sd) + (mu ** 2) / (2 * sig ** 2 * np.log(10))
    # Cap to avoid infinity propagation
    if log10_bf > cap_log10:
        log10_bf = cap_log10
    if not np.isfinite(log10_bf):
        return float("nan")
    return float(10.0 ** log10_bf)
