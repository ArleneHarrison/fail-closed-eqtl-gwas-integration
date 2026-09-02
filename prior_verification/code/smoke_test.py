"""Quick smoke test: verify NumPyro hcsMR fit runs end-to-end on a small scenario."""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

os.environ.setdefault("JAX_PLATFORMS", "cpu")
os.environ.setdefault("XLA_FLAGS", "--xla_force_host_platform_device_count=1")
os.environ.setdefault("OMP_NUM_THREADS", "4")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np

from hcsmr.simulate import SimulationScenario, simulate_scenario
from hcsmr.model import fit_hcsmr, bayes_factor_zero
from hcsmr.baselines import run_state_baselines


def main():
    print("=== hcsMR smoke test ===")
    # Tiny scenario: 3 cell types × 2 states each = 6 states, 10 IVs each.
    # State index 2 has a non-zero causal theta.
    s = SimulationScenario(
        name="smoke",
        n_types=3,
        states_per_type=2,
        causal_states=(2,),
        causal_theta=0.2,
        ivs_per_state=10,
        F_stat_mean=30.0,
        pleiotropy_fraction=0.1,
        cell_state_freq=0.05,
        rng_seed=42,
    )
    sim = simulate_scenario(s)
    print(f"simulated: J={sim['iv_states'].size} IVs, C={sim['theta_true'].size} states")
    print("theta_true:", sim["theta_true"])

    t0 = time.time()
    fit = fit_hcsmr(
        beta_hat=sim["beta_hat"],
        se_X=sim["se_X"],
        gamma_hat=sim["gamma_hat"],
        se_Y=sim["se_Y"],
        state_to_type=sim["state_to_type"],
        iv_in_state=sim["iv_in_state"],
        n_warmup=300,
        n_samples=300,
        n_chains=2,
        chain_method="sequential",
        progress_bar=True,
        seed=42,
    )
    dt = time.time() - t0
    print(f"fit done in {dt:.1f}s")

    theta = fit["samples"]["theta"]
    print(f"theta posterior mean: {theta.mean(axis=0)}")
    print(f"theta posterior 2.5/97.5: \n  {np.quantile(theta, 0.025, axis=0)}\n  {np.quantile(theta, 0.975, axis=0)}")

    bf = [bayes_factor_zero(theta[:, c]) for c in range(theta.shape[1])]
    print("BF₁₀ per state:", [f"{b:.2f}" for b in bf])

    # baselines
    bl = run_state_baselines(sim["beta_hat"], sim["se_X"], sim["gamma_hat"], sim["se_Y"], sim["iv_states"])
    print("\nIVW per state:")
    for c, r in enumerate(bl["ivw"]):
        if r:
            print(f"  state {c}: theta_hat={r['theta_hat']:.3f}  p={r['p']:.3g}  n_iv={r['n_iv']}")

    print("\nSMOKE TEST OK" if not np.any(np.isnan(theta.mean(axis=0))) else "\nSMOKE TEST FAIL")


if __name__ == "__main__":
    main()
