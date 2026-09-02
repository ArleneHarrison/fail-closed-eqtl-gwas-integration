"""Run a small, reproducible hcsMR fit and record convergence diagnostics.

This is intentionally a smoke/integration run.  It proves the documented
model can execute and emit posterior predictive diagnostics; it is not a
substitute for the locked, scenario-level validation study.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from hcsmr.model import fit_hcsmr
from hcsmr.simulate import SimulationScenario, simulate_scenario


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--warmup", type=int, default=300)
    parser.add_argument("--samples", type=int, default=500)
    parser.add_argument("--chains", type=int, default=2)
    args = parser.parse_args()
    simulated = simulate_scenario(
        SimulationScenario(
            name="two_state_smoke",
            n_types=1,
            states_per_type=2,
            causal_states=(0,),
            causal_theta=0.15,
            ivs_per_state=5,
            F_stat_mean=30.0,
            pleiotropy_fraction=0.1,
            rng_seed=20260815,
        )
    )
    fitted = fit_hcsmr(
        beta_hat=simulated["beta_hat"],
        se_X=simulated["se_X"],
        gamma_hat=simulated["gamma_hat"],
        se_Y=simulated["se_Y"],
        state_to_type=simulated["state_to_type"],
        iv_in_state=simulated["iv_in_state"],
        n_warmup=args.warmup,
        n_samples=args.samples,
        n_chains=args.chains,
        seed=20260815,
        progress_bar=False,
        chain_method="sequential",
    )
    payload = {
        "scenario": simulated["scenario"],
        "true_theta": simulated["theta_true"].tolist(),
        "theta_summary": fitted["summary"]["theta"],
        "diagnostics": fitted["diagnostics"],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
