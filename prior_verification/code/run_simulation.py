"""Run the hcsMR simulation study across the planned factorial design.

This script is the empirical engine for Figure 3 of the manuscript:
calibration of type-I error and power across factorial scenarios.

Usage on server:
    python run_simulation.py --out SERVER_ACCOUNT_ROOT/hcsmr-cvd/results/sim_v1 --replicates 50

For a full run with replicates=1000 expect a few hours on 112 CPU cores when
using --chain-method parallel.  For a fast sanity check use --replicates 5.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

# Make the hcsmr package importable when running as a script
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

# Force JAX to CPU before importing model.py
os.environ.setdefault("JAX_PLATFORMS", "cpu")
os.environ.setdefault("XLA_FLAGS", "--xla_force_host_platform_device_count=1")
os.environ.setdefault("OMP_NUM_THREADS", "2")

from hcsmr.simulate import SimulationScenario, simulate_scenario  # noqa: E402
from hcsmr.model import fit_hcsmr, bayes_factor_zero  # noqa: E402
from hcsmr.baselines import run_state_baselines, bulk_ivw  # noqa: E402


def make_factorial(replicates: int):
    """Return a list of (scenario_factory, replicate_id) pairs.

    The factorial design follows the spec § Methods M8.2:
      - ivs_per_state ∈ {5, 20, 50}
      - F_stat_mean ∈ {10, 30, 100}
      - pleiotropy_fraction ∈ {0.0, 0.1, 0.3}
      - cell_state_freq ∈ {0.05, 0.01, 0.001}
    × under-null and under-alternative for each.

    For tractability the full 3⁴ × 2 = 162 scenarios are pruned to a 3 × 3
    headline grid in this script (cell_state_freq × pleiotropy_fraction) with
    the other two factors fixed at defaults; the full factorial can be enabled
    via --full-factorial.
    """
    cell_freqs = [0.05, 0.01, 0.001]
    pleios = [0.0, 0.1, 0.3]
    ivs = [20]      # default for headline grid
    fstats = [30.0]
    is_causals = [False, True]

    scenarios = []
    for cf, pl, niv, F, is_causal in itertools.product(cell_freqs, pleios, ivs, fstats, is_causals):
        for rep in range(replicates):
            scenarios.append((cf, pl, niv, F, is_causal, rep))
    return scenarios


def run_one(cf, pl, niv, F, is_causal, rep_seed, *, n_warmup=500, n_samples=500, n_chains=2):
    causal_states = (3,) if is_causal else ()
    s = SimulationScenario(
        name=f"cf={cf}_pl={pl}_iv={niv}_F={F}_caus={is_causal}_rep={rep_seed}",
        n_types=4,
        states_per_type=3,
        causal_states=causal_states,
        causal_theta=0.15 if is_causal else 0.0,
        ivs_per_state=niv,
        F_stat_mean=F,
        pleiotropy_fraction=pl,
        cell_state_freq=cf,
        rng_seed=rep_seed,
    )
    sim = simulate_scenario(s)

    t0 = time.time()
    fit = fit_hcsmr(
        beta_hat=sim["beta_hat"],
        se_X=sim["se_X"],
        gamma_hat=sim["gamma_hat"],
        se_Y=sim["se_Y"],
        state_to_type=sim["state_to_type"],
        iv_in_state=sim["iv_in_state"],
        n_warmup=n_warmup,
        n_samples=n_samples,
        n_chains=n_chains,
        chain_method="sequential",
        progress_bar=False,
        seed=rep_seed,
    )
    elapsed = time.time() - t0

    theta_samples = fit["samples"]["theta"]  # (S, C)
    posterior_means = theta_samples.mean(axis=0)
    posterior_q025 = np.quantile(theta_samples, 0.025, axis=0)
    posterior_q975 = np.quantile(theta_samples, 0.975, axis=0)
    bf10 = np.array([
        bayes_factor_zero(theta_samples[:, c], prior_sd=1.0)
        for c in range(theta_samples.shape[1])
    ])
    rhat_theta = fit["summary"]["theta"].get("rhat") if isinstance(fit["summary"]["theta"], dict) else None

    # baselines
    baselines = run_state_baselines(
        sim["beta_hat"], sim["se_X"], sim["gamma_hat"], sim["se_Y"], sim["iv_states"]
    )
    bulk = bulk_ivw(
        sim["beta_hat"], sim["se_X"], sim["gamma_hat"], sim["se_Y"], sim["iv_states"], sim["state_to_type"]
    )

    # Decision rules: hcsMR positive  iff  BF₁₀ > 10 AND CrI excludes 0
    hcsmr_pos = (bf10 > 10) & ~((posterior_q025 < 0) & (posterior_q975 > 0))
    ivw_pos = np.array([
        (baselines["ivw"][c] or {}).get("p", 1.0) is not None
        and (baselines["ivw"][c] or {}).get("p", 1.0) < 0.05 / theta_samples.shape[1]
        for c in range(theta_samples.shape[1])
    ])

    return {
        "scenario": dict(cell_state_freq=cf, pleiotropy_fraction=pl, ivs_per_state=niv,
                         F_stat_mean=F, is_causal=is_causal, rep=rep_seed),
        "theta_true": sim["theta_true"].tolist(),
        "posterior_mean": posterior_means.tolist(),
        "posterior_q025": posterior_q025.tolist(),
        "posterior_q975": posterior_q975.tolist(),
        "bf10": bf10.tolist(),
        "rhat_theta": rhat_theta if not isinstance(rhat_theta, list) else [float(x) for x in rhat_theta],
        "hcsmr_pos": hcsmr_pos.astype(int).tolist(),
        "ivw_pos": ivw_pos.astype(int).tolist(),
        "bulk_ivw_theta_by_type": {int(k): v.get("theta_hat", float("nan")) for k, v in bulk.items()},
        "elapsed_sec": float(elapsed),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--replicates", type=int, default=20)
    p.add_argument("--n-warmup", type=int, default=500)
    p.add_argument("--n-samples", type=int, default=500)
    p.add_argument("--n-chains", type=int, default=2)
    p.add_argument("--max-scenarios", type=int, default=None,
                   help="cap on scenarios to run (for sanity tests)")
    args = p.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "run.log"
    res_path = out_dir / "results.jsonl"

    scenarios = make_factorial(args.replicates)
    if args.max_scenarios is not None:
        scenarios = scenarios[: args.max_scenarios]

    with res_path.open("w", encoding="utf-8") as fout, log_path.open("w", encoding="utf-8") as flog:
        flog.write(f"running {len(scenarios)} scenarios at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        flog.flush()
        for i, (cf, pl, niv, F, is_causal, rep) in enumerate(scenarios):
            try:
                rec = run_one(cf, pl, niv, F, is_causal, rep,
                              n_warmup=args.n_warmup,
                              n_samples=args.n_samples,
                              n_chains=args.n_chains)
                fout.write(json.dumps(rec) + "\n")
                fout.flush()
                if i % 5 == 0:
                    flog.write(f"[{i+1}/{len(scenarios)}] cf={cf} pl={pl} causal={is_causal} rep={rep} elapsed={rec['elapsed_sec']:.1f}s\n")
                    flog.flush()
            except Exception as e:  # noqa: BLE001
                flog.write(f"[{i+1}/{len(scenarios)}] FAIL cf={cf} pl={pl} causal={is_causal} rep={rep}: {e}\n")
                flog.flush()
        flog.write(f"finished at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")


if __name__ == "__main__":
    main()
