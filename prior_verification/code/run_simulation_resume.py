"""Memory-frugal resume of the simulation study.

Differences from run_simulation.py:
  - Smaller per-step memory footprint (n_chains=1, smaller warmup)
  - JAX preallocation disabled
  - Single-threaded XLA
  - Skips scenarios already in results.jsonl (resume on append)
  - Tunable cell-state-frequency subset via --skip-cf

Use:
    python run_simulation_resume.py \
        --out SERVER_ACCOUNT_ROOT/hcsmr-cvd/results/sim_v1 \
        --replicates 8 \
        --skip-cf 0.05      # already complete from prior run
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from pathlib import Path

# JAX memory hygiene
os.environ.setdefault("JAX_PLATFORMS", "cpu")
os.environ.setdefault("XLA_FLAGS", "--xla_force_host_platform_device_count=1")
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from hcsmr.simulate import SimulationScenario, simulate_scenario  # noqa: E402
from hcsmr.model import fit_hcsmr, bayes_factor_zero  # noqa: E402
from hcsmr.baselines import run_state_baselines, bulk_ivw  # noqa: E402


def already_done_keys(jsonl_path: Path) -> set:
    done = set()
    if not jsonl_path.exists():
        return done
    with jsonl_path.open(encoding="utf-8") as f:
        for line in f:
            try:
                r = json.loads(line)
                sc = r["scenario"]
                key = (sc["cell_state_freq"], sc["pleiotropy_fraction"], sc["ivs_per_state"],
                       sc["F_stat_mean"], sc["is_causal"], sc["rep"])
                done.add(key)
            except Exception:
                pass
    return done


def make_factorial(replicates: int, skip_cf=()):
    cell_freqs = [cf for cf in [0.05, 0.01, 0.001] if cf not in skip_cf]
    pleios = [0.0, 0.1, 0.3]
    ivs = [20]
    fstats = [30.0]
    is_causals = [False, True]
    scenarios = []
    for cf, pl, niv, F, is_causal in itertools.product(cell_freqs, pleios, ivs, fstats, is_causals):
        for rep in range(replicates):
            scenarios.append((cf, pl, niv, F, is_causal, rep))
    return scenarios


def run_one(cf, pl, niv, F, is_causal, rep_seed, *, n_warmup, n_samples, n_chains):
    causal_states = (3,) if is_causal else ()
    s = SimulationScenario(
        name=f"cf={cf}_pl={pl}_iv={niv}_F={F}_caus={is_causal}_rep={rep_seed}",
        n_types=4, states_per_type=3,
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
        beta_hat=sim["beta_hat"], se_X=sim["se_X"],
        gamma_hat=sim["gamma_hat"], se_Y=sim["se_Y"],
        state_to_type=sim["state_to_type"], iv_in_state=sim["iv_in_state"],
        n_warmup=n_warmup, n_samples=n_samples, n_chains=n_chains,
        chain_method="sequential", progress_bar=False, seed=rep_seed,
    )
    elapsed = time.time() - t0
    theta_samples = fit["samples"]["theta"]
    posterior_means = theta_samples.mean(axis=0)
    posterior_q025 = np.quantile(theta_samples, 0.025, axis=0)
    posterior_q975 = np.quantile(theta_samples, 0.975, axis=0)
    bf10 = np.array([
        bayes_factor_zero(theta_samples[:, c], prior_sd=1.0)
        for c in range(theta_samples.shape[1])
    ])
    baselines = run_state_baselines(
        sim["beta_hat"], sim["se_X"], sim["gamma_hat"], sim["se_Y"], sim["iv_states"]
    )
    bulk = bulk_ivw(
        sim["beta_hat"], sim["se_X"], sim["gamma_hat"], sim["se_Y"], sim["iv_states"], sim["state_to_type"]
    )
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
        "hcsmr_pos": hcsmr_pos.astype(int).tolist(),
        "ivw_pos": ivw_pos.astype(int).tolist(),
        "bulk_ivw_theta_by_type": {int(k): v.get("theta_hat", float("nan")) for k, v in bulk.items()},
        "elapsed_sec": float(elapsed),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--replicates", type=int, default=8)
    p.add_argument("--n-warmup", type=int, default=300)
    p.add_argument("--n-samples", type=int, default=300)
    p.add_argument("--n-chains", type=int, default=1)
    p.add_argument("--skip-cf", type=float, nargs="*", default=[])
    args = p.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    res_path = out_dir / "results.jsonl"
    log_path = out_dir / "run_resume.log"

    done = already_done_keys(res_path)
    scenarios = make_factorial(args.replicates, skip_cf=tuple(args.skip_cf))

    fresh = [s for s in scenarios if (s[0], s[1], s[2], s[3], s[4], s[5]) not in done]

    with res_path.open("a", encoding="utf-8") as fout, log_path.open("w", encoding="utf-8") as flog:
        flog.write(f"resume: {len(scenarios)} total, {len(scenarios) - len(fresh)} already done, "
                   f"{len(fresh)} to run; at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        flog.flush()
        for i, (cf, pl, niv, F, is_causal, rep) in enumerate(fresh):
            try:
                rec = run_one(cf, pl, niv, F, is_causal, rep,
                              n_warmup=args.n_warmup,
                              n_samples=args.n_samples,
                              n_chains=args.n_chains)
                fout.write(json.dumps(rec) + "\n")
                fout.flush()
                if i % 3 == 0:
                    flog.write(f"[{i+1}/{len(fresh)}] cf={cf} pl={pl} causal={is_causal} rep={rep} elapsed={rec['elapsed_sec']:.1f}s\n")
                    flog.flush()
            except Exception as e:
                flog.write(f"[{i+1}/{len(fresh)}] FAIL cf={cf} pl={pl} causal={is_causal} rep={rep}: {e}\n")
                flog.flush()
        flog.write(f"done at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")


if __name__ == "__main__":
    main()
