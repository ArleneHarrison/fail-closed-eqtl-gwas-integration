"""Run the locked v2 stress-test grid and write auditable CSV outputs."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hcsmr.reanalysis.calibrated_simulations import (
    run_calibrated_replicates,
    summarize_method_performance,
)
from hcsmr.reanalysis.simulation_protocol import required_scenarios


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--max-scenarios", type=int, default=None)
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    all_results, all_summaries = [], []
    scenarios = required_scenarios()
    if args.max_scenarios is not None:
        scenarios = scenarios[: args.max_scenarios]
    for scenario_index, scenario in enumerate(scenarios):
        scenario_results, scenario_summaries = [], []
        for true_theta in (0.0, 0.1):
            results = run_calibrated_replicates(
                scenario,
                replicates=args.replicates,
                true_theta=true_theta,
                seed_offset=1_000_000 * scenario_index + (0 if true_theta == 0 else 500_000),
            )
            summary = summarize_method_performance(results, true_theta=true_theta)
            for field, value in {
                "n_instruments": scenario.n_instruments,
                "pleiotropy": scenario.pleiotropy,
                "ld": scenario.ld,
                "sample_overlap": scenario.sample_overlap,
                "true_theta": true_theta,
            }.items():
                summary[field] = value
            scenario_results.append(results)
            scenario_summaries.append(summary)
        all_results.extend(scenario_results)
        all_summaries.extend(scenario_summaries)
        pd.concat(all_summaries, ignore_index=True).to_csv(
            out_dir / "performance_summary.partial.csv", index=False
        )
        print(
            f"completed scenario {scenario_index + 1}/{len(scenarios)}: "
            f"IVs={scenario.n_instruments}, pleiotropy={scenario.pleiotropy}, "
            f"LD={scenario.ld}, overlap={scenario.sample_overlap}",
            flush=True,
        )
    pd.concat(all_results, ignore_index=True).to_csv(out_dir / "replicate_results.csv.gz", index=False)
    pd.concat(all_summaries, ignore_index=True).to_csv(out_dir / "performance_summary.csv", index=False)


if __name__ == "__main__":
    main()
