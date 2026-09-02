"""Run the locked hcsMR pooling validation grid."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hcsmr.reanalysis.hierarchical_validation import (
    required_hierarchical_validation_scenarios,
    run_validation_replicate,
    summarize_hierarchical_validation,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--replicates", type=int, default=30)
    parser.add_argument("--warmup", type=int, default=200)
    parser.add_argument("--samples", type=int, default=300)
    parser.add_argument("--chains", type=int, default=2)
    parser.add_argument("--target-accept", type=float, default=0.98)
    parser.add_argument("--max-scenarios", type=int, default=None)
    args = parser.parse_args()
    if args.replicates < 1:
        raise ValueError("replicates must be positive")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    scenarios = list(required_hierarchical_validation_scenarios())
    if args.max_scenarios is not None:
        scenarios = scenarios[: args.max_scenarios]
    all_rows = []
    for scenario in scenarios:
        for replicate in range(args.replicates):
            all_rows.append(
                run_validation_replicate(
                    scenario,
                    replicate=replicate,
                    warmup=args.warmup,
                    samples=args.samples,
                    chains=args.chains,
                    target_accept_prob=args.target_accept,
                )
            )
        partial = pd.concat(all_rows, ignore_index=True)
        partial.to_csv(out_dir / "hierarchical_validation.partial.csv", index=False)
        print(f"completed {scenario.name}", flush=True)
    rows = pd.concat(all_rows, ignore_index=True)
    rows.to_csv(out_dir / "hierarchical_validation_replicates.csv.gz", index=False)
    summarize_hierarchical_validation(rows).to_csv(
        out_dir / "hierarchical_validation_summary.csv", index=False
    )


if __name__ == "__main__":
    main()
