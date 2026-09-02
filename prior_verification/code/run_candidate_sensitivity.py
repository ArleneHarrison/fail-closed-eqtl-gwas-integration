"""Run auditable candidate-level MR sensitivity diagnostics."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hcsmr.reanalysis.candidate_sensitivity import candidate_sensitivity


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--leave-one-out-out", required=True)
    args = parser.parse_args()
    summary, leave_one_out = candidate_sensitivity(pd.read_csv(args.variants))
    for path, table in (
        (Path(args.summary_out), summary),
        (Path(args.leave_one_out_out), leave_one_out),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(path, index=False)


if __name__ == "__main__":
    main()
