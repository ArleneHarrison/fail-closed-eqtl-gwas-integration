"""Create an all-or-nothing candidate gate from raw-GWAS QC rows."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hcsmr.reanalysis.raw_gwas_verification import summarize_candidate_qc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    result = summarize_candidate_qc(pd.read_csv(args.qc))
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
