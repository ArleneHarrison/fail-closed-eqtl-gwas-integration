"""Create schema-stable diagnostic rows from a candidate and variant table."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hcsmr.reanalysis.diagnostics import diagnostic_row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True, help="CSV with candidate_id and f_stat columns")
    parser.add_argument("--out", required=True, help="diagnostic CSV output")
    args = parser.parse_args()
    variants = pd.read_csv(args.variants)
    required = {"candidate_id", "f_stat"}
    missing = required.difference(variants.columns)
    if missing:
        raise ValueError(f"variants missing required columns: {sorted(missing)}")
    rows = [
        diagnostic_row(str(candidate_id), group)
        for candidate_id, group in variants.groupby("candidate_id", sort=False)
    ]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)


if __name__ == "__main__":
    main()
