"""Create diagnostic inputs from the prior exploratory hcsMR screen."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hcsmr.reanalysis.candidate_inputs import build_candidate_variants


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", required=True)
    parser.add_argument("--hf", required=True)
    parser.add_argument("--af", required=True)
    parser.add_argument("--is", dest="stroke", required=True)
    parser.add_argument("--cad", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    hits = pd.read_csv(args.hits, compression="infer")
    tables = {
        "HF": pd.read_csv(args.hf, compression="infer", low_memory=False),
        "AF": pd.read_csv(args.af, compression="infer", low_memory=False),
        "IS": pd.read_csv(args.stroke, compression="infer", low_memory=False),
        "CAD": pd.read_csv(args.cad, compression="infer", low_memory=False),
    }
    variants, status = build_candidate_variants(hits, tables)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    variants.to_csv(out_dir / "candidate_variants.csv", index=False)
    status.to_csv(out_dir / "candidate_input_status.csv", index=False)


if __name__ == "__main__":
    main()
