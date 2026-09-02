"""Write de-duplicated allele-aware LD requests from a harmonized summary table."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--column", default="target_variant")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    variants = pd.read_csv(args.input)[[args.column]].dropna().drop_duplicates()
    variants.columns = ["variant"]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    variants.to_csv(out, index=False)


if __name__ == "__main__":
    main()
