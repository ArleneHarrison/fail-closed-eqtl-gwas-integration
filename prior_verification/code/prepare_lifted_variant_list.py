"""Prepare a b37 candidate variant list and per-locus spans from liftover audit rows."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--liftover", required=True)
    parser.add_argument("--variants-out", required=True)
    parser.add_argument("--spans-out", required=True)
    args = parser.parse_args()
    lifted = pd.read_csv(args.liftover)
    mapped = lifted.loc[lifted["status"] == "mapped", ["target_variant"]].rename(
        columns={"target_variant": "variant"}
    )
    fields = mapped["variant"].str.split("_", n=3, expand=True)
    mapped["chrom"] = fields[0].str.removeprefix("chr")
    mapped["position"] = fields[1].astype(int)
    spans = (
        mapped.groupby("chrom", as_index=False)["position"]
        .agg(["min", "max"])
        .reset_index()
        .rename(columns={"min": "start", "max": "end"})
    )
    Path(args.variants_out).parent.mkdir(parents=True, exist_ok=True)
    mapped[["variant"]].drop_duplicates().to_csv(args.variants_out, index=False)
    spans.to_csv(args.spans_out, index=False)


if __name__ == "__main__":
    main()
