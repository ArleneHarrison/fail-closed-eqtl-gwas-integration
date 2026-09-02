"""Extract complete regional summary statistics with explicit source schemas."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hcsmr.reanalysis.regional_summary import filter_eqtl_region, filter_gwas_region


def _write_filtered(
    source: str,
    out: str,
    predicate,
    chunksize: int,
    **read_kwargs,
) -> None:
    path = Path(out)
    path.parent.mkdir(parents=True, exist_ok=True)
    first = True
    count = 0
    for chunk in pd.read_csv(source, chunksize=chunksize, **read_kwargs):
        filtered = predicate(chunk)
        if not filtered.empty:
            filtered.to_csv(path, index=False, mode="w" if first else "a", header=first)
            first = False
            count += len(filtered)
    if first:
        pd.DataFrame().to_csv(path, index=False)
    print(f"rows_written={count}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("eqtl", "gwas"), required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--end", required=True, type=int)
    parser.add_argument("--gene-id", default=None)
    parser.add_argument("--gwas-chrom-column", default="CHR")
    parser.add_argument("--gwas-position-column", default="BP")
    parser.add_argument("--compression", default="infer")
    parser.add_argument("--chunksize", type=int, default=200_000)
    args = parser.parse_args()
    if args.kind == "eqtl" and not args.gene_id:
        raise ValueError("--gene-id is required for eQTL extraction")
    read_kwargs = {"sep": "\t", "compression": args.compression, "low_memory": False}
    if args.kind == "eqtl":
        predicate = lambda frame: filter_eqtl_region(
            frame, args.gene_id, args.chrom, args.start, args.end
        )
    else:
        predicate = lambda frame: filter_gwas_region(
            frame,
            args.chrom,
            args.start,
            args.end,
            args.gwas_chrom_column,
            args.gwas_position_column,
        )
    _write_filtered(args.source, args.out, predicate, args.chunksize, **read_kwargs)


if __name__ == "__main__":
    main()
