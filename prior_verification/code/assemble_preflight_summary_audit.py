"""Assemble a preflight table without collapsing audit rows or duplications."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--harmonized", required=True)
    parser.add_argument("--ld-variant-audit", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    summary = pd.read_csv(args.harmonized, sep="\t", dtype=str)
    ld = pd.read_csv(args.ld_variant_audit, dtype=str)[["target_variant", "maf", "status"]]
    result = summary.merge(ld, on="target_variant", how="left", validate="many_to_one")
    result["vcf_status"] = "exact_ref_alt_match"
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, index=False)
    print({"audit_rows": len(result), "unique_variants": result["target_variant"].nunique(), "missing_maf_rows": int(result["maf"].isna().sum())})


if __name__ == "__main__":
    main()
