"""Harmonize lifted eQTL rows with coordinate-based GWAS summary statistics."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hcsmr.reanalysis.harmonization import (
    build_coordinate_coloc_audit,
    build_coordinate_coloc_table,
)
from hcsmr.reanalysis.gwas_prefilter import resolve_delimiter


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eqtl", required=True)
    parser.add_argument("--liftover", required=True)
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--gwas-sep", default=",", help="CSV delimiter or `tab`")
    parser.add_argument("--out", required=True)
    parser.add_argument("--audit-out", required=True)
    parser.add_argument("--audit-rows-out", required=True,
                        help="write every source row, including failed or ambiguous rows")
    parser.add_argument("--gwas-chrom-column", required=True)
    parser.add_argument("--gwas-position-column", required=True)
    parser.add_argument("--gwas-effect-column", required=True)
    parser.add_argument("--gwas-other-column", required=True)
    parser.add_argument("--gwas-beta-column", required=True)
    parser.add_argument("--gwas-se-column", required=True)
    parser.add_argument("--gwas-n-column", required=True)
    args = parser.parse_args()
    eqtl = pd.read_csv(args.eqtl, sep="\t")
    lifted = pd.read_csv(args.liftover)
    gwas = pd.read_csv(args.gwas, sep=resolve_delimiter(args.gwas_sep))
    audit = build_coordinate_coloc_audit(
        eqtl,
        lifted,
        gwas,
        gwas_chrom_column=args.gwas_chrom_column,
        gwas_position_column=args.gwas_position_column,
        gwas_effect_column=args.gwas_effect_column,
        gwas_other_column=args.gwas_other_column,
        gwas_beta_column=args.gwas_beta_column,
        gwas_se_column=args.gwas_se_column,
        gwas_n_column=args.gwas_n_column,
    )
    merged = audit.loc[audit["alignment_status"].isin(["aligned", "flipped"])].copy()
    Path(args.audit_out).parent.mkdir(parents=True, exist_ok=True)
    audit[["alignment_status"]].value_counts().rename("n").reset_index().to_csv(args.audit_out, index=False)
    audit.to_csv(args.audit_rows_out, index=False)
    merged.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
