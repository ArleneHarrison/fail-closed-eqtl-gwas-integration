"""Join lifted eQTL rows to GWAS records with auditable effect directions."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hcsmr.reanalysis.harmonization import align_outcome_to_exposure


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eqtl", required=True)
    parser.add_argument("--liftover", required=True)
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--audit-out", required=True)
    args = parser.parse_args()
    eqtl = pd.read_csv(args.eqtl, sep="\t")
    lifted = pd.read_csv(args.liftover).query("status == 'mapped'")
    gwas = pd.read_csv(args.gwas)
    merged = (
        eqtl.merge(lifted, left_on="variant", right_on="source_variant", how="inner")
        .merge(gwas, left_on="rsid", right_on="SNP", how="inner", suffixes=("_eqtl", "_gwas"))
    )
    alignments = [
        align_outcome_to_exposure(row.ref, row.alt, row.A1, row.A2, row.b)
        for row in merged.itertuples(index=False)
    ]
    merged["alignment_status"] = [result.status for result in alignments]
    merged["gwas_beta_aligned"] = [result.beta for result in alignments]
    Path(args.audit_out).parent.mkdir(parents=True, exist_ok=True)
    merged[["alignment_status"]].value_counts().rename("n").reset_index().to_csv(args.audit_out, index=False)
    retained = merged.loc[merged["alignment_status"].isin(["aligned", "flipped"])].copy()
    retained.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
