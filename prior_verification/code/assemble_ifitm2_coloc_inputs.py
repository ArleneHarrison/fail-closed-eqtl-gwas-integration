"""Assemble ordered coloc inputs, retaining an audit for duplicate source rows."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--harmonized", required=True)
    parser.add_argument("--ld-npz", required=True)
    parser.add_argument("--ld-variant-audit", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    harmonized = pd.read_csv(args.harmonized, sep="\t")
    ld_audit = pd.read_csv(args.ld_variant_audit)
    archive = np.load(args.ld_npz, allow_pickle=False)
    ld_variants = [str(x) for x in archive["variant_ids"]]
    ld = archive["ld"]
    if ld.shape != (len(ld_variants), len(ld_variants)):
        raise ValueError("LD matrix does not match variant order")
    compare = ["eqtl_beta", "eqtl_se", "eqtl_pvalue", "eqtl_an", "gwas_beta_aligned_to_eqtl", "gwas_se", "gwas_pvalue", "gwas_n"]
    duplicate_rows = []
    retained_rows = []
    for variant, group in harmonized.groupby("target_variant", sort=False):
        identical = group[compare].nunique(dropna=False).eq(1).all()
        if len(group) == 1:
            row = group.iloc[0].copy()
            row["duplicate_status"] = "unique"
            retained_rows.append(row)
        elif identical:
            row = group.iloc[0].copy()
            row["duplicate_status"] = "duplicate_identical_source_rows_collapsed"
            retained_rows.append(row)
            duplicate_rows.append({"target_variant": variant, "source_rows": len(group), "status": row["duplicate_status"]})
        else:
            duplicate_rows.append({"target_variant": variant, "source_rows": len(group), "status": "duplicate_nonidentical_excluded"})
    deduplicated = pd.DataFrame(retained_rows)
    maf = ld_audit.loc[ld_audit["status"] == "retained_mean_imputed_missing_genotypes", ["target_variant", "maf"]]
    data = deduplicated.merge(maf, on="target_variant", how="inner", validate="one_to_one")
    by_variant = data.set_index("target_variant")
    ordered = by_variant.loc[[v for v in ld_variants if v in by_variant.index]].reset_index()
    positions = {variant: index for index, variant in enumerate(ld_variants)}
    order_index = [positions[v] for v in ordered["target_variant"]]
    ordered_ld = ld[np.ix_(order_index, order_index)]
    summary = pd.DataFrame({
        "target_variant": ordered["target_variant"], "beta": ordered["eqtl_beta"], "se_eqtl": ordered["eqtl_se"],
        "gwas_beta_aligned": ordered["gwas_beta_aligned_to_eqtl"], "se_gwas": ordered["gwas_se"],
        "N": ordered["gwas_n"], "an": ordered["eqtl_an"], "BP": ordered["position"], "maf": ordered["maf"],
    })
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out / "summary.csv", index=False)
    pd.DataFrame(ordered_ld, index=summary["target_variant"], columns=summary["target_variant"]).to_csv(out / "ld.csv")
    pd.DataFrame(duplicate_rows).to_csv(out / "duplicate_variant_audit.csv", index=False)
    manifest = {"harmonized_rows": len(harmonized), "unique_harmonized_variants": len(deduplicated), "ld_variants": len(ld_variants), "coloc_variants": len(summary), "excluded_nonidentical_duplicates": sum(x["status"] == "duplicate_nonidentical_excluded" for x in duplicate_rows), "case_fraction": "not embedded; derive from retained raw CAD cases/n and record before R run"}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
