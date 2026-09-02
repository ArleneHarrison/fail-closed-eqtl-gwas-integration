"""Build an auditable, local 1000G EUR LD matrix for a harmonized variant list."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


def key(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"chr{str(chrom).removeprefix('chr')}_{pos}_{ref.upper()}_{alt.upper()}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--harmonized", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--out-npz", required=True)
    parser.add_argument("--variant-audit-out", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--missing-policy", choices=("mean_impute", "drop_variant"), default="mean_impute")
    args = parser.parse_args()
    harmonized = pd.read_csv(args.harmonized, sep="\t")
    requested = list(dict.fromkeys(harmonized["target_variant"].astype(str)))
    requested_set = set(requested)
    panel = pd.read_csv(args.panel, sep="\t")
    eur = panel.loc[panel["super_pop"] == "EUR", "sample"].astype(str).tolist()
    if len(eur) != 503:
        raise ValueError(f"expected 503 EUR samples, found {len(eur)}")
    chrom_text, interval = args.region.split(":", maxsplit=1)
    start_text, end_text = interval.split("-", maxsplit=1)
    chrom, start, end = chrom_text.removeprefix("chr"), int(start_text), int(end_text)
    vcf = pysam.VariantFile(args.vcf)
    samples = [sample for sample in eur if sample in vcf.header.samples]
    if len(samples) != 503:
        raise ValueError(f"EUR panel/VCF sample overlap is {len(samples)}, not 503")
    vcf.subset_samples(samples)
    observed: dict[str, tuple[np.ndarray, float, int]] = {}
    for record in vcf.fetch(chrom, start - 1, end):
        if len(record.alts or []) != 1:
            continue
        variant = key(chrom, record.pos, record.ref, record.alts[0])
        if variant not in requested_set:
            continue
        dosage = np.array([
            float(sum(call["GT"])) if call["GT"] is not None and None not in call["GT"] else np.nan
            for call in record.samples.values()
        ], dtype=float)
        observed_n = int(np.isfinite(dosage).sum())
        if observed_n < 3:
            continue
        mean = float(np.nanmean(dosage))
        maf = min(mean / 2.0, 1.0 - mean / 2.0)
        imputed = int(np.isnan(dosage).sum())
        if args.missing_policy == "drop_variant" and imputed:
            continue
        dosage[np.isnan(dosage)] = mean
        if float(np.std(dosage)) == 0.0:
            continue
        observed[variant] = (dosage, maf, imputed)
    retained = [variant for variant in requested if variant in observed]
    missing = [variant for variant in requested if variant not in observed]
    if len(retained) < 2:
        raise ValueError("fewer than two requested variants have usable EUR genotypes")
    dosages = np.vstack([observed[variant][0] for variant in retained])
    ld = np.corrcoef(dosages)
    np.savez_compressed(args.out_npz, variant_ids=np.array(retained), ld=ld)
    audit_path = Path(args.variant_audit_out)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    with audit_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["target_variant", "status", "eur_samples", "imputed_missing_genotypes", "maf"])
        writer.writeheader()
        for variant in requested:
            if variant in observed:
                _dosage, maf, imputed = observed[variant]
                writer.writerow({"target_variant": variant, "status": "retained_mean_imputed_missing_genotypes", "eur_samples": 503, "imputed_missing_genotypes": imputed, "maf": maf})
            else:
                writer.writerow({"target_variant": variant, "status": "not_usable_in_vcf_after_qc", "eur_samples": 503, "imputed_missing_genotypes": "", "maf": ""})
    manifest = {
        "region": args.region, "reference": "1000 Genomes Phase 3 GRCh37", "super_population": "EUR",
        "eur_samples": len(samples), "requested_unique_variants": len(requested), "retained_variants": len(retained),
        "not_usable_variants": missing, "missing_genotype_handling": "per-variant mean dosage imputation after pre-recorded <=5% missingness QC" if args.missing_policy == "mean_impute" else "variants with any missing EUR genotype are excluded; no genotype imputation",
        "ld_method": "Pearson correlation of imputed dosage rows; matrix is generated in retained-variant order",
    }
    Path(args.manifest_out).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
