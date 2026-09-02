"""Create compact EUR LD matrices from indexed 1000 Genomes Phase 3 VCFs."""
from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from hcsmr.reanalysis.ld_reference import (
    build_variant_key,
    compute_ld,
    parse_region,
    select_requested_variants,
)


DEFAULT_PANEL_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)
DEFAULT_VCF_TEMPLATE = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)


def eur_samples(panel_url: str) -> list[str]:
    with urllib.request.urlopen(panel_url) as response:
        panel = pd.read_csv(response, sep="\t")
    return panel.loc[panel["super_pop"] == "EUR", "sample"].tolist()


def extract_ld(
    region: str,
    panel_url: str,
    vcf_template: str,
    index_dir: Path | None = None,
    requested_variants: set[str] | None = None,
) -> tuple[list[str], np.ndarray, int, list[str]]:
    chrom, start, end = parse_region(region)
    requested_in_region = None
    if requested_variants is not None:
        requested_in_region = {
            variant
            for variant in requested_variants
            if (parts := str(variant).split("_", maxsplit=3))
            and len(parts) == 4
            and parts[0].lower() == f"chr{chrom}".lower()
            and start <= int(parts[1]) <= end
        }
    vcf_url = vcf_template.format(chrom=chrom)
    index_filename = (
        str(index_dir / f"chr{chrom}.vcf.gz.tbi") if index_dir is not None else None
    )
    vcf = pysam.VariantFile(vcf_url, index_filename=index_filename)
    samples = [sample for sample in eur_samples(panel_url) if sample in vcf.header.samples]
    if len(samples) < 50:
        raise ValueError("fewer than 50 EUR samples available in the reference panel")
    vcf.subset_samples(samples)
    keys, rows = [], []
    for record in vcf.fetch(chrom, start - 1, end):
        if len(record.alts or []) != 1 or len(record.ref) != 1 or len(record.alts[0]) != 1:
            continue
        key = build_variant_key(chrom, record.pos, record.ref, record.alts[0])
        if requested_in_region is not None and key not in requested_in_region:
            continue
        dosage = np.array(
            [
                float(sum(call["GT"])) if call["GT"] is not None and None not in call["GT"] else np.nan
                for call in record.samples.values()
            ]
        )
        if np.mean(np.isfinite(dosage)) < 0.95:
            continue
        maf = np.nanmean(dosage) / 2.0
        maf = min(maf, 1.0 - maf)
        if maf < 0.01:
            continue
        keys.append(key)
        rows.append(dosage)
    if len(rows) < 2:
        raise ValueError(f"insufficient reference variants in {region}")
    selected, missing = select_requested_variants(keys, requested_in_region or set())
    if requested_in_region is not None:
        index = {key: position for position, key in enumerate(keys)}
        rows = [rows[index[key]] for key in selected]
        keys = selected
    if len(rows) < 2:
        raise ValueError(f"fewer than two requested variants matched in {region}; missing={missing}")
    return keys, compute_ld(np.vstack(rows)), len(samples), missing


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--region", action="append", required=True, help="label=chr:start-end")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--panel-url", default=DEFAULT_PANEL_URL)
    parser.add_argument("--vcf-template", default=DEFAULT_VCF_TEMPLATE)
    parser.add_argument("--index-dir", default=None, help="directory containing chrN.vcf.gz.tbi")
    parser.add_argument(
        "--variant-list",
        default=None,
        help="CSV containing a `variant` column; limits LD to requested allele-aware variants",
    )
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    requested_variants = None
    if args.variant_list:
        requested_variants = set(pd.read_csv(args.variant_list)["variant"].dropna())
    for item in args.region:
        label, region = item.split("=", maxsplit=1)
        variants, ld, n_samples, missing = extract_ld(
            region,
            args.panel_url,
            args.vcf_template,
            Path(args.index_dir) if args.index_dir else None,
            requested_variants,
        )
        np.savez_compressed(out_dir / f"{label}.npz", variant_ids=variants, ld=ld)
        (out_dir / f"{label}.json").write_text(
            json.dumps(
                {
                    "label": label,
                    "region": region,
                    "build": "GRCh37",
                    "reference": "1000 Genomes Phase 3",
                    "super_population": "EUR",
                    "n_samples": n_samples,
                    "n_variants": len(variants),
                    "unmatched_requested_variants": missing,
                    "panel_url": args.panel_url,
                    "vcf_template": args.vcf_template,
                },
                indent=2,
            ),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
