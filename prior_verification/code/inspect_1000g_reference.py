"""Inspect reference alleles at requested eQTL variant coordinates."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pysam


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True)
    parser.add_argument("--vcf-template", required=True)
    parser.add_argument("--index-dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    variants = pd.read_csv(args.variants)
    rows = []
    for variant in variants["variant"].dropna().unique():
        chrom, position, ref, alt = str(variant).split("_", maxsplit=3)
        chrom_number = chrom.removeprefix("chr")
        vcf = pysam.VariantFile(
            args.vcf_template.format(chrom=chrom_number),
            index_filename=str(Path(args.index_dir) / f"chr{chrom_number}.vcf.gz.tbi"),
        )
        matching = list(vcf.fetch(chrom_number, int(position) - 1, int(position)))
        if not matching:
            rows.append({"requested": variant, "reference": "", "status": "missing_position"})
        else:
            for record in matching:
                reference = f"chr{chrom_number}_{record.pos}_{record.ref}_{','.join(record.alts or [])}"
                status = "same" if record.ref == ref and record.alts == (alt,) else "allele_mismatch"
                rows.append({"requested": variant, "reference": reference, "status": status})
    pd.DataFrame(rows).to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
