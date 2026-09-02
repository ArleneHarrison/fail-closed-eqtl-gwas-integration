"""Lift inherited GRCh38 eQTL candidate variants to GRCh37 with full audit rows."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from pyliftover import LiftOver

from hcsmr.reanalysis.liftover_contract import transformed_alleles


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True, help="CSV containing a `variant` column")
    parser.add_argument("--chain", required=True, help="UCSC hg38ToHg19 chain file")
    parser.add_argument("--out", required=True)
    parser.add_argument("--sep", default=",", help="input delimiter, e.g. tab")
    args = parser.parse_args()
    converter = LiftOver(args.chain)
    rows = []
    separator = "\t" if args.sep == "tab" else args.sep
    for variant in pd.read_csv(args.variants, sep=separator)["variant"].dropna().drop_duplicates():
        chrom, position, ref, alt = str(variant).split("_", maxsplit=3)
        mappings = converter.convert_coordinate(chrom, int(position) - 1)
        if len(mappings) != 1:
            rows.append(
                {
                    "source_variant": variant,
                    "status": "unmapped" if not mappings else "ambiguous_mapping",
                    "target_variant": "",
                    "strand": "",
                }
            )
            continue
        target_chrom, target_zero_based, strand, _ = mappings[0]
        target_ref, target_alt = transformed_alleles(ref, alt, strand)
        target_variant = f"chr{target_chrom.removeprefix('chr')}_{target_zero_based + 1}_{target_ref}_{target_alt}"
        rows.append(
            {
                "source_variant": variant,
                "status": "mapped",
                "target_variant": target_variant,
                "strand": strand,
            }
        )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)


if __name__ == "__main__":
    main()
