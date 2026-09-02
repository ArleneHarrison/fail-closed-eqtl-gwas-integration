"""Create locked, identically ordered CSV inputs for a coloc-SuSiE run."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from hcsmr.reanalysis.coloc_input_preparation import align_harmonized_to_ld


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--harmonized", required=True)
    parser.add_argument("--ld-npz", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--region-id", required=True)
    parser.add_argument(
        "--molecular-trait-id",
        help="Select one molecular trait when a source has multiple probes per gene.",
    )
    args = parser.parse_args()
    archive = np.load(args.ld_npz, allow_pickle=False)
    variants = [str(value) for value in archive["variant_ids"]]
    aligned, ld, dropped = align_harmonized_to_ld(
        pd.read_csv(args.harmonized),
        variants,
        archive["ld"],
        molecular_trait_id=args.molecular_trait_id,
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    aligned.to_csv(out_dir / "summary.csv", index=False)
    pd.DataFrame(ld, index=aligned["target_variant"], columns=aligned["target_variant"]).to_csv(
        out_dir / "ld.csv"
    )
    (out_dir / "manifest.json").write_text(
        json.dumps(
            {
                "region_id": args.region_id,
                "molecular_trait_id": args.molecular_trait_id,
                "n_variants": len(aligned),
                "dropped_harmonized_variants": dropped,
                "eqtl_sample_size": float(aligned["an"].median() / 2.0),
                "gwas_sample_size": float(aligned["N"].median()),
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
