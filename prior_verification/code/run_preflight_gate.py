"""Command-line entry point for the fail-closed regional-analysis gate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from hcsmr.reanalysis.preflight_gate import run_preflight


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True, help="CSV/TSV summary table with auditable per-variant fields")
    parser.add_argument("--summary-separator", default=",")
    parser.add_argument("--ld-npz", required=True, help="NPZ containing an array named 'ld'")
    parser.add_argument("--provenance-json", required=True, help="JSON array of immutable-input records")
    parser.add_argument("--eqtl-trait-type", required=True, choices=("quant", "cc", "other"))
    parser.add_argument("--gwas-trait-type", required=True, choices=("quant", "cc", "other"))
    parser.add_argument("--case-fraction", type=float, required=True)
    parser.add_argument(
        "--case-fraction-policy",
        choices=("require_constant_per_row", "study_level"),
        default="require_constant_per_row",
        help="Use study_level when a declared study-level fraction is authoritative and per-row N/case fields are diagnostic.",
    )
    parser.add_argument("--max-condition-number", type=float, required=True)
    parser.add_argument(
        "--ld-rank-policy", choices=("require_full_rank", "diagnostic"),
        default="require_full_rank",
        help="A named downstream policy; diagnostic records rank deficiency without treating it as universal ineligibility.",
    )
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    summary = pd.read_csv(args.summary, sep=args.summary_separator, dtype=str)
    for column in summary.columns:
        summary[column] = summary[column].where(summary[column].notna(), None)
    with Path(args.provenance_json).open(encoding="utf-8") as handle:
        provenance = json.load(handle)
    archive = np.load(args.ld_npz, allow_pickle=False)
    result = run_preflight(
        rows=summary.to_dict(orient="records"),
        ld=archive["ld"], ld_variant_ids=[str(value) for value in archive["variant_ids"]] if "variant_ids" in archive.files else None,
        provenance=provenance,
        eqtl_trait_type=args.eqtl_trait_type,
        gwas_trait_type=args.gwas_trait_type,
        case_fraction=args.case_fraction,
        case_fraction_policy=args.case_fraction_policy,
        ld_rank_policy=args.ld_rank_policy,
        max_condition_number=args.max_condition_number,
    )
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "status_codes": result["status_codes"]}, sort_keys=True))


if __name__ == "__main__":
    main()
