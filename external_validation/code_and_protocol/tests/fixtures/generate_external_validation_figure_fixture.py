#!/usr/bin/env python3
"""Create a tiny, explicitly synthetic input bundle for figure smoke tests."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


TISSUES = ["QTD000131", "QTD000136", "QTD000251", "QTD000256"]
OUTCOMES = ["CAD", "HF"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    regions = []
    for outcome in OUTCOMES:
        for index in range(1, 3):
            position = 10_000_000 * index + (0 if outcome == "CAD" else 1_000_000)
            regions.append({
                "region_id": f"{outcome}_L{index:02d}", "anchor_outcome": outcome,
                "lead_variant": f"{index}:A:G", "chromosome": str(index),
                "lead_position": position, "region_start": position - 500_000,
                "region_end": position + 500_000, "_synthetic_fixture": "TRUE",
            })
    genes = []
    for index, region in enumerate(regions):
        selected = index < 3
        genes.append({
            "region_id": region["region_id"], "gene_id": f"ENSG_TEST_{index:03d}" if selected else "",
            "gene_selection_status": "SELECTED" if selected else "INELIGIBLE",
            "gene_selection_code": "OK" if selected else "E_NO_ELIGIBLE_GENE",
            "n_shared_complete_candidates": 2 if selected else 0,
            "_synthetic_fixture": "TRUE",
        })

    attempts = []
    reruns = []
    unit_index = 0
    for region, gene in zip(regions, genes):
        for tissue in TISSUES:
            for outcome in OUTCOMES:
                unit_index += 1
                unit_id = f"{region['region_id']}__{tissue}__{outcome}"
                selected = gene["gene_selection_status"] == "SELECTED"
                ready = selected and unit_index % 5 != 0
                attempts.append({
                    "unit_id": unit_id, "region_id": region["region_id"],
                    "anchor_outcome": region["anchor_outcome"], "tissue_id": tissue,
                    "outcome_id": outcome, "pre_gate_status": "PENDING" if selected else "INELIGIBLE",
                    "pre_gate_code": "PENDING" if selected else "E_NO_ELIGIBLE_GENE",
                    "gate_status": "READY" if ready else "INELIGIBLE",
                    "error_codes": "OK" if ready else "E_SYNTHETIC_STOP",
                    "variant_overlap_n": 80 + unit_index * 7 if selected else "",
                    "ld_dimension": 80 + unit_index * 7 if selected else "",
                    "ld_rank": min(502, 75 + unit_index * 6) if selected else "",
                    "posterior_generated": "TRUE" if ready and unit_index % 7 == 0 else "FALSE",
                    "_synthetic_fixture": "TRUE",
                })
                rerun_pass = unit_index != 2
                reruns.append({
                    "unit_id": unit_id, "status": "PASS" if rerun_pass else "FAIL",
                    "byte_identical": str(rerun_pass).upper(),
                    "semantic_identical": str(rerun_pass).upper(),
                    "reason": "" if rerun_pass else "synthetic_gate_record_mismatch",
                    "_synthetic_fixture": "TRUE",
                })

    pd.DataFrame(regions).to_csv(args.out_dir / "regions.tsv", sep="\t", index=False)
    pd.DataFrame(genes).to_csv(args.out_dir / "genes.tsv", sep="\t", index=False)
    pd.DataFrame(attempts).to_csv(args.out_dir / "summary.tsv", sep="\t", index=False)
    pd.DataFrame(reruns).to_csv(args.out_dir / "rerun_agreement.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
