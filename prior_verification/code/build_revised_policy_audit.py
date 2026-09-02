"""Re-audit stored regional artifacts under the corrected, explicitly scoped gate."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd

from hcsmr.reanalysis.preflight_gate import prepare_analysis_rows, run_preflight


def read_rows(path: Path) -> list[dict]:
    frame = pd.read_csv(path, sep="\t" if path.suffix == ".tsv" else ",", dtype=str)
    frame = frame.where(frame.notna(), None)
    return frame.to_dict(orient="records")


def audit_case(label: str, summary_path: Path, ld_path: Path, provenance: list[dict], out_dir: Path) -> dict:
    rows = read_rows(summary_path)
    archive = np.load(ld_path, allow_pickle=False)
    ld_ids = [str(value) for value in archive["variant_ids"]]
    analysis_rows, row_audit = prepare_analysis_rows(rows, ld_variant_ids=ld_ids)
    result = run_preflight(
        rows=analysis_rows,
        ld=archive["ld"],
        ld_variant_ids=ld_ids,
        provenance=provenance,
        eqtl_trait_type="quant",
        gwas_trait_type="cc",
        case_fraction=181_522 / 1_165_690,
        case_fraction_policy="study_level",
        max_condition_number=1e12,
        ld_rank_policy="diagnostic",
    )
    case_dir = out_dir / label
    case_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(row_audit).to_csv(case_dir / "analysis_view_row_audit.tsv", sep="\t", index=False)
    (case_dir / "revised_preflight_status.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    decisions = pd.Series([item["decision"] for item in row_audit]).value_counts().to_dict()
    diagnostics = result["ld_diagnostics"]
    return {
        "case_id": label,
        "raw_summary_rows": len(rows),
        "analysis_rows": len(analysis_rows),
        "ld_variants": len(ld_ids),
        "exact_duplicates_excluded": decisions.get("excluded_exact_duplicate", 0),
        "rows_not_in_locked_ld_excluded": decisions.get("excluded_not_in_locked_ld", 0),
        "discordant_duplicates_retained": decisions.get("retained_discordant_duplicate_for_gate_failure", 0),
        "numerical_rank": diagnostics["numerical_rank"],
        "rank_deficient": diagnostics["rank_deficient"],
        "condition_number_2": diagnostics["condition_number_2"],
        "rank_policy": diagnostics["rank_policy"],
        "case_fraction_policy": result["case_fraction_policy"],
        "gate_status": result["status"],
        "status_codes": ";".join(result["status_codes"]),
        "model_eligibility": result["model_eligibility"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    coordinate = root / "paper/reanalysis_outputs/new_manuscript_20260826/coordinate_control"
    provenance = json.loads((coordinate / "input_provenance.json").read_text(encoding="utf-8"))
    summaries = [audit_case(
        "coordinate_case",
        coordinate / "preflight_summary.csv",
        coordinate / "ld_mean_impute.npz",
        provenance,
        out_dir,
    )]
    regions = root / "paper/reanalysis_outputs/new_manuscript_20260828/multiregion_transfer/executed/regions"
    for index in range(1, 13):
        label = f"R{index:02d}"
        region = regions / label
        if not (region / "ld_mean_impute.npz").is_file():
            summaries.append({
                "case_id": label, "raw_summary_rows": 0, "analysis_rows": 0, "ld_variants": 0,
                "exact_duplicates_excluded": 0, "rows_not_in_locked_ld_excluded": 0,
                "discordant_duplicates_retained": 0, "numerical_rank": None,
                "rank_deficient": None, "condition_number_2": None, "rank_policy": "diagnostic",
                "case_fraction_policy": "study_level", "gate_status": "INELIGIBLE",
                "status_codes": "E_NO_ELIGIBLE_GENE", "model_eligibility": None,
            })
            continue
        summaries.append(audit_case(
            label,
            region / "preflight_summary.tsv",
            region / "ld_mean_impute.npz",
            provenance,
            out_dir,
        ))
    output = out_dir / "revised_policy_audit_summary.tsv"
    pd.DataFrame(summaries).to_csv(output, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    print(output)


if __name__ == "__main__":
    main()
