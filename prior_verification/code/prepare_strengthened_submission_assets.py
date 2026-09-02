"""Assemble the strengthened manuscript figure/table package from locked outputs."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT / "paper" / "revision_20260830" / "final_package"
REV = ROOT / "paper" / "revision_20260831_strengthened"
PACKAGE = REV / "final_package"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if not PACKAGE.exists():
        shutil.copytree(OLD, PACKAGE)
    figures = PACKAGE / "figures"
    tables = PACKAGE / "tables"
    manuscript_dir = PACKAGE / "manuscript"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    manuscript_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(
        ROOT / "paper" / "revision_20260830" / "manuscript_submission.md",
        manuscript_dir / "manuscript_submission_strengthened.md",
    )

    # Figure 6 is generated directly by the pinned R stress-test script.
    stress = REV / "downstream_stress"
    for suffix in ("png", "svg"):
        shutil.copy2(stress / f"figure_6_downstream_stress_test.{suffix}", figures)
    figure6_manifest = {
        "schema_version": "3.0",
        "figure_number": "6",
        "caption": (
            "Figure 6. Downstream-consumer stress test across 100 deterministic replicates per condition. "
            "Panels show SuSiE-RSS posterior-inclusion-probability (PIP) drift relative to correct input, "
            "PIP assigned to the simulated shared causal variant, and regional PP.H4 from coloc.abf. "
            "Summary-order and declared ancestry/LD substitutions are contract-detectable; internally "
            "corrupted signs with otherwise self-consistent metadata define an intentionally unprotected boundary."
        ),
        "figure_contract": (
            "synthetic downstream-consumer sensitivity only; no empirical association, biological, "
            "therapeutic, clinical, or superiority claim"
        ),
        "sources": [
            {
                "path": str(stress / "downstream_stress_test_replicates.tsv"),
                "sha256": sha256(stress / "downstream_stress_test_replicates.tsv"),
            },
            {
                "path": str(REV / "ancestry_ld" / "manifest.json"),
                "sha256": sha256(REV / "ancestry_ld" / "manifest.json"),
            },
        ],
        "outputs": {
            "png": {
                "path": str(figures / "figure_6_downstream_stress_test.png"),
                "sha256": sha256(figures / "figure_6_downstream_stress_test.png"),
                "dpi": 600,
            },
            "svg": {
                "path": str(figures / "figure_6_downstream_stress_test.svg"),
                "sha256": sha256(figures / "figure_6_downstream_stress_test.svg"),
            },
        },
    }
    (figures / "figure_6_downstream_stress_test_manifest.json").write_text(
        json.dumps(figure6_manifest, indent=2) + "\n", encoding="utf-8"
    )

    table2_path = tables / "table_2_validation_summary.tsv"
    table2 = pd.read_csv(table2_path, sep="\t")
    table2.loc[table2["Validation component"] == "Local test suite", ["Cases", "Observed result"]] = [
        "93 tests", "93/93 passed"
    ]
    downstream_mask = table2["Validation component"] == "Downstream model validation"
    table2.loc[downstream_mask, :] = [
        "Synthetic downstream consumer stress test",
        "5 conditions x 100 replicates",
        "500/500 SuSiE-RSS and coloc.abf evaluations completed",
        "Selected input corruptions produce quantifiable consumer-output drift",
    ]
    table2 = table2.loc[table2["Validation component"] != "External biological validation"].copy()
    table2.loc[len(table2)] = [
        "External biological validation", "Not performed", "No result",
        "No association, causal, biological, therapeutic, or clinical claim",
    ]
    table2.to_csv(table2_path, sep="\t", index=False)

    # Update the evidence table because a synthetic downstream stress test now exists.
    table4 = pd.DataFrame(
        [
            ["Enumerated software checks", "All fixed expectations matched", "Named implementation paths execute as specified", "Absence of untested defects"],
            ["Independent numerical path", "R/SVD and gate agreed within tolerance", "Rank and condition-number calculations are numerically concordant", "External end-to-end workflow validation"],
            ["Permutation testing", "All diagnostics invariant", "Numerical output is invariant to simultaneous relabeling", "Correctness under mismatched summary/LD ordering"],
            ["Stored artifacts", "C2 and 11 windows conformed post hoc", "Traceable behavior under one named policy", "Prevalence, generalizability, or model eligibility"],
            ["Cross-ancestry LD transfer", "Five 1000G panels rebuilt for fixed variants", "LD-panel differences are measurable", "Multi-ancestry eQTL/GWAS replication"],
            ["Synthetic downstream consumers", "500/500 SuSiE-RSS and coloc.abf runs completed", "Selected corruptions can change model output", "Biological truth, empirical calibration, or clinical validity"],
            ["Content-level sign corruption", "Large output drift despite self-consistent metadata", "A documented fail-closed boundary", "Guaranteed detection of internally false statistics"],
        ],
        columns=["Evidence component", "Observed result", "Supported claim", "Not supported"],
    )
    table4.to_csv(tables / "table_4_claim_boundaries.tsv", sep="\t", index=False)

    pop = pd.read_csv(REV / "ancestry_ld" / "population_summary.tsv", sep="\t")
    transfer = pd.read_csv(REV / "ancestry_ld" / "ld_transfer_summary.tsv", sep="\t")
    ancestry = pop.merge(
        transfer[
            [
                "population", "common_variants", "normalized_frobenius_distance",
                "median_absolute_delta_r", "p95_absolute_delta_r",
                "sign_discordance_informative_pairs",
            ]
        ],
        on="population",
        how="left",
    )
    ancestry = ancestry[
        [
            "population", "samples", "retained_variants", "common_variants",
            "normalized_frobenius_distance", "median_absolute_delta_r",
            "p95_absolute_delta_r", "sign_discordance_informative_pairs",
        ]
    ]
    ancestry.columns = [
        "Panel", "Samples", "Retained of 234", "Common benchmark variants",
        "Normalized Frobenius distance vs EUR", "Median absolute delta r",
        "95th percentile absolute delta r", "Sign discordance in informative pairs",
    ]
    for column in ancestry.columns[4:]:
        ancestry[column] = ancestry[column].map(lambda value: f"{value:.3f}")
    ancestry.to_csv(tables / "table_6_cross_ancestry_ld.tsv", sep="\t", index=False)

    downstream = pd.read_csv(stress / "downstream_stress_test_summary.tsv", sep="\t")
    downstream = downstream[
        [
            "condition", "completed", "contract_detectable", "median_causal_pip",
            "causal_in_95_cs_rate", "lead_changed_rate_vs_correct",
            "median_pip_l1_drift_vs_correct", "median_max_pp_h4",
        ]
    ]
    downstream.columns = [
        "Input condition", "Completed of 100", "Detectable from contract metadata",
        "Median causal PIP", "Causal variant in 95% CS", "Lead changed vs correct",
        "Median PIP L1 drift", "Median regional PP.H4",
    ]
    downstream["Detectable from contract metadata"] = downstream["Detectable from contract metadata"].map(
        {True: "Yes", False: "No"}
    )
    downstream["Input condition"] = downstream["Input condition"].map(
        {
            "correct": "Correct contract",
            "summary_order_permuted": "Summary order permuted",
            "ten_percent_sign_flipped": "10% signs reversed",
            "EAS_LD_substituted": "EAS LD substituted",
            "AFR_LD_substituted": "AFR LD substituted",
        }
    )
    for column in downstream.columns[3:]:
        downstream[column] = downstream[column].map(lambda value: f"{value:.3f}")
    downstream.to_csv(tables / "table_7_downstream_stress.tsv", sep="\t", index=False)

    resource = pd.read_csv(REV / "cross_resource" / "cross_resource_status_summary.tsv", sep="\t")
    resource.columns = ["GWAS outcome", "Molecular-QTL resource family", "Workflow status", "Cases"]
    resource["Workflow status"] = resource["Workflow status"].replace(
        {"completed_no_comparable_credible_sets": "completed; no comparable CS pairs"}
    )
    resource.to_csv(tables / "table_8_cross_resource_status.tsv", sep="\t", index=False)

    manifest_path = tables / "table_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["schema_version"] = "4.0"
    manifest["tables"][3]["caption"] = (
        "The synthetic downstream stress test quantifies model sensitivity but does not establish "
        "empirical association, biological validity, clinical utility, or comparative superiority."
    )
    manifest["tables"] = [
        item for item in manifest["tables"] if str(item.get("number")) not in {"6", "7", "8"}
    ]
    manifest["tables"].extend(
        [
            {
                "number": "6", "title": "Cross-ancestry LD-panel portability for fixed regional variants",
                "path": "table_6_cross_ancestry_ld.tsv",
                "column_widths_inches": [0.45, 0.55, 0.65, 0.75, 1.1, 0.75, 0.9, 1.15],
                "caption": (
                    "Matrices were reconstructed from 1000 Genomes Phase 3 GRCh37 for the same fixed variants. "
                    "This is an LD-reference portability test, not multi-ancestry molecular-QTL or GWAS validation."
                ),
            },
            {
                "number": "7", "title": "Synthetic downstream-consumer stress test",
                "path": "table_7_downstream_stress.tsv",
                "column_widths_inches": [1.15, 0.65, 0.95, 0.7, 0.75, 0.8, 0.75, 0.8],
                "caption": (
                    "One shared causal variant was simulated under EUR LD in each replicate. SuSiE-RSS PIP metrics "
                    "and coloc.abf regional PP.H4 describe sensitivity to controlled input corruption only."
                ),
            },
            {
                "number": "8", "title": "Retrospective cross-resource workflow-status audit",
                "path": "table_8_cross_resource_status.tsv",
                "column_widths_inches": [0.8, 1.5, 2.25, 0.55],
                "caption": (
                    "Eight archived cases span two GWAS outcomes, three molecular-QTL resource families, and seven "
                    "contexts. Posterior values are deliberately omitted; statuses are not biological results."
                ),
            },
        ]
    )
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    checksum_path = PACKAGE / "SHA256SUMS_strengthened.txt"
    files = [path for path in PACKAGE.rglob("*") if path.is_file() and path != checksum_path]
    checksum_rows = [f"{sha256(path)}  {path.relative_to(PACKAGE).as_posix()}" for path in sorted(files)]
    checksum_path.write_text("\n".join(checksum_rows) + "\n", encoding="utf-8")
    print(json.dumps({"files": len(files), "package": str(PACKAGE)}, sort_keys=True))


if __name__ == "__main__":
    main()
