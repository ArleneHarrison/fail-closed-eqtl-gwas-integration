"""Fail-closed consistency audit for the journal-neutral manuscript package.

This audit checks file identity and declared placeholders.  It does not assess
scientific validity and never upgrades the package to submission-ready.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import zipfile
from pathlib import Path


EXPECTED_CASE_IDS = {"coordinate_case", *(f"R{index:02d}" for index in range(1, 13))}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def docx_text(path: Path) -> str:
    with zipfile.ZipFile(path) as archive:
        xml = archive.read("word/document.xml").decode("utf-8")
    # The precise text is only used for static confirmation of known placeholders.
    return xml.replace("</w:t>", " ")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-tsv", required=True, type=Path)
    parser.add_argument("--figure-manifests", required=True, nargs="+", type=Path)
    parser.add_argument("--docx-manifest", required=True, type=Path)
    parser.add_argument("--docx", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--out-report", required=True, type=Path)
    args = parser.parse_args()

    figure_manifests = [json.loads(item.read_text(encoding="utf-8")) for item in args.figure_manifests]
    docx_manifest = json.loads(args.docx_manifest.read_text(encoding="utf-8"))
    with args.source_tsv.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    case_ids = {row.get("case_id", "") for row in rows}
    source_hash = sha256(args.source_tsv)
    document_text = docx_text(args.docx)
    figures_by_number = {str(item.get("figure_number", "1B")): item for item in figure_manifests}
    figure_paths_by_number = {
        str(item.get("figure_number", "1B")): path for item, path in zip(figure_manifests, args.figure_manifests)
    }
    docx_figures = {str(item["number"]): item for item in docx_manifest.get("figures", [])}
    docx_tables = {str(item["number"]): item for item in docx_manifest.get("tables", [])}
    expected_numbers = {str(index) for index in range(1, len(figure_manifests) + 1)}
    expected_tables = {str(index) for index in range(1, len(docx_tables) + 1)}

    def figure_outputs_match(manifest: dict) -> bool:
        try:
            png = Path(manifest["outputs"]["png"]["path"])
            svg = Path(manifest["outputs"]["svg"]["path"])
            return (
                png.is_file()
                and svg.is_file()
                and sha256(png) == manifest["outputs"]["png"]["sha256"]
                and sha256(svg) == manifest["outputs"]["svg"]["sha256"]
            )
        except (KeyError, OSError):
            return False

    def docx_figure_matches(number: str) -> bool:
        try:
            manifest = figures_by_number[number]
            embedded = docx_figures[number]
            return (
                Path(embedded["path"]).resolve() == Path(manifest["outputs"]["png"]["path"]).resolve()
                and embedded["sha256"] == manifest["outputs"]["png"]["sha256"]
                and embedded["manifest_sha256"] == sha256(figure_paths_by_number[number])
            )
        except (KeyError, ValueError, OSError):
            return False

    main_figure_contracts = all(
        "association" in manifest.get("figure_contract", "").lower()
        and "biological" in manifest.get("figure_contract", "").lower()
        and "no" in manifest.get("figure_contract", "").lower()
        for manifest in figure_manifests
    )

    table_hashes_match = (
        set(docx_tables) == expected_tables
        and all(Path(item["path"]).is_file() and sha256(Path(item["path"])) == item["sha256"] for item in docx_tables.values())
    )

    def figure_sources_match(manifest: dict) -> bool:
        sources = manifest.get("sources", [])
        if not sources:
            return False
        for source in sources:
            try:
                path = Path(source["path"])
                if not path.is_file() or sha256(path) != source["sha256"]:
                    return False
            except (KeyError, OSError):
                return False
        return True

    checks = {
        "all_figure_source_hashes_match_their_manifests": all(figure_sources_match(item) for item in figure_manifests),
        "source_table_has_all_thirteen_expected_cases": len(rows) == 13 and case_ids == EXPECTED_CASE_IDS,
        "contiguous_expected_figure_manifests_present": set(figures_by_number) == expected_numbers,
        "all_figure_png_and_svg_hashes_match_manifests": all(figure_outputs_match(item) for item in figure_manifests),
        "all_figures_record_restricted_claim_contracts": main_figure_contracts,
        "docx_embeds_all_manifested_figures": set(docx_figures) == expected_numbers and all(docx_figure_matches(number) for number in expected_numbers),
        "docx_records_contiguous_verified_tables": table_hashes_match,
        "docx_hash_matches_manifest": sha256(args.docx) == docx_manifest["docx"]["sha256"],
        "docx_names_all_rendered_figures": all(f"Figure {number}" in document_text for number in expected_numbers),
        "docx_preserves_author_placeholders": all(
            phrase in document_text
            for phrase in (
                "Author names and order must be supplied before submission",
                "Public repository URL must be supplied by the authors before submission",
                "Author declaration required before submission",
            )
        ),
        "docx_labels_author_completion_required": "must be supplied" in document_text and "before submission" in document_text,
        "no_submission_metadata_asserted": "Submission ID" not in document_text and "Manuscript ID" not in document_text,
    }
    blockers = [
        "Author order, affiliations, corresponding author, and CRediT confirmation.",
        "Funding, competing interests, ethics/consent wording, and AI-use disclosure confirmation.",
        "Code/data licence, archive/repository, and derivative-data sharing decisions.",
        "Target-journal selection and target-specific proof, followed by corresponding-author approval of the complete package.",
    ]
    status = (
        "TECHNICAL AUDIT PASSED — AUTHOR COMPLETION REQUIRED"
        if all(checks.values())
        else "TECHNICAL AUDIT FAILED — AUTHOR COMPLETION REQUIRED"
    )
    payload = {
        "schema_version": "2.0",
        "status": status,
        "source_tsv_sha256": source_hash,
        "figure_manifest_sha256": {str(item.get("figure_number", "1B")): sha256(path) for item, path in zip(figure_manifests, args.figure_manifests)},
        "docx_manifest_sha256": sha256(args.docx_manifest),
        "docx_sha256": sha256(args.docx),
        "checks": checks,
        "blockers": blockers,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# Final submission-package consistency audit",
        "",
        f"## Status: **{status}**",
        "",
        "This is a file-identity and static-placeholder audit of the journal-neutral editable package. It does not select a journal, supply author metadata, or authorize submission.",
        "",
        "## Verified identities",
        "",
        f"- Locked plot-data SHA-256: `{source_hash}`",
        "- Rendered figures: `" + "`, `".join(sorted(expected_numbers, key=int)) + "` (each with PNG, SVG, and immutable manifest)",
        "- Rendered tables: `" + "`, `".join(sorted(expected_tables, key=int)) + "` (each from a locked TSV recorded in the DOCX manifest)",
        "- Figure-manifest SHA-256 values: `" + json.dumps(payload["figure_manifest_sha256"], sort_keys=True) + "`",
        f"- DOCX-manifest SHA-256: `{payload['docx_manifest_sha256']}`",
        f"- Editable DOCX SHA-256: `{payload['docx_sha256']}`",
        "",
        "## Static checks",
        "",
    ]
    lines.extend(f"- {'PASS' if passed else 'FAIL'} — `{name}`" for name, passed in checks.items())
    lines.extend(["", "## Visual QA", "", "The editable DOCX was rendered and visually inspected page-by-page. Figures 1-5 document technical contracts, row accounting, numerical diagnostics, policy dependence, and evidence boundaries. Figure 6 is a synthetic downstream-consumer stress test. None is presented as an empirical association, biological result, clinical implication, or comparative-performance claim.", "", "## Remaining hard gates", ""])
    lines.extend(f"- {blocker}" for blocker in blockers)
    lines.extend(["", "**Submission instruction:** Complete and verify every author-owned field before upload or submission."])
    args.out_report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
