"""Audit indexed 1000 Genomes regional records before constructing EUR LD.

This program deliberately writes an audit row for every tabix-returned VCF
record.  It does not calculate LD and cannot silently discard multiallelic,
non-PASS, or high-missingness records.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def classify_record(fields: list[str], selected_sample_indices: list[int], max_missing: float) -> dict[str, object]:
    """Classify one VCF record using only the specified sample columns."""
    chrom, position, _identifier, ref, alt, _qual, filt, _info, fmt = fields[:9]
    alt_alleles = alt.split(",")
    format_fields = fmt.split(":")
    gt_index = format_fields.index("GT") if "GT" in format_fields else None
    selected_calls = [fields[9 + i] for i in selected_sample_indices]
    missing = 0
    if gt_index is None:
        missing = len(selected_calls)
    else:
        for call in selected_calls:
            values = call.split(":")
            genotype = values[gt_index] if gt_index < len(values) else "."
            if "." in genotype:
                missing += 1
    missing_fraction = missing / len(selected_calls) if selected_calls else float("nan")
    reasons: list[str] = []
    if len(alt_alleles) != 1:
        reasons.append("multiallelic")
    if filt not in {"PASS", "."}:
        reasons.append(f"filter_{filt}")
    if gt_index is None:
        reasons.append("format_missing_GT")
    elif missing_fraction > max_missing:
        reasons.append("missingness_exceeds_threshold")
    return {
        "chrom": chrom,
        "position": position,
        "ref": ref,
        "alt": alt,
        "filter": filt,
        "eur_n": len(selected_calls),
        "eur_missing_n": missing,
        "eur_missing_fraction": missing_fraction,
        "is_biallelic": len(alt_alleles) == 1,
        "include_for_ld": not reasons,
        "exclusion_reason": ";".join(reasons),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tabix", required=True, help="tabix executable")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--region", required=True, help="indexed GRCh37 region, e.g. 11:1-1303655")
    parser.add_argument("--panel", required=True)
    parser.add_argument("--expected-eur-n", type=int, default=503)
    parser.add_argument("--max-missing", type=float, default=0.05)
    parser.add_argument("--variant-audit-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--selected-samples-out", required=True)
    args = parser.parse_args()

    panel_rows = list(csv.DictReader(Path(args.panel).open(encoding="utf-8"), delimiter="\t"))
    eur_samples = [row["sample"] for row in panel_rows if row["super_pop"] == "EUR"]
    if len(eur_samples) != args.expected_eur_n:
        raise ValueError(f"expected {args.expected_eur_n} EUR panel samples, found {len(eur_samples)}")

    result = subprocess.run(
        [args.tabix, "-h", args.vcf, args.region], capture_output=True, text=True, check=True
    )
    header = next((line for line in result.stdout.splitlines() if line.startswith("#CHROM")), None)
    if header is None:
        raise ValueError("tabix output lacks a VCF #CHROM header")
    vcf_samples = header.split("\t")[9:]
    sample_to_index = {sample: index for index, sample in enumerate(vcf_samples)}
    missing_panel_samples = [sample for sample in eur_samples if sample not in sample_to_index]
    if missing_panel_samples:
        raise ValueError(f"EUR panel samples missing from VCF: {len(missing_panel_samples)}")
    selected_indices = [sample_to_index[sample] for sample in eur_samples]

    rows = []
    for line in result.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9 + len(vcf_samples):
            raise ValueError("tabix returned malformed VCF record")
        rows.append(classify_record(fields, selected_indices, args.max_missing))

    for target in (args.variant_audit_out, args.summary_out, args.selected_samples_out):
        Path(target).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.variant_audit_out).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [
            "chrom", "position", "ref", "alt", "filter", "eur_n", "eur_missing_n",
            "eur_missing_fraction", "is_biallelic", "include_for_ld", "exclusion_reason",
        ])
        writer.writeheader()
        writer.writerows(rows)
    with Path(args.selected_samples_out).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "super_pop"])
        writer.writerows((sample, "EUR") for sample in eur_samples)
    summary = {
        "region": args.region,
        "eur_samples": len(eur_samples),
        "records_returned": len(rows),
        "biallelic_records": sum(bool(row["is_biallelic"]) for row in rows),
        "multiallelic_records": sum(not bool(row["is_biallelic"]) for row in rows),
        "records_eligible_for_ld": sum(bool(row["include_for_ld"]) for row in rows),
        "records_excluded": sum(not bool(row["include_for_ld"]) for row in rows),
        "max_missing_fraction": args.max_missing,
        "missing_genotype_handling": "records above threshold are retained in audit and excluded; no imputation occurs in this audit",
        "status": "audit_complete_no_ld_calculated",
    }
    with Path(args.summary_out).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary))
        writer.writeheader()
        writer.writerow(summary)


if __name__ == "__main__":
    main()
