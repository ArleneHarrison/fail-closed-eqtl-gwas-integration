"""Build auditable EUR LD using a complete local gzip VCF scan.

This is the persistent fallback for environments without pysam/tabix.  It
streams the verified local VCF, uses exactly the panel-declared EUR samples,
and records every requested variant that is absent or fails the immutable VCF
eligibility contract.  It never accesses a partial input or modifies LD.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from audit_1000g_eur_reference import classify_record


def variant_key(chrom: str, position: str, ref: str, alt: str) -> str:
    return f"chr{str(chrom).removeprefix('chr')}_{position}_{ref.upper()}_{alt.upper()}"


def dosage(call: str, gt_index: int) -> float:
    fields = call.split(":")
    gt = fields[gt_index] if gt_index < len(fields) else "."
    if "." in gt:
        return float("nan")
    alleles = gt.replace("|", "/").split("/")
    try:
        return float(sum(int(value) for value in alleles))
    except ValueError:
        return float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--harmonized", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--out-npz", required=True)
    parser.add_argument("--variant-audit-out", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--missing-policy", choices=("mean_impute", "drop_variant"), required=True)
    parser.add_argument("--max-missing", type=float, default=0.05)
    args = parser.parse_args()

    requested = list(dict.fromkeys(pd.read_csv(args.harmonized, sep="\t")["target_variant"].astype(str)))
    requested_set = set(requested)
    panel = list(csv.DictReader(Path(args.panel).open(encoding="utf-8"), delimiter="\t"))
    eur = [row["sample"] for row in panel if row["super_pop"] == "EUR"]
    if len(eur) != 503:
        raise ValueError(f"expected 503 EUR samples, found {len(eur)}")
    chrom, span = args.region.split(":", 1)
    start, end = map(int, span.split("-", 1))
    selected_indices: list[int] | None = None
    observed: dict[str, tuple[np.ndarray, float, int]] = {}
    with gzip.open(args.vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                all_samples = line.rstrip("\n").split("\t")[9:]
                sample_index = {sample: index for index, sample in enumerate(all_samples)}
                absent = [sample for sample in eur if sample not in sample_index]
                if absent:
                    raise ValueError(f"EUR panel/VCF sample overlap is {len(eur) - len(absent)}, not 503")
                selected_indices = [sample_index[sample] for sample in eur]
                continue
            if line.startswith("#"):
                continue
            if selected_indices is None:
                raise ValueError("VCF lacks #CHROM header")
            fields = line.rstrip("\n").split("\t")
            if fields[0].removeprefix("chr") != chrom.removeprefix("chr"):
                continue
            position = int(fields[1])
            if position < start:
                continue
            if position > end:
                break
            if len(fields[4].split(",")) != 1:
                continue
            key = variant_key(fields[0], fields[1], fields[3], fields[4])
            if key not in requested_set:
                continue
            audit = classify_record(fields, selected_indices, args.max_missing)
            if not audit["include_for_ld"]:
                continue
            gt_index = fields[8].split(":").index("GT")
            calls = [fields[9 + index] for index in selected_indices]
            values = np.array([dosage(call, gt_index) for call in calls], dtype=float)
            missing = int(np.isnan(values).sum())
            if args.missing_policy == "drop_variant" and missing:
                continue
            mean = float(np.nanmean(values))
            values[np.isnan(values)] = mean
            if not np.isfinite(values).all() or np.std(values) == 0:
                continue
            maf = min(mean / 2.0, 1.0 - mean / 2.0)
            if not 0 < maf < 0.5:
                continue
            observed[key] = (values, maf, missing)
    retained = [key for key in requested if key in observed]
    if len(retained) < 2:
        raise ValueError("fewer than two requested variants have usable EUR genotypes")
    matrix = np.corrcoef(np.vstack([observed[key][0] for key in retained]))
    np.savez_compressed(args.out_npz, variant_ids=np.array(retained), ld=matrix)
    audit_path, manifest_path = Path(args.variant_audit_out), Path(args.manifest_out)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    with audit_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["target_variant", "status", "eur_samples", "imputed_missing_genotypes", "maf"])
        writer.writeheader()
        for key in requested:
            if key in observed:
                _, maf, missing = observed[key]
                writer.writerow({"target_variant": key, "status": f"retained_{args.missing_policy}", "eur_samples": 503, "imputed_missing_genotypes": missing, "maf": maf})
            else:
                writer.writerow({"target_variant": key, "status": "not_usable_in_vcf_after_qc", "eur_samples": 503, "imputed_missing_genotypes": "", "maf": ""})
    manifest = {
        "region": args.region, "access_method": "complete_gzip_stream_no_index", "index_not_used_reason": "persistent tabix/htslib executable unavailable",
        "reference": "1000 Genomes Phase 3 GRCh37", "super_population": "EUR", "eur_samples": 503,
        "requested_unique_variants": len(requested), "retained_variants": len(retained), "not_usable_variants": [key for key in requested if key not in observed],
        "max_missing_fraction": args.max_missing, "missing_genotype_handling": args.missing_policy,
        "ld_method": "Pearson correlation of dosage rows in retained-variant order", "status": "ld_complete",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
