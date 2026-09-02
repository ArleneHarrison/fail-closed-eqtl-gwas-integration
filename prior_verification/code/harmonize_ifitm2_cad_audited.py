"""Create a full retained-row eQTL--CAD harmonisation audit for the locked region."""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def truth(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def allele(value: str) -> str:
    """Use uppercase only for textual allele comparison; retain raw source fields."""
    return str(value).upper()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eqtl-lifted", required=True)
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--ld-audit", required=True)
    parser.add_argument("--out-audit", required=True)
    parser.add_argument("--out-harmonized", required=True)
    args = parser.parse_args()
    gwas_by_coordinate: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    with Path(args.gwas).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (str(row["chromosome"]).removeprefix("chr"), row["base_pair_location"])
            gwas_by_coordinate[key].append(row)
    ld_by_variant: dict[tuple[str, str, str, str], dict[str, str]] = {}
    with Path(args.ld_audit).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            ld_by_variant[(str(row["chrom"]).removeprefix("chr"), row["position"], allele(row["ref"]), allele(row["alt"]))] = row

    audit_path = Path(args.out_audit)
    harmonized_path = Path(args.out_harmonized)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    harmonized_path.parent.mkdir(parents=True, exist_ok=True)
    status_counts: Counter[str] = Counter()
    audit_fields: list[str] | None = None
    harmonized_fields = [
        "eqtl_variant", "target_variant", "chromosome", "position", "eqtl_effect_allele",
        "eqtl_other_allele", "eqtl_beta", "eqtl_se", "eqtl_pvalue", "eqtl_an",
        "gwas_markername", "gwas_effect_allele", "gwas_other_allele", "gwas_beta_aligned_to_eqtl",
        "gwas_se", "gwas_pvalue", "gwas_n", "gwas_cases", "alignment_status",
    ]
    with Path(args.eqtl_lifted).open(encoding="utf-8", newline="") as source, audit_path.open(
        "w", encoding="utf-8", newline=""
    ) as audit_handle, harmonized_path.open("w", encoding="utf-8", newline="") as harmonized_handle:
        reader = csv.DictReader(source, delimiter="\t")
        audit_fields = list(reader.fieldnames or []) + [
            "vcf_status", "ld_eligibility", "gwas_coordinate_row_count", "gwas_markername",
            "gwas_effect_allele", "gwas_other_allele", "alignment_status", "retained_for_coloc",
        ]
        audit_writer = csv.DictWriter(audit_handle, fieldnames=audit_fields, delimiter="\t")
        harmonized_writer = csv.DictWriter(harmonized_handle, fieldnames=harmonized_fields, delimiter="\t")
        audit_writer.writeheader()
        harmonized_writer.writeheader()
        for row in reader:
            row.update(vcf_status="not_checked", ld_eligibility="not_checked", gwas_coordinate_row_count="0",
                       gwas_markername="", gwas_effect_allele="", gwas_other_allele="",
                       alignment_status="", retained_for_coloc="False")
            if row.get("liftover_status") != "mapped":
                row["alignment_status"] = "liftover_" + str(row.get("liftover_status"))
                status_counts[row["alignment_status"]] += 1
                audit_writer.writerow(row)
                continue
            key = (row["target_chromosome"], row["target_position"], allele(row["target_ref"]), allele(row["target_alt"]))
            ld_row = ld_by_variant.get(key)
            if ld_row is None:
                position_present = any(k[:2] == key[:2] for k in ld_by_variant)
                row["vcf_status"] = "position_present_alleles_mismatch" if position_present else "position_absent"
                row["alignment_status"] = "vcf_" + row["vcf_status"]
                status_counts[row["alignment_status"]] += 1
                audit_writer.writerow(row)
                continue
            row["vcf_status"] = "exact_ref_alt_match"
            row["ld_eligibility"] = "eligible" if truth(ld_row["include_for_ld"]) else "ineligible:" + ld_row["exclusion_reason"]
            matches = gwas_by_coordinate.get((row["target_chromosome"], row["target_position"]), [])
            row["gwas_coordinate_row_count"] = str(len(matches))
            if not matches:
                row["alignment_status"] = "gwas_coordinate_absent"
                status_counts[row["alignment_status"]] += 1
                audit_writer.writerow(row)
                continue
            allele_matches: list[tuple[str, dict[str, str]]] = []
            for gwas in matches:
                if allele(row["target_alt"]) == allele(gwas["effect_allele"]) and allele(row["target_ref"]) == allele(gwas["other_allele"]):
                    allele_matches.append(("aligned", gwas))
                elif allele(row["target_alt"]) == allele(gwas["other_allele"]) and allele(row["target_ref"]) == allele(gwas["effect_allele"]):
                    allele_matches.append(("flipped", gwas))
            if len(allele_matches) != 1:
                row["alignment_status"] = "gwas_allele_no_exact_match" if not allele_matches else "gwas_allele_ambiguous_multiple_matches"
                status_counts[row["alignment_status"]] += 1
                audit_writer.writerow(row)
                continue
            alignment, gwas = allele_matches[0]
            row.update(gwas_markername=gwas["markername"], gwas_effect_allele=gwas["effect_allele"],
                       gwas_other_allele=gwas["other_allele"], alignment_status=alignment)
            if row["ld_eligibility"] != "eligible":
                status_counts["ld_ineligible_after_" + alignment] += 1
                audit_writer.writerow(row)
                continue
            beta = float(gwas["beta"])
            aligned_beta = beta if alignment == "aligned" else -beta
            row["retained_for_coloc"] = "True"
            status_counts[alignment] += 1
            audit_writer.writerow(row)
            harmonized_writer.writerow({
                "eqtl_variant": row["variant"], "target_variant": row["target_variant"],
                "chromosome": row["target_chromosome"], "position": row["target_position"],
                "eqtl_effect_allele": row["target_alt"], "eqtl_other_allele": row["target_ref"],
                "eqtl_beta": row["beta"], "eqtl_se": row["se"], "eqtl_pvalue": row["pvalue"], "eqtl_an": row["an"],
                "gwas_markername": gwas["markername"], "gwas_effect_allele": gwas["effect_allele"],
                "gwas_other_allele": gwas["other_allele"], "gwas_beta_aligned_to_eqtl": aligned_beta,
                "gwas_se": gwas["standard_error"], "gwas_pvalue": gwas["p_value"], "gwas_n": gwas["n"],
                "gwas_cases": gwas.get("cases", ""),
                "alignment_status": alignment,
            })
    print(";".join(f"{name}={count}" for name, count in sorted(status_counts.items())))


if __name__ == "__main__":
    main()
