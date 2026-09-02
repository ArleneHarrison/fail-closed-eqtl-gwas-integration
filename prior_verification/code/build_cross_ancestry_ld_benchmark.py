"""Build matched 1000 Genomes LD panels and quantify cross-ancestry LD transfer.

The script operates on a fixed, harmonized regional variant list.  It writes one
matrix per requested super-population, a common-variant matrix set, and an
auditable summary.  It does not infer association or biological effects.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"chr{str(chrom).removeprefix('chr')}_{pos}_{ref.upper()}_{alt.upper()}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_matrix(path: Path, variants: list[str], matrix: np.ndarray) -> None:
    frame = pd.DataFrame(matrix, columns=variants)
    frame.insert(0, "variant_id", variants)
    frame.to_csv(path, sep="\t", index=False, float_format="%.10g")


def matrix_diagnostics(matrix: np.ndarray) -> dict[str, float | int]:
    eigenvalues = np.linalg.eigvalsh(matrix)
    positive = eigenvalues[eigenvalues > 1e-10]
    return {
        "rank": int(np.linalg.matrix_rank(matrix, tol=1e-10)),
        "minimum_eigenvalue": float(eigenvalues.min()),
        "maximum_eigenvalue": float(eigenvalues.max()),
        "condition_number_positive_spectrum": (
            float(positive.max() / positive.min()) if len(positive) else float("inf")
        ),
    }


def main() -> None:
    import pysam

    parser = argparse.ArgumentParser()
    parser.add_argument("--harmonized", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--populations", nargs="+", default=["EUR", "AFR", "EAS", "SAS", "AMR"])
    parser.add_argument("--reference-population", default="EUR")
    parser.add_argument("--maf-min", type=float, default=0.01)
    parser.add_argument("--max-common-variants", type=int, default=100)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    harmonized_path = Path(args.harmonized)
    vcf_path = Path(args.vcf)
    panel_path = Path(args.panel)
    harmonized = pd.read_csv(harmonized_path, sep="\t")
    requested = list(dict.fromkeys(harmonized["target_variant"].astype(str)))
    requested_set = set(requested)
    panel = pd.read_csv(panel_path, sep="\t")

    chrom_text, interval = args.region.split(":", maxsplit=1)
    start_text, end_text = interval.split("-", maxsplit=1)
    chrom = chrom_text.removeprefix("chr")
    start, end = int(start_text), int(end_text)
    vcf = pysam.VariantFile(str(vcf_path))
    vcf_samples = set(vcf.header.samples)

    sample_map: dict[str, list[str]] = {}
    for population in args.populations:
        population_samples = panel.loc[panel["super_pop"] == population, "sample"].astype(str).tolist()
        overlap = [sample for sample in population_samples if sample in vcf_samples]
        if not overlap or len(overlap) != len(population_samples):
            raise ValueError(
                f"{population}: panel has {len(population_samples)} samples but VCF overlap is {len(overlap)}"
            )
        sample_map[population] = overlap

    # Read the indexed interval once and retain only requested biallelic records.
    records: dict[str, tuple[str, int, str, str]] = {}
    for record in vcf.fetch(chrom, start - 1, end):
        if len(record.alts or []) != 1:
            continue
        key = variant_key(chrom, record.pos, record.ref, record.alts[0])
        if key in requested_set:
            records[key] = (chrom, int(record.pos), str(record.ref), str(record.alts[0]))
    vcf.close()

    matrices: dict[str, tuple[list[str], np.ndarray]] = {}
    audit_rows: list[dict[str, object]] = []
    population_summaries: list[dict[str, object]] = []
    for population in args.populations:
        pop_vcf = pysam.VariantFile(str(vcf_path))
        pop_vcf.subset_samples(sample_map[population])
        observed: dict[str, tuple[np.ndarray, float, int]] = {}
        for record in pop_vcf.fetch(chrom, start - 1, end):
            if len(record.alts or []) != 1:
                continue
            key = variant_key(chrom, record.pos, record.ref, record.alts[0])
            if key not in requested_set:
                continue
            dosage = np.array(
                [
                    float(sum(call["GT"]))
                    if call["GT"] is not None and None not in call["GT"]
                    else np.nan
                    for call in record.samples.values()
                ],
                dtype=float,
            )
            observed_n = int(np.isfinite(dosage).sum())
            if observed_n < 3:
                continue
            mean = float(np.nanmean(dosage))
            maf = min(mean / 2.0, 1.0 - mean / 2.0)
            missing = int(np.isnan(dosage).sum())
            dosage[np.isnan(dosage)] = mean
            if maf < args.maf_min or float(np.std(dosage)) == 0.0:
                audit_rows.append(
                    {
                        "population": population,
                        "variant_id": key,
                        "status": "excluded_low_frequency_or_monomorphic",
                        "maf": maf,
                        "missing_genotypes_imputed": missing,
                    }
                )
                continue
            observed[key] = (dosage, maf, missing)
        pop_vcf.close()
        retained = [variant for variant in requested if variant in observed]
        if len(retained) < 2:
            raise ValueError(f"{population}: fewer than two usable variants")
        dosage_matrix = np.vstack([observed[variant][0] for variant in retained])
        ld = np.corrcoef(dosage_matrix)
        matrices[population] = (retained, ld)
        np.savez_compressed(out_dir / f"ld_{population}.npz", variant_ids=np.array(retained), ld=ld)
        write_matrix(out_dir / f"ld_{population}.tsv", retained, ld)
        diag = matrix_diagnostics(ld)
        population_summaries.append(
            {
                "population": population,
                "samples": len(sample_map[population]),
                "requested_variants": len(requested),
                "retained_variants": len(retained),
                **diag,
            }
        )
        retained_set = set(retained)
        for variant in requested:
            if variant in retained_set:
                _, maf, missing = observed[variant]
                audit_rows.append(
                    {
                        "population": population,
                        "variant_id": variant,
                        "status": "retained",
                        "maf": maf,
                        "missing_genotypes_imputed": missing,
                    }
                )
            elif variant not in records:
                audit_rows.append(
                    {
                        "population": population,
                        "variant_id": variant,
                        "status": "not_observed_as_requested_biallelic_record",
                        "maf": "",
                        "missing_genotypes_imputed": "",
                    }
                )

    common_set = set.intersection(*(set(variants) for variants, _ in matrices.values()))
    common_all = [variant for variant in requested if variant in common_set]
    if len(common_all) < 10:
        raise ValueError(f"only {len(common_all)} variants are common across populations")
    if len(common_all) > args.max_common_variants:
        selected_indices = np.linspace(0, len(common_all) - 1, args.max_common_variants, dtype=int)
        common = [common_all[index] for index in selected_indices]
    else:
        common = common_all
    pd.DataFrame({"variant_id": common}).to_csv(out_dir / "common_variants.tsv", sep="\t", index=False)

    common_matrices: dict[str, np.ndarray] = {}
    for population, (variants, matrix) in matrices.items():
        index = {variant: idx for idx, variant in enumerate(variants)}
        take = [index[variant] for variant in common]
        common_matrix = matrix[np.ix_(take, take)]
        common_matrices[population] = common_matrix
        write_matrix(out_dir / f"common_ld_{population}.tsv", common, common_matrix)
        np.savez_compressed(
            out_dir / f"common_ld_{population}.npz", variant_ids=np.array(common), ld=common_matrix
        )

    if args.reference_population not in common_matrices:
        raise ValueError("reference population is not among requested populations")
    reference = common_matrices[args.reference_population]
    transfer_rows: list[dict[str, object]] = []
    upper = np.triu_indices_from(reference, k=1)
    for population, matrix in common_matrices.items():
        delta = np.abs(matrix[upper] - reference[upper])
        informative = (np.abs(reference[upper]) >= 0.2) | (np.abs(matrix[upper]) >= 0.2)
        sign_discordance = (
            float(np.mean(np.sign(reference[upper][informative]) != np.sign(matrix[upper][informative])))
            if informative.any()
            else 0.0
        )
        transfer_rows.append(
            {
                "population": population,
                "reference_population": args.reference_population,
                "common_variants": len(common),
                "normalized_frobenius_distance": float(
                    np.linalg.norm(matrix - reference, ord="fro") / np.linalg.norm(reference, ord="fro")
                ),
                "median_absolute_delta_r": float(np.median(delta)),
                "p95_absolute_delta_r": float(np.quantile(delta, 0.95)),
                "max_absolute_delta_r": float(delta.max()),
                "sign_discordance_informative_pairs": sign_discordance,
                "informative_pairs": int(informative.sum()),
                **matrix_diagnostics(matrix),
            }
        )

    pd.DataFrame(audit_rows).to_csv(out_dir / "variant_audit.tsv", sep="\t", index=False)
    pd.DataFrame(population_summaries).to_csv(out_dir / "population_summary.tsv", sep="\t", index=False)
    pd.DataFrame(transfer_rows).to_csv(out_dir / "ld_transfer_summary.tsv", sep="\t", index=False)
    manifest = {
        "analysis_type": "cross-ancestry LD-panel portability benchmark",
        "claim_boundary": (
            "This benchmark compares LD references for the same fixed regional variants; "
            "it is not multi-ancestry eQTL or GWAS association validation."
        ),
        "region": args.region,
        "reference": "1000 Genomes Project Phase 3 GRCh37",
        "populations": args.populations,
        "reference_population": args.reference_population,
        "maf_threshold": args.maf_min,
        "requested_unique_variants": len(requested),
        "common_before_deterministic_thinning": len(common_all),
        "common_used_for_transfer_and_downstream_stress_test": len(common),
        "deterministic_selection": "evenly spaced indices in harmonized coordinate order",
        "inputs": {
            "harmonized_sha256": sha256(harmonized_path),
            "panel_sha256": sha256(panel_path),
            "vcf_sha256": sha256(vcf_path),
        },
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
