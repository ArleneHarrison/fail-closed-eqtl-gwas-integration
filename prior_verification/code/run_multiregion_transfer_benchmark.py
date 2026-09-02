"""Run the locked chr11 multi-region workflow-only transfer benchmark.

The runner never calls a biological/statistical model.  It preserves every
selected eQTL row, records all reconciliation failures, builds unmodified EUR
dosage-correlation LD, and delegates the final decision to the existing
fail-closed preflight gate.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from pyliftover import LiftOver
from scipy.stats import beta as beta_distribution

from audit_1000g_eur_reference import classify_record
from hcsmr.reanalysis.preflight_gate import run_preflight


EXPECTED_SHA256 = {
    "qtd": "8666ab0c55109a09d7df8b4aa9129fac46300c559f71fbfba5602814cb4092c1",
    "cad": "90d3ef0a4a3783e374e35ae341383cf3df96b114d236299c366eb50c08695e97",
    "vcf": "64b9f63f6d3b5056f60cb297016d9da169e0668083620e4819c934259cd5b6d1",
    "chain": "14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1",
}
EXPECTED_BYTES = {
    "qtd": 3_161_468_217,
    "cad": 3_250_443_774,
    "vcf": 734_583_405,
    "chain": 1_246_411,
}
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, span = region.split(":", 1)
    start, end = map(int, span.split("-", 1))
    return chrom.removeprefix("chr"), start, end


def allele(value: Any) -> str:
    return str(value).upper()


def truth(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def transform_allele(value: str, strand: str) -> str:
    return value.translate(COMPLEMENT) if strand == "-" else value


def variant_key(chrom: Any, position: Any, ref: Any, alt: Any) -> str:
    return f"chr{str(chrom).removeprefix('chr')}_{position}_{allele(ref)}_{allele(alt)}"


def dosage(call: str, gt_index: int) -> float:
    values = call.split(":")
    genotype = values[gt_index] if gt_index < len(values) else "."
    if "." in genotype:
        return float("nan")
    try:
        return float(sum(int(value) for value in genotype.replace("|", "/").split("/")))
    except ValueError:
        return float("nan")


def exact_binomial_interval(successes: int, trials: int, alpha: float = 0.05) -> tuple[float, float]:
    if trials <= 0:
        return float("nan"), float("nan")
    lower = 0.0 if successes == 0 else float(beta_distribution.ppf(alpha / 2, successes, trials - successes + 1))
    upper = 1.0 if successes == trials else float(beta_distribution.ppf(1 - alpha / 2, successes + 1, trials - successes))
    return lower, upper


def verify_inputs(paths: dict[str, Path], tabix: str) -> dict[str, Any]:
    verification: dict[str, Any] = {}
    for label, path in paths.items():
        if label not in EXPECTED_SHA256:
            continue
        observed_bytes = path.stat().st_size
        observed_sha256 = sha256_file(path)
        verification[label] = {
            "path": str(path),
            "byte_count": observed_bytes,
            "expected_byte_count": EXPECTED_BYTES[label],
            "sha256": observed_sha256,
            "expected_sha256": EXPECTED_SHA256[label],
            "bytes_match": observed_bytes == EXPECTED_BYTES[label],
            "sha256_match": observed_sha256 == EXPECTED_SHA256[label],
        }
        if not verification[label]["bytes_match"] or not verification[label]["sha256_match"]:
            raise ValueError(f"immutable input verification failed for {label}")
    for label in ("qtd", "vcf", "chain"):
        completed = subprocess.run(["gzip", "-t", str(paths[label])], check=False)
        verification[label]["gzip_test_exit_code"] = completed.returncode
        if completed.returncode != 0:
            raise ValueError(f"gzip integrity failed for {label}")
    tabix_test = subprocess.run([tabix, "-l", str(paths["qtd"])], capture_output=True, text=True, check=True)
    verification["qtd"]["tabix_contigs_sha256"] = hashlib.sha256(tabix_test.stdout.encode()).hexdigest()
    vcf_index_test = subprocess.run([tabix, "-l", str(paths["vcf"])], capture_output=True, text=True, check=True)
    verification["vcf"]["tabix_contigs_sha256"] = hashlib.sha256(vcf_index_test.stdout.encode()).hexdigest()
    return verification


def qtd_header(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def extract_and_liftover_regions(
    *, qtd: Path, tabix: str, chain: Path, protocol: dict[str, Any], out_dir: Path
) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, str]]], list[str]]:
    header = qtd_header(qtd)
    required = {"gene_id", "molecular_trait_id", "chromosome", "position", "variant", "ref", "alt"}
    missing = required.difference(header)
    if missing:
        raise ValueError(f"QTD header missing required fields: {sorted(missing)}")
    width = len(header)
    converter = LiftOver(str(chain))
    region_states: list[dict[str, Any]] = []
    lifted_by_region: dict[str, list[dict[str, str]]] = {}
    for spec in protocol["regions"]:
        started = time.perf_counter()
        region_id, source_region = spec["region_id"], spec["source_region"]
        region_dir = out_dir / "regions" / region_id
        region_dir.mkdir(parents=True, exist_ok=True)
        chrom, start, end = parse_region(source_region)
        completed = subprocess.run(
            [tabix, str(qtd), source_region], capture_output=True, text=True, check=True
        )
        raw_lines = [line for line in completed.stdout.splitlines() if line]
        malformed = 0
        region_rows: list[dict[str, str]] = []
        for line in raw_lines:
            values = line.split("\t")
            if len(values) != width:
                malformed += 1
                continue
            row = dict(zip(header, values, strict=True))
            try:
                position = int(row["position"])
            except ValueError:
                malformed += 1
                continue
            if str(row["chromosome"]).removeprefix("chr") == chrom and start <= position <= end:
                region_rows.append(row)
        counts = Counter(row["gene_id"] for row in region_rows if row.get("gene_id"))
        eligible = sorted(gene for gene, count in counts.items() if count >= 2)
        state: dict[str, Any] = {
            "region_id": region_id,
            "source_region": source_region,
            "tabix_rows": len(raw_lines),
            "complete_region_rows": len(region_rows),
            "malformed_rows": malformed,
            "selected_gene_id": eligible[0] if eligible else None,
            "selection_rule": protocol["selection_rule"],
            "source_status": "selected" if eligible else "E_NO_ELIGIBLE_GENE",
        }
        if not eligible:
            state["runtime_extract_liftover_seconds"] = time.perf_counter() - started
            region_states.append(state)
            lifted_by_region[region_id] = []
            write_json(region_dir / "eqtl_extract_audit.json", state)
            continue
        selected_gene = eligible[0]
        selected = [dict(row) for row in region_rows if row["gene_id"] == selected_gene]
        write_tsv(region_dir / "eqtl_selected_gene.tsv", header, selected)
        lifted: list[dict[str, str]] = []
        liftover_counts: Counter[str] = Counter()
        target_positions: list[int] = []
        lifted_fields = header + [
            "liftover_status", "liftover_mapping_count", "target_chromosome",
            "target_position", "target_ref", "target_alt", "target_variant",
            "liftover_strand", "liftover_query_chromosome",
        ]
        for source_row in selected:
            row = dict(source_row)
            source_chrom = str(row["chromosome"]).removeprefix("chr")
            try:
                source_pos = int(row["position"])
            except ValueError:
                row.update(liftover_status="invalid_source_position", liftover_mapping_count="0")
                liftover_counts["invalid_source_position"] += 1
                lifted.append(row)
                continue
            query_chrom = f"chr{source_chrom}"
            mappings = converter.convert_coordinate(query_chrom, source_pos - 1) or []
            if not mappings:
                query_chrom = source_chrom
                mappings = converter.convert_coordinate(query_chrom, source_pos - 1) or []
            row["liftover_query_chromosome"] = query_chrom
            row["liftover_mapping_count"] = str(len(mappings))
            if len(mappings) != 1:
                status = "unmapped" if not mappings else "ambiguous_mapping"
                row["liftover_status"] = status
                liftover_counts[status] += 1
                lifted.append(row)
                continue
            target_chrom, target_zero, strand, _score = mappings[0]
            target_pos = int(target_zero) + 1
            target_chromosome = str(target_chrom).removeprefix("chr")
            target_ref = transform_allele(str(row["ref"]), strand)
            target_alt = transform_allele(str(row["alt"]), strand)
            row.update(
                liftover_status="mapped",
                target_chromosome=target_chromosome,
                target_position=str(target_pos),
                target_ref=target_ref,
                target_alt=target_alt,
                target_variant=variant_key(target_chromosome, target_pos, target_ref, target_alt),
                liftover_strand=strand,
            )
            liftover_counts["mapped"] += 1
            target_positions.append(target_pos)
            lifted.append(row)
        write_tsv(region_dir / "eqtl_lifted.tsv", lifted_fields, lifted)
        state.update(
            selected_gene_rows=len(selected),
            liftover_status_counts=dict(sorted(liftover_counts.items())),
            target_position_min=min(target_positions) if target_positions else None,
            target_position_max=max(target_positions) if target_positions else None,
            target_region=(
                f"11:{min(target_positions)}-{max(target_positions)}" if target_positions else None
            ),
            runtime_extract_liftover_seconds=time.perf_counter() - started,
        )
        write_json(region_dir / "eqtl_extract_audit.json", state)
        region_states.append(state)
        lifted_by_region[region_id] = lifted
    return region_states, lifted_by_region, header


def scan_cad_once(
    cad: Path, region_states: list[dict[str, Any]], out_dir: Path
) -> tuple[dict[str, list[dict[str, str]]], list[str], dict[str, Any]]:
    active = {
        state["region_id"]: (int(state["target_position_min"]), int(state["target_position_max"]))
        for state in region_states
        if state.get("target_position_min") is not None and state.get("target_position_max") is not None
    }
    rows_by_region: dict[str, list[dict[str, str]]] = {state["region_id"]: [] for state in region_states}
    scanned = malformed = 0
    started = time.perf_counter()
    with cad.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        required = {
            "chromosome", "base_pair_location", "markername", "effect_allele",
            "other_allele", "beta", "standard_error", "n",
        }
        missing = required.difference(header)
        if missing:
            raise ValueError(f"CAD header missing fields: {sorted(missing)}")
        width = len(header)
        chromosome_index = header.index("chromosome")
        position_index = header.index("base_pair_location")
        for values in reader:
            scanned += 1
            if len(values) != width:
                malformed += 1
                continue
            if str(values[chromosome_index]).removeprefix("chr") != "11":
                continue
            try:
                position = int(values[position_index])
            except ValueError:
                malformed += 1
                continue
            for region_id, (start, end) in active.items():
                if start <= position <= end:
                    rows_by_region[region_id].append(dict(zip(header, values, strict=True)))
    for state in region_states:
        region_id = state["region_id"]
        write_tsv(out_dir / "regions" / region_id / "cad_region.tsv", header, rows_by_region[region_id])
        state["cad_rows"] = len(rows_by_region[region_id])
    audit = {
        "source_rows_scanned_excluding_header": scanned,
        "malformed_width_or_position_rows": malformed,
        "active_target_regions": len(active),
        "single_complete_stream_scan": True,
        "runtime_seconds": time.perf_counter() - started,
    }
    write_json(out_dir / "cad_complete_scan_audit.json", audit)
    return rows_by_region, header, audit


def read_eur_panel(panel: Path) -> list[str]:
    rows = list(csv.DictReader(panel.open("r", encoding="utf-8", newline=""), delimiter="\t"))
    samples = [row["sample"] for row in rows if row["super_pop"] == "EUR"]
    if len(samples) != 503:
        raise ValueError(f"expected 503 EUR panel samples, found {len(samples)}")
    return samples


def indexed_vcf_region(
    *, tabix: str, vcf: Path, region: str, eur_samples: list[str], max_missing: float
) -> tuple[list[dict[str, Any]], dict[tuple[str, str, str, str], dict[str, Any]]]:
    process = subprocess.Popen(
        [tabix, "-h", str(vcf), region], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, encoding="utf-8",
    )
    if process.stdout is None:
        raise RuntimeError("tabix stdout unavailable")
    selected_indices: list[int] | None = None
    audits: list[dict[str, Any]] = []
    records: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    for line in process.stdout:
        if line.startswith("#CHROM"):
            all_samples = line.rstrip("\n").split("\t")[9:]
            index = {sample: offset for offset, sample in enumerate(all_samples)}
            missing = [sample for sample in eur_samples if sample not in index]
            if missing:
                process.kill()
                raise ValueError(f"{len(missing)} EUR samples absent from VCF")
            selected_indices = [index[sample] for sample in eur_samples]
            continue
        if line.startswith("#"):
            continue
        if selected_indices is None:
            process.kill()
            raise ValueError("VCF output lacks #CHROM header")
        fields = line.rstrip("\n").split("\t")
        audit = classify_record(fields, selected_indices, max_missing)
        audits.append(audit)
        key = (str(fields[0]).removeprefix("chr"), fields[1], allele(fields[3]), allele(fields[4]))
        if key not in records:
            records[key] = {"fields": fields, "audit": audit, "sample_indices": selected_indices}
        else:
            records[key]["duplicate_exact_vcf_record"] = True
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"tabix failed ({return_code}): {stderr.strip()}")
    return audits, records


def harmonize_region(
    *, lifted: list[dict[str, str]], cad_rows: list[dict[str, str]],
    vcf_records: dict[tuple[str, str, str, str], dict[str, Any]], out_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], Counter[str]]:
    gwas_by_coordinate: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in cad_rows:
        gwas_by_coordinate[(str(row["chromosome"]).removeprefix("chr"), row["base_pair_location"])].append(row)
    vcf_positions = {(key[0], key[1]) for key in vcf_records}
    audit_rows: list[dict[str, Any]] = []
    harmonized: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    for source in lifted:
        row: dict[str, Any] = dict(source)
        row.update(
            vcf_status="not_checked", ld_eligibility="not_checked",
            gwas_coordinate_row_count="0", gwas_markername="",
            gwas_effect_allele="", gwas_other_allele="",
            alignment_status="", retained_for_regional_input="False",
        )
        if row.get("liftover_status") != "mapped":
            status = "liftover_" + str(row.get("liftover_status"))
            row["alignment_status"] = status
            counts[status] += 1
            audit_rows.append(row)
            continue
        key = (
            str(row["target_chromosome"]).removeprefix("chr"), row["target_position"],
            allele(row["target_ref"]), allele(row["target_alt"]),
        )
        vcf = vcf_records.get(key)
        if vcf is None:
            row["vcf_status"] = (
                "position_present_alleles_mismatch" if key[:2] in vcf_positions else "position_absent"
            )
            status = "vcf_" + row["vcf_status"]
            row["alignment_status"] = status
            counts[status] += 1
            audit_rows.append(row)
            continue
        row["vcf_status"] = "exact_ref_alt_match"
        include_for_ld = bool(vcf["audit"]["include_for_ld"])
        row["ld_eligibility"] = "eligible" if include_for_ld else "ineligible:" + str(vcf["audit"]["exclusion_reason"])
        matches = gwas_by_coordinate.get((key[0], key[1]), [])
        row["gwas_coordinate_row_count"] = str(len(matches))
        if not matches:
            row["alignment_status"] = "gwas_coordinate_absent"
            counts["gwas_coordinate_absent"] += 1
            audit_rows.append(row)
            continue
        allele_matches: list[tuple[str, dict[str, str]]] = []
        for gwas in matches:
            if allele(row["target_alt"]) == allele(gwas["effect_allele"]) and allele(row["target_ref"]) == allele(gwas["other_allele"]):
                allele_matches.append(("aligned", gwas))
            elif allele(row["target_alt"]) == allele(gwas["other_allele"]) and allele(row["target_ref"]) == allele(gwas["effect_allele"]):
                allele_matches.append(("flipped", gwas))
        if len(allele_matches) != 1:
            status = "gwas_allele_no_exact_match" if not allele_matches else "gwas_allele_ambiguous_multiple_matches"
            row["alignment_status"] = status
            counts[status] += 1
            audit_rows.append(row)
            continue
        alignment, gwas = allele_matches[0]
        row.update(
            gwas_markername=gwas["markername"],
            gwas_effect_allele=gwas["effect_allele"],
            gwas_other_allele=gwas["other_allele"],
            alignment_status=alignment,
        )
        if not include_for_ld:
            counts["ld_ineligible_after_" + alignment] += 1
            audit_rows.append(row)
            continue
        raw_beta = float(gwas["beta"])
        oriented_beta = raw_beta if alignment == "aligned" else -raw_beta
        row["retained_for_regional_input"] = "True"
        counts[alignment] += 1
        audit_rows.append(row)
        harmonized.append({
            "eqtl_variant": row["variant"],
            "target_variant": row["target_variant"],
            "chromosome": row["target_chromosome"],
            "position": row["target_position"],
            "eqtl_effect_allele": row["target_alt"],
            "eqtl_other_allele": row["target_ref"],
            "eqtl_beta": row.get("beta", ""),
            "eqtl_se": row.get("se", ""),
            "eqtl_pvalue": row.get("pvalue", ""),
            "eqtl_an": row.get("an", ""),
            "gwas_markername": gwas["markername"],
            "gwas_effect_allele": gwas["effect_allele"],
            "gwas_other_allele": gwas["other_allele"],
            "gwas_beta_raw": gwas["beta"],
            "gwas_beta_aligned_to_eqtl": oriented_beta,
            "gwas_se": gwas["standard_error"],
            "gwas_pvalue": gwas.get("p_value", ""),
            "gwas_n": gwas["n"],
            "gwas_cases": gwas.get("cases", ""),
            "alignment_status": alignment,
            "vcf_status": "exact_ref_alt_match",
        })
    audit_fields = list(audit_rows[0]) if audit_rows else ["alignment_status"]
    harmonized_fields = list(harmonized[0]) if harmonized else [
        "eqtl_variant", "target_variant", "chromosome", "position", "eqtl_effect_allele",
        "eqtl_other_allele", "eqtl_beta", "eqtl_se", "eqtl_pvalue", "eqtl_an",
        "gwas_markername", "gwas_effect_allele", "gwas_other_allele", "gwas_beta_raw",
        "gwas_beta_aligned_to_eqtl", "gwas_se", "gwas_pvalue", "gwas_n", "gwas_cases",
        "alignment_status", "vcf_status",
    ]
    write_tsv(out_dir / "all_rows_audit.tsv", audit_fields, audit_rows)
    write_tsv(out_dir / "harmonized.tsv", harmonized_fields, harmonized)
    return audit_rows, harmonized, counts


def build_ld(
    *, harmonized: list[dict[str, Any]], vcf_records: dict[tuple[str, str, str, str], dict[str, Any]],
    out_dir: Path,
) -> tuple[np.ndarray, list[str], dict[str, float], dict[str, Any]]:
    requested = list(dict.fromkeys(str(row["target_variant"]) for row in harmonized))
    observed: dict[str, dict[str, Any]] = {}
    for key, record in vcf_records.items():
        target = variant_key(*key)
        if target not in requested or not bool(record["audit"]["include_for_ld"]):
            continue
        fields = record["fields"]
        fmt = fields[8].split(":")
        if "GT" not in fmt:
            continue
        gt_index = fmt.index("GT")
        values = np.array(
            [dosage(fields[9 + index], gt_index) for index in record["sample_indices"]], dtype=float
        )
        missing = int(np.isnan(values).sum())
        mean = float(np.nanmean(values)) if missing < len(values) else float("nan")
        if not np.isfinite(mean):
            continue
        imputed = values.copy()
        imputed[np.isnan(imputed)] = mean
        if not np.isfinite(imputed).all() or float(np.std(imputed)) == 0:
            continue
        maf = min(mean / 2.0, 1.0 - mean / 2.0)
        if not 0 < maf < 0.5:
            continue
        observed[target] = {"dosage": imputed, "maf": maf, "missing": missing}
    retained_mean = [variant for variant in requested if variant in observed]
    retained_complete = [variant for variant in retained_mean if observed[variant]["missing"] == 0]

    def correlation(variants: list[str]) -> np.ndarray:
        if not variants:
            return np.empty((0, 0), dtype=float)
        if len(variants) == 1:
            return np.array([[1.0]], dtype=float)
        return np.corrcoef(np.vstack([observed[variant]["dosage"] for variant in variants]))

    mean_ld = correlation(retained_mean)
    complete_ld = correlation(retained_complete)
    np.savez_compressed(out_dir / "ld_mean_impute.npz", variant_ids=np.array(retained_mean), ld=mean_ld)
    np.savez_compressed(out_dir / "ld_complete_genotype_only.npz", variant_ids=np.array(retained_complete), ld=complete_ld)
    variant_audit = []
    for variant in requested:
        record = observed.get(variant)
        variant_audit.append({
            "target_variant": variant,
            "mean_impute_status": "retained" if record else "not_usable_in_vcf_after_qc",
            "complete_genotype_status": "retained" if record and record["missing"] == 0 else "excluded",
            "eur_samples": 503,
            "missing_genotypes": record["missing"] if record else "",
            "maf": record["maf"] if record else "",
        })
    write_tsv(
        out_dir / "ld_variant_audit.tsv",
        ["target_variant", "mean_impute_status", "complete_genotype_status", "eur_samples", "missing_genotypes", "maf"],
        variant_audit,
    )
    maf = {variant: float(observed[variant]["maf"]) for variant in retained_mean}
    sensitivity = {
        "requested_unique_variants": len(requested),
        "mean_impute_retained_variants": len(retained_mean),
        "complete_genotype_only_retained_variants": len(retained_complete),
        "variants_with_realized_mean_imputation": sum(observed[variant]["missing"] > 0 for variant in retained_mean),
        "same_variant_order": retained_mean == retained_complete,
        "matrices_byte_equal": bool(
            retained_mean == retained_complete and mean_ld.shape == complete_ld.shape and np.array_equal(mean_ld, complete_ld)
        ),
    }
    write_json(out_dir / "ld_sensitivity.json", sensitivity)
    return mean_ld, retained_mean, maf, sensitivity


def run_region_gate(
    *, harmonized: list[dict[str, Any]], ld: np.ndarray, ld_variant_ids: list[str], maf: dict[str, float],
    provenance: list[dict[str, Any]], protocol: dict[str, Any], out_dir: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    for source in harmonized:
        row = dict(source)
        row["maf"] = maf.get(str(row["target_variant"]))
        rows.append(row)
    fieldnames = list(rows[0]) if rows else [
        "target_variant", "eqtl_se", "gwas_se", "maf", "eqtl_an", "gwas_n",
        "alignment_status", "vcf_status", "gwas_cases",
    ]
    write_tsv(out_dir / "preflight_summary.tsv", fieldnames, rows)
    result = run_preflight(
        rows=rows,
        ld=ld,
        ld_variant_ids=ld_variant_ids,
        provenance=provenance,
        eqtl_trait_type="quant",
        gwas_trait_type="cc",
        case_fraction=float(protocol["cad_cases"]) / float(protocol["cad_total_n"]),
        max_condition_number=float(protocol["max_condition_number"]),
    )
    write_json(out_dir / "preflight_status.json", result)
    return result, rows


def artifact_manifest(root: Path, command: list[str], verification: dict[str, Any]) -> dict[str, Any]:
    files = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "artifact_manifest.json":
            files.append({
                "path": str(path.relative_to(root)).replace(os.sep, "/"),
                "byte_count": path.stat().st_size,
                "sha256": sha256_file(path),
            })
    return {
        "schema_version": "1.0",
        "workflow_only": True,
        "command": command,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "input_verification": verification,
        "artifacts": files,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qtd", required=True)
    parser.add_argument("--qtd-index", required=True)
    parser.add_argument("--cad", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--chain", required=True)
    parser.add_argument("--provenance-json", required=True)
    parser.add_argument("--regions-json", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--tabix", default="tabix")
    parser.add_argument("--max-missing", type=float, default=0.05)
    args = parser.parse_args()
    started = time.perf_counter()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "qtd": Path(args.qtd), "cad": Path(args.cad), "vcf": Path(args.vcf),
        "chain": Path(args.chain), "qtd_index": Path(args.qtd_index), "panel": Path(args.panel),
    }
    protocol = json.loads(Path(args.regions_json).read_text(encoding="utf-8"))
    provenance = json.loads(Path(args.provenance_json).read_text(encoding="utf-8"))
    verification = verify_inputs(paths, args.tabix)
    verification["qtd_index"] = {
        "path": str(paths["qtd_index"]), "byte_count": paths["qtd_index"].stat().st_size,
        "sha256": sha256_file(paths["qtd_index"]),
    }
    verification["panel"] = {
        "path": str(paths["panel"]), "byte_count": paths["panel"].stat().st_size,
        "sha256": sha256_file(paths["panel"]),
    }
    verification["regions_protocol"] = {
        "path": str(Path(args.regions_json)), "sha256": sha256_file(Path(args.regions_json)),
    }
    write_json(out_dir / "input_verification.json", verification)
    region_states, lifted_by_region, _qtd_fields = extract_and_liftover_regions(
        qtd=paths["qtd"], tabix=args.tabix, chain=paths["chain"], protocol=protocol, out_dir=out_dir
    )
    cad_by_region, _cad_fields, cad_audit = scan_cad_once(paths["cad"], region_states, out_dir)
    eur_samples = read_eur_panel(paths["panel"])
    summary_rows: list[dict[str, Any]] = []
    for state in region_states:
        region_started = time.perf_counter()
        region_id = state["region_id"]
        region_dir = out_dir / "regions" / region_id
        if not state.get("target_region"):
            result = {
                "schema_version": "1.0", "status": "INELIGIBLE",
                "status_codes": [state.get("source_status", "E_NO_MAPPED_ROWS")],
                "failures": [{"code": state.get("source_status", "E_NO_MAPPED_ROWS"), "message": "no target region was available", "affected_rows": None}],
                "actions": ["Do not run a model or substitute a new region."],
            }
            write_json(region_dir / "preflight_status.json", result)
            summary_rows.append({
                **state, "harmonized_rows": 0, "unique_harmonized_variants": 0,
                "ld_variants": 0, "gate_status": "INELIGIBLE",
                "status_codes": ";".join(result["status_codes"]),
                "runtime_region_seconds": time.perf_counter() - region_started,
            })
            continue
        vcf_audit, vcf_records = indexed_vcf_region(
            tabix=args.tabix, vcf=paths["vcf"], region=state["target_region"],
            eur_samples=eur_samples, max_missing=args.max_missing,
        )
        vcf_fields = [
            "chrom", "position", "ref", "alt", "filter", "eur_n", "eur_missing_n",
            "eur_missing_fraction", "is_biallelic", "include_for_ld", "exclusion_reason",
        ]
        write_tsv(region_dir / "vcf_variant_audit.tsv", vcf_fields, vcf_audit)
        audit_rows, harmonized, alignment_counts = harmonize_region(
            lifted=lifted_by_region[region_id], cad_rows=cad_by_region[region_id],
            vcf_records=vcf_records, out_dir=region_dir,
        )
        ld, ld_variant_ids, maf, sensitivity = build_ld(
            harmonized=harmonized, vcf_records=vcf_records, out_dir=region_dir,
        )
        result, preflight_rows = run_region_gate(
            harmonized=harmonized, ld=ld, ld_variant_ids=ld_variant_ids, maf=maf,
            provenance=provenance, protocol=protocol, out_dir=region_dir,
        )
        diagnostics = result.get("ld_diagnostics", {})
        summary_rows.append({
            **state,
            "vcf_records": len(vcf_audit),
            "vcf_biallelic_eligible": sum(bool(row["include_for_ld"]) for row in vcf_audit),
            "alignment_status_counts": json.dumps(dict(sorted(alignment_counts.items())), sort_keys=True),
            "full_audit_rows": len(audit_rows),
            "harmonized_rows": len(harmonized),
            "unique_harmonized_variants": len({row["target_variant"] for row in harmonized}),
            "ld_variants": len(ld_variant_ids),
            "complete_genotype_only_variants": sensitivity["complete_genotype_only_retained_variants"],
            "realized_mean_imputation_variants": sensitivity["variants_with_realized_mean_imputation"],
            "gate_status": result["status"],
            "status_codes": ";".join(result["status_codes"]),
            "ld_min_eigenvalue": diagnostics.get("min_eigenvalue"),
            "ld_max_eigenvalue": diagnostics.get("max_eigenvalue"),
            "ld_condition_number": diagnostics.get("condition_number_2"),
            "ld_symmetry_max_error": diagnostics.get("symmetric_max_abs_error"),
            "ld_diagonal_max_error": diagnostics.get("diagonal_max_abs_error_from_one"),
            "runtime_region_seconds": time.perf_counter() - region_started,
        })
        print(json.dumps({
            "region_id": region_id, "selected_gene_id": state.get("selected_gene_id"),
            "status": result["status"], "status_codes": result["status_codes"],
            "harmonized_rows": len(preflight_rows), "ld_variants": len(ld_variant_ids),
        }, sort_keys=True), flush=True)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "benchmark_summary.tsv", sep="\t", index=False)
    attempts = len(summary_rows)
    ready = sum(row.get("gate_status") == "READY" for row in summary_rows)
    ineligible = attempts - ready
    ready_ci = exact_binomial_interval(ready, attempts)
    code_counts: Counter[str] = Counter()
    for row in summary_rows:
        for code in str(row.get("status_codes", "")).split(";"):
            if code:
                code_counts[code] += 1
    code_summary = {}
    for code, count in sorted(code_counts.items()):
        low, high = exact_binomial_interval(count, attempts)
        code_summary[code] = {"regions": count, "proportion": count / attempts, "exact_95_ci": [low, high]}
    aggregate = {
        "schema_version": "1.0",
        "scope": "QTD000216--GCST90132314--1000G-EUR chr11 workflow-only transfer benchmark",
        "attempted_regions": attempts,
        "ready_regions": ready,
        "ineligible_regions": ineligible,
        "ready_proportion": ready / attempts,
        "ready_exact_95_ci": list(ready_ci),
        "status_code_distribution": code_summary,
        "all_regions_retained": attempts == len(protocol["regions"]),
        "no_model_run": True,
        "cad_complete_scan": cad_audit,
        "total_runtime_seconds": time.perf_counter() - started,
        "interpretation_boundary": "Technical transfer results only; no biological, association, causal, clinical, or method-superiority inference.",
    }
    write_json(out_dir / "benchmark_summary.json", aggregate)
    command = [sys.executable, *sys.argv]
    write_json(out_dir / "artifact_manifest.json", artifact_manifest(out_dir, command, verification))
    print(json.dumps(aggregate, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
