#!/usr/bin/env python3
"""Prepare all frozen external-validation units directly from immutable raw data.

The implementation is deliberately fail closed.  It does not strand-flip,
guess identifiers, choose among discordant duplicate keys, condition LD, or
impute any summary-statistic field.  The only imputation is the prospectively
declared per-variant mean-dosage handling for reference-panel genotypes after a
2% missingness audit; complete-case LD is written as a required sensitivity.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

try:  # Kept optional at import time so pure helper tests run on build hosts.
    import certifi
    import pysam
    from pyliftover import LiftOver
except ImportError:  # pragma: no cover - checked explicitly by main()
    certifi = None
    pysam = None
    LiftOver = None


TISSUES = ("QTD000131", "QTD000136", "QTD000251", "QTD000256")
OUTCOMES = ("CAD", "HF")
SUMMARY_FIELDS = (
    "target_variant", "position", "eqtl_beta", "eqtl_se", "gwas_beta",
    "gwas_se", "maf", "eqtl_maf", "gwas_maf", "ld_maf", "eqtl_n",
    "gwas_n", "gwas_cases", "gwas_case_fraction_row_audit",
    "gwas_case_fraction_locked_study", "alignment_status",
    "vcf_status", "gene_id", "tissue_id", "outcome_id",
)
VCF_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5()  # noqa: S324 - provider-supplied integrity checksum, not security use
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json_sha256(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def ordered_hash(values: Iterable[str]) -> str:
    return hashlib.sha256(("\n".join(values) + "\n").encode("utf-8")).hexdigest()


def canonical_row(row: Mapping[str, Any]) -> str:
    return json.dumps(dict(row), sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def register_unique(
    registry: dict[Any, Mapping[str, Any]], key: Any, row: Mapping[str, Any], *, label: str
) -> str:
    """Register a key, collapsing exact duplicates and stopping on disagreement."""
    if key not in registry:
        registry[key] = row
        return "new"
    if canonical_row(registry[key]) == canonical_row(row):
        return "exact_duplicate"
    raise RuntimeError(f"E_DISCORDANT_DUPLICATE_KEY:{label}:{key}")


def register_eqtl_analytic(
    registry: dict[Any, Mapping[str, Any]], key: Any, row: Mapping[str, Any], *,
    label: str, annotation_audit: list[dict[str, Any]],
) -> str:
    """Collapse only analytically identical rows; retain rsID aliases in audit."""
    analytical = {field: value for field, value in row.items() if field != "rsid"}
    analytical["_annotation_rsid"] = str(row.get("rsid", ""))
    if key not in registry:
        registry[key] = analytical
        return "new"
    old_analysis = {field: value for field, value in registry[key].items() if field != "_annotation_rsid"}
    new_analysis = {field: value for field, value in analytical.items() if field != "_annotation_rsid"}
    if canonical_row(old_analysis) != canonical_row(new_analysis):
        raise RuntimeError(f"E_DISCORDANT_DUPLICATE_KEY:{label}:{key}")
    aliases = sorted(set(str(registry[key].get("_annotation_rsid", "")).split(";")) | {str(row.get("rsid", ""))} - {""})
    registry[key] = {**old_analysis, "_annotation_rsid": ";".join(aliases)}
    annotation_audit.append({
        "source": label,
        "analytic_key": repr(key),
        "annotation_field": "rsid",
        "annotation_values": ";".join(aliases),
        "decision": "collapsed_analytic_exact_duplicate_annotation_alias",
    })
    return "analytic_exact_duplicate"


def validate_protocol_freeze(protocol_path: Path, freeze_path: Path) -> dict[str, Any]:
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    checks = {
        "state": freeze.get("state") == "READ_ONLY_FROZEN",
        "protocol_id": freeze.get("protocol_id") == protocol.get("protocol_id"),
        "frozen_at_utc": freeze.get("frozen_at_utc") == protocol.get("frozen_at_utc"),
        "byte_count": int(freeze.get("protocol_byte_count", -1)) == protocol_path.stat().st_size,
        "file_sha256": freeze.get("protocol_file_sha256") == sha256_file(protocol_path),
        "canonical_sha256": freeze.get("protocol_canonical_json_sha256") == canonical_json_sha256(protocol),
    }
    failed = [key for key, passed in checks.items() if not passed]
    if failed:
        raise RuntimeError("E_PROTOCOL_FREEZE_MISMATCH:" + ",".join(failed))
    return protocol


def validate_metadata_lock(metadata: Mapping[str, Any], protocol: Mapping[str, Any]) -> dict[str, float]:
    if metadata.get("state") != "LOCKED" or metadata.get("protocol_canonical_json_sha256") != canonical_json_sha256(protocol):
        raise RuntimeError("E_METADATA_NOT_LOCKED")
    expected = {row["outcome_id"]: row["dataset_id"] for row in protocol["outcomes"]}
    fractions = {}
    for outcome in OUTCOMES:
        record = metadata.get("outcomes", {}).get(outcome, {})
        try:
            cases, controls, sample_size = int(record["cases"]), int(record["controls"]), int(record["sample_size"])
            fraction = float(record["case_fraction"])
        except (KeyError, TypeError, ValueError) as error:
            raise RuntimeError(f"E_METADATA_FIELDS_INVALID:{outcome}") from error
        if record.get("dataset_id") != expected[outcome] or cases <= 0 or controls <= 0 or cases + controls != sample_size:
            raise RuntimeError(f"E_METADATA_COUNTS_INVALID:{outcome}")
        if not np.isclose(fraction, cases / sample_size, rtol=0, atol=1e-15):
            raise RuntimeError(f"E_METADATA_CASE_FRACTION_INVALID:{outcome}")
        for field in ("ancestry", "genome_build", "metadata_url", "source_citation"):
            if not record.get(field):
                raise RuntimeError(f"E_METADATA_FIELDS_INVALID:{outcome}:{field}")
        fractions[outcome] = fraction
    return fractions


def validate_regions(
    regions: list[dict[str, str]], regions_path: Path, audit_path: Path,
    freeze_path: Path, observed_source_hashes: Mapping[str, str], protocol: Mapping[str, Any],
) -> None:
    if len(regions) != 24 or len({row["region_id"] for row in regions}) != 24:
        raise RuntimeError(f"E_REGION_REGISTRY_COUNT:{len(regions)}")
    for outcome in OUTCOMES:
        subset = [row for row in regions if row.get("anchor_outcome") == outcome]
        if len(subset) != 12:
            raise RuntimeError(f"E_REGION_OUTCOME_COUNT:{outcome}:{len(subset)}")
        if {int(row["selection_rank"]) for row in subset} != set(range(1, 13)):
            raise RuntimeError(f"E_REGION_RANK_SET:{outcome}")
        for row in subset:
            lead, start, end = int(row["lead_position"]), int(row["region_start"]), int(row["region_end"])
            if row["selection_basis"] != "GWAS_chromosome_position_and_p_value_only":
                raise RuntimeError("E_REGION_SELECTION_BASIS")
            if start != max(1, lead - 500_000) or end != lead + 500_000 or float(row["lead_p_value"]) > 5e-8:
                raise RuntimeError(f"E_REGION_WINDOW_OR_THRESHOLD:{row['region_id']}")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("n_regions") != 24 or audit.get("selection_basis") != "GWAS_chromosome_position_and_p_value_only":
        raise RuntimeError("E_REGION_AUDIT_INVALID")
    selection_inputs = audit.get("selection_inputs", {})
    cad_audit_hash = audit.get("CAD_sha256", selection_inputs.get("CAD", {}).get("sha256"))
    hf_audit_hash = audit.get("HF_sha256", selection_inputs.get("HF", {}).get("sha256"))
    if cad_audit_hash != observed_source_hashes["GCST005195"] or hf_audit_hash != observed_source_hashes["HERMES2_EUR_Pheno1"]:
        raise RuntimeError("E_REGION_SOURCE_HASH_MISMATCH")
    if audit.get("downstream_results_consulted") is not False:
        raise RuntimeError("E_REGION_AUDIT_DOWNSTREAM_FLAG")
    if audit.get("protocol_sha256") not in (None, canonical_json_sha256(protocol)):
        raise RuntimeError("E_REGION_AUDIT_PROTOCOL_MISMATCH")
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    freeze_checks = {
        "state": freeze.get("state") == "READ_ONLY_FROZEN",
        "regions_sha256": freeze.get("regions_sha256") == sha256_file(regions_path),
        "audit_sha256": freeze.get("audit_sha256") == sha256_file(audit_path),
        "protocol_sha256": freeze.get("protocol_canonical_json_sha256") == canonical_json_sha256(protocol),
        "n_regions": int(freeze.get("n_regions", -1)) == 24,
    }
    failed = [field for field, passed in freeze_checks.items() if not passed]
    if failed:
        raise RuntimeError("E_REGION_FREEZE_MISMATCH:" + ",".join(failed))


def verify_manifest_assets(
    root: Path, assets: dict[str, dict[str, Any]], regions: list[dict[str, str]]
) -> list[dict[str, Any]]:
    ids = set(TISSUES) | {f"{item}_tbi" for item in TISSUES} | {
        "GCST005195", "HERMES2_EUR_outer", "HERMES2_EUR_Pheno1",
        "UCSC_hg19ToHg38_chain", "UCSC_hg38ToHg19_chain", "1000G_panel",
    }
    ids |= {f"1000G_chr{normalize_chromosome(row['chromosome'])}_tbi_cache" for row in regions}
    records = []
    for asset_id in sorted(ids):
        if asset_id not in assets:
            raise RuntimeError(f"E_MANIFEST_ASSET_MISSING:{asset_id}")
        asset = assets[asset_id]
        path = root / asset["project_relative_path"]
        if not path.is_file() or path.stat().st_size != int(asset["byte_size"]):
            raise RuntimeError(f"E_MANIFEST_SIZE_MISMATCH:{asset_id}")
        observed = sha256_file(path)
        expected = str(asset.get("sha256") or "")
        if expected and observed != expected:
            raise RuntimeError(f"E_MANIFEST_SHA256_MISMATCH:{asset_id}")
        if asset_id == "HERMES2_EUR_Pheno1":
            md5_path = path.with_name(path.name + ".md5")
            expected_md5 = md5_path.read_text(encoding="utf-8").split()[0]
            if md5_file(path) != expected_md5:
                raise RuntimeError("E_HERMES_MEMBER_MD5_MISMATCH")
            asset["sha256"] = observed
        records.append({"asset_id": asset_id, "path": str(path), "byte_count": path.stat().st_size, "sha256": observed, "status": "PASS"})
    hf_tbi = root / "data/hf/primary/FORMAT-METAL_Pheno1_EUR.tsv.gz.tbi"
    if not hf_tbi.is_file():
        raise RuntimeError("E_HERMES_TABIX_INDEX_MISSING")
    runtime_freeze = json.loads((root / "config/runtime_indices.freeze.json").read_text(encoding="utf-8"))
    expected_hf_tbi = runtime_freeze.get("indices", {}).get("HERMES2_EUR_Pheno1_tbi", {})
    observed_hf_tbi_sha = sha256_file(hf_tbi)
    if (
        runtime_freeze.get("state") != "READ_ONLY_FROZEN"
        or expected_hf_tbi.get("project_relative_path") != "data/hf/primary/FORMAT-METAL_Pheno1_EUR.tsv.gz.tbi"
        or int(expected_hf_tbi.get("byte_count", -1)) != hf_tbi.stat().st_size
        or expected_hf_tbi.get("sha256") != observed_hf_tbi_sha
    ):
        raise RuntimeError("E_HERMES_TABIX_INDEX_FREEZE_MISMATCH")
    records.append({"asset_id": "HERMES2_EUR_Pheno1_tbi_runtime", "path": str(hf_tbi), "byte_count": hf_tbi.stat().st_size, "sha256": observed_hf_tbi_sha, "status": "PASS_FROZEN_RUNTIME_INDEX_HASH"})
    return records


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")


def write_json_atomic(path: Path, payload: Any) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_runtime() -> None:
    missing = []
    if certifi is None:
        missing.append("certifi")
    if pysam is None:
        missing.append("pysam")
    if LiftOver is None:
        missing.append("pyliftover")
    if missing:
        raise RuntimeError("E_RUNTIME_DEPENDENCY_MISSING:" + ",".join(missing))


def normalize_chromosome(value: Any) -> str:
    text = str(value).strip().removeprefix("chr").lstrip("0")
    return text or "0"


class StrictRoundTrip:
    """Unique, positive-strand, exact point round trips between GRCh37/38."""

    def __init__(self, chain_37_to_38: Path, chain_38_to_37: Path) -> None:
        self.forward = LiftOver(str(chain_37_to_38))
        self.reverse = LiftOver(str(chain_38_to_37))

    @staticmethod
    def _hits(lifter: Any, chrom: str, pos_1based: int) -> list[tuple[Any, ...]]:
        primary = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
        return [
            hit for hit in lifter.convert_coordinate("chr" + normalize_chromosome(chrom), pos_1based - 1)
            if hit[0] in primary
        ]

    def point(self, chrom: str, pos: int, direction: str) -> tuple[str, int] | None:
        if direction == "37to38":
            first, back = self.forward, self.reverse
        elif direction == "38to37":
            first, back = self.reverse, self.forward
        else:
            raise ValueError(direction)
        source_chrom = normalize_chromosome(chrom)
        hits = self._hits(first, source_chrom, pos)
        if len(hits) != 1 or hits[0][2] != "+":
            return None
        target_chrom, target_pos = normalize_chromosome(hits[0][0]), int(hits[0][1]) + 1
        reverse_hits = self._hits(back, target_chrom, target_pos)
        if len(reverse_hits) != 1 or reverse_hits[0][2] != "+":
            return None
        if normalize_chromosome(reverse_hits[0][0]) != source_chrom or int(reverse_hits[0][1]) + 1 != pos:
            return None
        return target_chrom, target_pos

    def window_37_to_38(self, chrom: str, start: int, end: int) -> tuple[str, int, int]:
        left = self.point(chrom, start, "37to38")
        right = self.point(chrom, end, "37to38")
        if left is None or right is None or left[0] != right[0] or right[1] < left[1]:
            raise RuntimeError(f"E_WINDOW_LIFTOVER:{chrom}:{start}-{end}")
        return left[0], left[1], right[1]


class TabixReader:
    def __init__(self, path: Path) -> None:
        self.path = path
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            self.header = handle.readline().rstrip("\n").split("\t")
        self.handle = pysam.TabixFile(str(path))

    def fetch(self, chrom: str, start: int, end: int) -> list[dict[str, str]]:
        rows = []
        for line in self.handle.fetch(normalize_chromosome(chrom), start - 1, end):
            values = line.rstrip("\n").split("\t")
            if len(values) != len(self.header):
                raise RuntimeError(f"E_TABIX_ROW_WIDTH:{self.path}:{len(values)}:{len(self.header)}")
            row = dict(zip(self.header, values, strict=True))
            chrom_field = "chromosome" if "chromosome" in row else "chr"
            pos_field = "position" if "position" in row else "pos_b37"
            try:
                observed_chrom, observed_pos = normalize_chromosome(row[chrom_field]), int(row[pos_field])
            except (KeyError, TypeError, ValueError) as error:
                raise RuntimeError(f"E_TABIX_COORDINATE_FIELDS:{self.path}") from error
            if observed_chrom != normalize_chromosome(chrom) or not start <= observed_pos <= end:
                raise RuntimeError(f"E_TABIX_INDEX_ROW_OUTSIDE_QUERY:{self.path}:{observed_chrom}:{observed_pos}")
            rows.append(row)
        return rows

    def close(self) -> None:
        self.handle.close()


def finite(value: Any) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def complete_eqtl(row: Mapping[str, str]) -> bool:
    required = ("gene_id", "chromosome", "position", "variant", "ref", "alt", "beta", "se", "an", "maf")
    if any(row.get(field, "") in {"", "NA", "NaN", "nan"} for field in required):
        return False
    if not all(finite(row[field]) for field in ("position", "beta", "se", "an", "maf")):
        return False
    an = int(float(row["an"]))
    return an > 0 and an % 2 == 0 and float(row["se"]) > 0 and 0 < float(row["maf"]) <= 0.5


def parse_alleles(first: str, second: str) -> frozenset[str] | None:
    a, b = first.upper(), second.upper()
    allowed = set("ACGT")
    if not a or not b or a == b or not set(a) <= allowed or not set(b) <= allowed:
        return None
    return frozenset((a, b))


def reference_record_position_decision(
    record_chrom: str, record_pos: int, record_stop: int, query_chrom: str, start: int, end: int
) -> str:
    if normalize_chromosome(record_chrom) != normalize_chromosome(query_chrom):
        raise RuntimeError(f"E_REMOTE_REFERENCE_CONTIG_MISMATCH:{record_chrom}:{record_pos}")
    if start <= record_pos <= end:
        return "include"
    if record_pos < start and record_stop >= start:
        return "exclude_interval_overlap_outside_POS_window"
    raise RuntimeError(f"E_REMOTE_REFERENCE_RECORD_OUTSIDE_QUERY:{record_chrom}:{record_pos}")


def register_reference_key(
    registry: dict[Any, dict[str, Any]], ambiguous: set[Any], details: dict[Any, list[dict[str, Any]]],
    key: Any, row: dict[str, Any],
) -> str:
    """Maintain a unique reference map with permanent discordance quarantine."""
    descriptor = row["audit_descriptor"]
    if key in ambiguous:
        details[key].append(descriptor)
        return "ambiguous_additional_record_quarantined"
    if key not in registry:
        registry[key] = row
        return "new"
    old = registry[key]
    identical = (
        old["ref"] == row["ref"] and old["alt"] == row["alt"]
        and np.array_equal(old["dosage"], row["dosage"], equal_nan=True)
    )
    if identical:
        return "exact_duplicate"
    registry.pop(key)
    ambiguous.add(key)
    details[key] = [old["audit_descriptor"], descriptor]
    return "discordant_key_permanently_quarantined"


def scan_cad_once(
    path: Path, regions: list[dict[str, str]], sentinel_path: Path, verified_sha256: str
) -> tuple[dict[str, list[dict[str, str]]], dict[str, Any]]:
    intervals: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for region in regions:
        intervals[normalize_chromosome(region["chromosome"])].append(
            (int(region["region_start"]), int(region["region_end"]), region["region_id"])
        )
    allocated = {region["region_id"]: [] for region in regions}
    source_rows = assigned_rows = malformed = 0
    expected_sentinels = {
        "rs11806316": ("G", 0.0354124),
        "rs60154123": ("C", -0.0438974),
        "rs699": ("A", -0.028766),
    }
    sentinel_rows = read_tsv(sentinel_path)
    if {row.get("variant_id") for row in sentinel_rows} != set(expected_sentinels) or any(
        row.get("status") != "PASS_EFFECT_DIRECTION_CONSISTENT" for row in sentinel_rows
    ):
        raise RuntimeError("E_CAD_SIGNED_SENTINEL_AUDIT_INVALID")
    observed_sentinels: dict[str, tuple[str, float]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        required = {"chr", "bp", "a1", "a2", "beta", "se", "pval", "N", "af", "uniqid"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise RuntimeError("E_CAD_HEADER_MISMATCH")
        for row in reader:
            source_rows += 1
            try:
                chrom, pos = normalize_chromosome(row["chr"]), int(row["bp"])
            except (TypeError, ValueError):
                malformed += 1
                continue
            old_id = str(row.get("oldID", ""))
            if old_id in expected_sentinels:
                observed_sentinels[old_id] = (row["a1"].upper(), float(row["beta"]))
            for start, end, region_id in intervals.get(chrom, ()):
                if start <= pos <= end:
                    allocated[region_id].append(row)
                    assigned_rows += 1
    for variant, (allele, beta) in expected_sentinels.items():
        observed = observed_sentinels.get(variant)
        if observed is None or observed[0] != allele or not np.isclose(observed[1], beta, rtol=0, atol=5e-10):
            raise RuntimeError(f"E_CAD_SIGNED_SENTINEL_SOURCE_MISMATCH:{variant}:{observed}")
    return allocated, {
        "scan_count": 1,
        "source_rows": source_rows,
        "assigned_region_rows_including_overlap": assigned_rows,
        "malformed_coordinate_rows": malformed,
        "source_sha256": verified_sha256,
        "signed_sentinel_audit_sha256": sha256_file(sentinel_path),
        "signed_sentinels_verified_against_source": sorted(observed_sentinels),
    }


def eur_samples(panel_path: Path) -> list[str]:
    with panel_path.open(encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        samples = [row["sample"] for row in rows if row.get("super_pop") == "EUR"]
    if len(samples) != 503 or len(samples) != len(set(samples)):
        raise RuntimeError(f"E_EUR_SAMPLE_REGISTRY:{len(samples)}")
    return samples


def fetch_reference_region(
    chrom: str, start: int, end: int, samples: list[str], verified_local_index: Path
) -> tuple[dict[tuple[int, frozenset[str]], dict[str, Any]], set[tuple[int, frozenset[str]]], dict[str, Any]]:
    url = VCF_URL.format(chrom=normalize_chromosome(chrom))
    os.environ["CURL_CA_BUNDLE"] = certifi.where()
    os.environ["SSL_CERT_FILE"] = certifi.where()
    variants: dict[tuple[int, frozenset[str]], dict[str, Any]] = {}
    ambiguous_keys: set[tuple[int, frozenset[str]]] = set()
    ambiguous_details: dict[tuple[int, frozenset[str]], list[dict[str, Any]]] = {}
    exact_duplicates = 0
    overlapping_records_excluded_by_pos = 0
    with pysam.VariantFile(url, index_filename=str(verified_local_index)) as source:
        absent = [sample for sample in samples if sample not in source.header.samples]
        if absent:
            raise RuntimeError(f"E_REFERENCE_SAMPLES_ABSENT:{len(absent)}")
        for record in source.fetch(normalize_chromosome(chrom), start - 1, end):
            position_decision = reference_record_position_decision(
                record.contig, int(record.pos), int(record.stop or record.pos), chrom, start, end
            )
            if position_decision != "include":
                overlapping_records_excluded_by_pos += 1
                continue
            if len(record.alts or ()) != 1:
                continue
            ref, alt = str(record.ref).upper(), str(record.alts[0]).upper()
            allele_key = parse_alleles(ref, alt)
            if allele_key is None:
                continue
            dosage = np.full(len(samples), np.nan, dtype=np.float64)
            for index, sample in enumerate(samples):
                gt = record.samples[sample].get("GT")
                if gt is not None and len(gt) == 2 and all(value is not None and value in (0, 1) for value in gt):
                    dosage[index] = float(gt[0] + gt[1])
            dosage_codes = np.where(np.isnan(dosage), 255, dosage).astype(np.uint8).tobytes()
            descriptor = {
                "ref": ref, "alt": alt, "rsid": str(record.id or ""),
                "AC": list(record.info.get("AC", ())), "AF": list(record.info.get("AF", ())),
                "dosage_sha256": hashlib.sha256(dosage_codes).hexdigest(),
            }
            row = {"position": int(record.pos), "ref": ref, "alt": alt, "rsid": str(record.id or ""), "dosage": dosage, "audit_descriptor": descriptor}
            key = (int(record.pos), allele_key)
            decision = register_reference_key(variants, ambiguous_keys, ambiguous_details, key, row)
            exact_duplicates += decision == "exact_duplicate"
    content_digest = hashlib.sha256()
    content_bytes = 0
    for (position, _), row in sorted(variants.items(), key=lambda item: (item[0][0], sorted(item[0][1]))):
        prefix = f"{position}\t{row['ref']}\t{row['alt']}\t".encode("ascii")
        encoded_dosage = np.where(np.isnan(row["dosage"]), 255, row["dosage"]).astype(np.uint8).tobytes()
        payload = prefix + encoded_dosage + b"\n"
        content_digest.update(payload)
        content_bytes += len(payload)
    return variants, ambiguous_keys, {
        "url": url,
        "region": f"{normalize_chromosome(chrom)}:{start}-{end}",
        "n_biallelic_records": len(variants),
        "n_exact_duplicate_keys": exact_duplicates,
        "n_discordant_duplicate_keys_excluded": len(ambiguous_keys),
        "discordant_duplicate_key_audit": [
            {"position": key[0], "alleles": sorted(key[1]), "code": "E_REFERENCE_KEY_AMBIGUOUS_EXCLUDED", "decision": "key_level_fail_closed_exclusion_no_record_selected", "records": ambiguous_details[key]}
            for key in sorted(ambiguous_keys, key=lambda item: (item[0], sorted(item[1])))
        ],
        "n_interval_overlapping_records_excluded_by_pos": overlapping_records_excluded_by_pos,
        "n_eur_samples": len(samples),
        "transport": "pysam_htslib_remote_range_fetch_with_certifi_CA",
        "index_path": str(verified_local_index),
        "index_sha256": sha256_file(verified_local_index),
        "decoded_content_sha256": content_digest.hexdigest(),
        "decoded_canonical_byte_count": content_bytes,
        "fingerprint_scope": "sorted_position_ref_alt_plus_EUR503_uint8_dosage;255_is_missing;not_remote_BGZF_file_checksum",
    }


def correlation_matrices(raw_columns: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, int]:
    n_variants = len(raw_columns)
    if n_variants == 0:
        empty = np.empty((0, 0), dtype=np.float64)
        return empty, empty.copy(), 0
    raw = np.column_stack(raw_columns).astype(np.float64, copy=False)
    completed = raw.copy()
    for index in range(n_variants):
        missing = np.isnan(completed[:, index])
        completed[missing, index] = np.nanmean(completed[:, index])
    if n_variants == 1:
        primary = np.ones((1, 1), dtype=np.float64)
    else:
        primary = np.corrcoef(completed, rowvar=False)
    complete = raw[~np.isnan(raw).any(axis=1)]
    if complete.shape[0] < 2:
        sensitivity = np.full((n_variants, n_variants), np.nan, dtype=np.float64)
    elif n_variants == 1:
        sensitivity = np.ones((1, 1), dtype=np.float64)
    else:
        sensitivity = np.corrcoef(complete, rowvar=False)
    return primary, sensitivity, int(complete.shape[0])


def source_record(asset: Mapping[str, Any], label: str) -> dict[str, Any]:
    return {
        "label": label,
        "url": asset["stable_source_url"],
        "version_or_build": asset["genome_build"],
        "license_or_access": "public_aggregate_or_reference_data_original_provider_terms_apply",
        "retrieval_date": "2026-09-01",
        "sha256": asset.get("sha256") or "embedded_member_md5_verified_see_parent_archive",
        "byte_count": int(asset["byte_size"]),
        "structural_check": asset["integrity_status"],
    }


def hardlink_or_copy(source: Path, target: Path) -> str:
    if target.exists():
        target.unlink()
    try:
        os.link(source, target)
        return "hardlink"
    except OSError:
        shutil.copy2(source, target)
        return "copy"


def validate_ld_cache(
    primary_path: Path, complete_path: Path, variant_ids: list[str],
    primary: np.ndarray, complete: np.ndarray, digest: str,
) -> None:
    if not complete_path.exists():
        raise RuntimeError(f"E_LD_COMPLETE_CASE_CACHE_MISSING:{digest}")
    cached = np.load(primary_path, allow_pickle=False)
    cached_complete = np.load(complete_path, allow_pickle=False)
    if (
        list(cached["variant_ids"].astype(str)) != variant_ids
        or not np.array_equal(cached["ld"], primary, equal_nan=True)
        or list(cached_complete["variant_ids"].astype(str)) != variant_ids
        or not np.array_equal(cached_complete["ld"], complete, equal_nan=True)
    ):
        raise RuntimeError(f"E_LD_CACHE_HASH_COLLISION:{digest}")


def build_gwas_map(rows: list[dict[str, str]], outcome: str, region_id: str) -> tuple[dict[Any, Mapping[str, str]], dict[str, int]]:
    registry: dict[Any, Mapping[str, str]] = {}
    exact = invalid = 0
    for row in rows:
        if outcome == "CAD":
            fields = ("bp", "a1", "a2", "beta", "se", "N", "af")
            pos_field, a1_field, a2_field = "bp", "a1", "a2"
        else:
            fields = ("pos_b37", "A1", "A2", "A1_beta", "se", "N_case", "N_total", "A1_freq")
            pos_field, a1_field, a2_field = "pos_b37", "A1", "A2"
        if any(not finite(row.get(field)) for field in fields if field not in {a1_field, a2_field}):
            invalid += 1
            continue
        alleles = parse_alleles(row[a1_field], row[a2_field])
        if alleles is None:
            invalid += 1
            continue
        key = (int(row[pos_field]), alleles)
        decision = register_unique(registry, key, row, label=f"{outcome}:{region_id}")
        exact += decision == "exact_duplicate"
    return registry, {"source_rows": len(rows), "valid_unique_keys": len(registry), "exact_duplicates_collapsed": exact, "invalid_rows_excluded": invalid}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--limit-regions", type=int)
    parser.add_argument("--expected-protocol-id", default="cv_external_validation_20260901_v2")
    args = parser.parse_args()
    require_runtime()

    root, out = args.root.resolve(), args.out.resolve()
    if out.exists() and any(out.iterdir()):
        raise RuntimeError(f"E_OUTPUT_DIRECTORY_NOT_EMPTY:{out}")
    out.mkdir(parents=True, exist_ok=True)
    protocol_path = root / "config/protocol.v2.json"
    protocol = validate_protocol_freeze(protocol_path, root / "config/protocol.v2.freeze.json")
    if protocol.get("protocol_id") != args.expected_protocol_id or not protocol.get("frozen"):
        raise RuntimeError("E_PROTOCOL_NOT_FROZEN_V2")
    metadata_lock_path = root / "config/metadata.locked.v2.json"
    metadata_lock = json.loads(metadata_lock_path.read_text(encoding="utf-8"))
    locked_case_fractions = validate_metadata_lock(metadata_lock, protocol)
    manifest = json.loads((root / "audit/data_manifest.json").read_text(encoding="utf-8"))
    assets = {row["asset_id"]: row for row in manifest["assets"]}
    required_assets = set(TISSUES) | {"GCST005195", "HERMES2_EUR_outer", "HERMES2_EUR_Pheno1"}
    if not required_assets.issubset(assets):
        raise RuntimeError("E_SOURCE_MANIFEST_ASSETS_MISSING:" + ",".join(sorted(required_assets - set(assets))))

    eqtl_paths = {tissue: root / f"data/eqtl/{tissue}.all.tsv.gz" for tissue in TISSUES}
    cad_path = root / "data/cad/GCST005195_CAD_UKBIOBANK.gz"
    hf_path = root / "data/hf/primary/FORMAT-METAL_Pheno1_EUR.tsv.gz"
    panel_path = root / "data/ld/integrated_call_samples_v3.20130502.ALL.panel"
    regions_all = read_tsv(args.regions)
    source_verification = verify_manifest_assets(root, assets, regions_all)
    write_tsv(out / "source_verification.tsv", source_verification, list(source_verification[0]))
    observed_hashes = {row["asset_id"]: row["sha256"] for row in source_verification}
    validate_regions(
        regions_all, args.regions, args.regions.with_suffix(".audit.json"),
        args.regions.with_suffix(".freeze.json"), observed_hashes, protocol,
    )
    regions = regions_all[: args.limit_regions] if args.limit_regions else regions_all
    smoke = args.limit_regions is not None
    lifter = StrictRoundTrip(
        root / "data/liftover/hg19ToHg38.over.chain.gz",
        root / "data/liftover/hg38ToHg19.over.chain.gz",
    )
    eqtl_readers = {key: TabixReader(path) for key, path in eqtl_paths.items()}
    hf_reader = TabixReader(hf_path)
    samples = eur_samples(panel_path)
    cad_by_region, cad_scan_audit = scan_cad_once(
        cad_path, regions, root / "audit/cad_a1_sentinel_sign_audit.tsv", observed_hashes["GCST005195"]
    )
    write_json(out / "cad_single_scan_audit.json", cad_scan_audit)

    genes: list[dict[str, Any]] = []
    attempts: list[dict[str, Any]] = []
    region_payloads: dict[str, dict[str, Any]] = {}
    global_duplicates: list[dict[str, Any]] = []
    annotation_duplicates: list[dict[str, Any]] = []
    ld_cache = out / "ld_cache"
    ld_cache.mkdir(exist_ok=True)
    prepared_root = out / "prepared"
    prepared_root.mkdir(exist_ok=True)

    try:
        for region_index, region in enumerate(regions, start=1):
            region_id = region["region_id"]
            chrom = normalize_chromosome(region["chromosome"])
            start, end = int(region["region_start"]), int(region["region_end"])
            target_chrom, target_start, target_end = lifter.window_37_to_38(chrom, start, end)
            eqtl_raw = {
                tissue: reader.fetch(target_chrom, target_start, target_end)
                for tissue, reader in eqtl_readers.items()
            }
            candidates: dict[str, set[str]] = {}
            dedup_by_tissue: dict[str, dict[Any, Mapping[str, str]]] = {}
            candidate_counts: dict[str, dict[str, int]] = {}
            for tissue, rows in eqtl_raw.items():
                registry: dict[Any, Mapping[str, str]] = {}
                counts: dict[str, int] = defaultdict(int)
                exact = 0
                for row in rows:
                    if not complete_eqtl(row):
                        continue
                    key = (row["gene_id"], normalize_chromosome(row["chromosome"]), int(row["position"]), row["ref"].upper(), row["alt"].upper())
                    decision = register_eqtl_analytic(
                        registry, key, row, label=f"eQTL:{tissue}:{region_id}", annotation_audit=annotation_duplicates
                    )
                    if decision == "new":
                        counts[row["gene_id"]] += 1
                    else:
                        exact += 1
                dedup_by_tissue[tissue] = registry
                candidate_counts[tissue] = dict(counts)
                candidates[tissue] = {gene for gene, count in counts.items() if count >= 20}
                global_duplicates.append({"region_id": region_id, "source": tissue, "exact_duplicates_collapsed": exact, "discordant_duplicates": 0})
            shared = set.intersection(*(candidates[tissue] for tissue in TISSUES))
            gene = sorted(shared)[0] if shared else ""
            gene_status = "SELECTED" if gene else "INELIGIBLE"
            gene_code = "OK" if gene else "E_NO_ELIGIBLE_GENE"
            genes.append({
                "region_id": region_id,
                "gene_id": gene,
                "gene_selection_status": gene_status,
                "gene_selection_code": gene_code,
                "n_shared_complete_candidates": len(shared),
                **{f"n_candidates_{tissue}": len(candidates[tissue]) for tissue in TISSUES},
            })

            hf_rows = hf_reader.fetch(chrom, start, end)
            gwas_maps: dict[str, dict[Any, Mapping[str, str]]] = {}
            gwas_audits: dict[str, Any] = {}
            for outcome, rows in (("CAD", cad_by_region[region_id]), ("HF", hf_rows)):
                gwas_maps[outcome], gwas_audits[outcome] = build_gwas_map(rows, outcome, region_id)
                global_duplicates.append({
                    "region_id": region_id,
                    "source": outcome,
                    "exact_duplicates_collapsed": gwas_audits[outcome]["exact_duplicates_collapsed"],
                    "discordant_duplicates": 0,
                })
            reference, ambiguous_reference_keys, reference_audit = fetch_reference_region(
                chrom, start, end, samples,
                root / f"data/ld/remote_index_cache/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi",
            )
            global_duplicates.append({"region_id": region_id, "source": "1000G", "exact_duplicates_collapsed": reference_audit["n_exact_duplicate_keys"], "discordant_duplicates": len(ambiguous_reference_keys)})
            region_payloads[region_id] = {
                "region": region,
                "liftover": {
                    "source_build": "GRCh37", "source_region": [chrom, start, end],
                    "target_build": "GRCh38", "target_region": [target_chrom, target_start, target_end],
                    "endpoint_status": "PASS_UNIQUE_POSITIVE_STRAND_EXACT_ROUNDTRIP",
                },
                "gene": gene,
                "reference": reference_audit,
                "gwas": gwas_audits,
            }

            for tissue in TISSUES:
                eqtl_map: dict[Any, Mapping[str, str]] = {}
                eqtl_liftover_fail = eqtl_outside = exact = 0
                if gene:
                    for key, row in dedup_by_tissue[tissue].items():
                        if row["gene_id"] != gene:
                            continue
                        mapped = lifter.point(row["chromosome"], int(row["position"]), "38to37")
                        if mapped is None:
                            eqtl_liftover_fail += 1
                            continue
                        mapped_chrom, mapped_pos = mapped
                        if mapped_chrom != chrom or not start <= mapped_pos <= end:
                            eqtl_outside += 1
                            continue
                        allele_key = parse_alleles(row["ref"], row["alt"])
                        if allele_key is None:
                            continue
                        decision = register_unique(eqtl_map, (mapped_pos, allele_key), row, label=f"eQTL_lifted:{tissue}:{region_id}")
                        exact += decision == "exact_duplicate"
                global_duplicates.append({"region_id": region_id, "source": tissue + "_lifted", "exact_duplicates_collapsed": exact, "discordant_duplicates": 0})

                for outcome in OUTCOMES:
                    unit_id = f"{region_id}__{tissue}__{outcome}"
                    unit_dir = prepared_root / unit_id
                    unit_dir.mkdir(parents=True, exist_ok=True)
                    ambiguous_matched = set(eqtl_map) & set(gwas_maps[outcome]) & ambiguous_reference_keys if gene else set()
                    unit_pre_status = "PENDING"
                    unit_pre_code = "PENDING"
                    if not gene:
                        unit_pre_status, unit_pre_code = "INELIGIBLE", "E_NO_ELIGIBLE_GENE"
                    if unit_pre_status == "INELIGIBLE":
                        write_summary(unit_dir / "summary.csv", [])
                        variant_ids: list[str] = []
                        primary, complete_case, n_complete = correlation_matrices([])
                        missing_trace: list[dict[str, Any]] = []
                    else:
                        gwas_map = gwas_maps[outcome]
                        raw_keys = sorted(set(eqtl_map) & set(gwas_map) & set(reference), key=lambda item: (item[0], sorted(item[1])))
                        rows_out: list[dict[str, Any]] = []
                        raw_columns: list[np.ndarray] = []
                        missing_trace = []
                        variant_ids = []
                        for variant_key in raw_keys:
                            eqtl, gwas, refdata = eqtl_map[variant_key], gwas_map[variant_key], reference[variant_key]
                            dosage = np.asarray(refdata["dosage"], dtype=np.float64)
                            missing_fraction = float(np.mean(np.isnan(dosage)))
                            trace = {
                                "candidate_key": f"{chrom}:{variant_key[0]}:{refdata['ref']}:{refdata['alt']}",
                                "missing_fraction": missing_fraction,
                                "threshold": 0.02,
                            }
                            if missing_fraction > 0.02:
                                trace["decision"] = "excluded_missingness_gt_0.02"
                                missing_trace.append(trace)
                                continue
                            if np.all(np.isnan(dosage)) or float(np.nanstd(dosage)) == 0.0:
                                trace["decision"] = "excluded_monomorphic_or_all_missing"
                                missing_trace.append(trace)
                                continue
                            ref, alt = refdata["ref"], refdata["alt"]
                            eqtl_alt = eqtl["alt"].upper()
                            if eqtl_alt not in {ref, alt}:
                                raise RuntimeError(f"E_EQTL_EFFECT_ALLELE_UNRESOLVED:{unit_id}:{variant_key}")
                            if outcome == "CAD":
                                effect_allele = gwas["a1"].upper()
                                gwas_beta, gwas_se = float(gwas["beta"]), float(gwas["se"])
                                gwas_n, gwas_cases, gwas_af = int(float(gwas["N"])), 34541, float(gwas["af"])
                            else:
                                effect_allele = gwas["A1"].upper()
                                gwas_beta, gwas_se = float(gwas["A1_beta"]), float(gwas["se"])
                                gwas_n, gwas_cases = int(float(gwas["N_total"])), int(float(gwas["N_case"]))
                                gwas_af = float(gwas["A1_freq"])
                            if effect_allele not in {ref, alt}:
                                raise RuntimeError(f"E_GWAS_EFFECT_ALLELE_UNRESOLVED:{unit_id}:{variant_key}")
                            eqtl_beta = float(eqtl["beta"]) * (1.0 if eqtl_alt == alt else -1.0)
                            gwas_beta *= 1.0 if effect_allele == alt else -1.0
                            ld_af = float(np.nanmean(dosage) / 2.0)
                            ld_maf = min(ld_af, 1.0 - ld_af)
                            eqtl_maf = min(float(eqtl["maf"]), 1.0 - float(eqtl["maf"]))
                            gwas_maf = min(gwas_af, 1.0 - gwas_af)
                            if not (0 < ld_maf <= 0.5 and 0 < eqtl_maf <= 0.5 and 0 < gwas_maf <= 0.5):
                                trace["decision"] = "excluded_invalid_source_or_reference_maf"
                                missing_trace.append(trace)
                                continue
                            target_variant = f"{chrom}:{variant_key[0]}:{ref}:{alt}"
                            variant_ids.append(target_variant)
                            raw_columns.append(dosage)
                            trace["decision"] = "retained_mean_imputation_for_primary_LD"
                            missing_trace.append(trace)
                            rows_out.append({
                                "target_variant": target_variant,
                                "position": variant_key[0],
                                "eqtl_beta": eqtl_beta,
                                "eqtl_se": float(eqtl["se"]),
                                "gwas_beta": gwas_beta,
                                "gwas_se": gwas_se,
                                "maf": ld_maf,
                                "eqtl_maf": eqtl_maf,
                                "gwas_maf": gwas_maf,
                                "ld_maf": ld_maf,
                                "eqtl_n": int(float(eqtl["an"])) // 2,
                                "gwas_n": gwas_n,
                                "gwas_cases": gwas_cases,
                                "gwas_case_fraction_row_audit": gwas_cases / gwas_n,
                                "gwas_case_fraction_locked_study": locked_case_fractions[outcome],
                                "alignment_status": "aligned" if eqtl_alt == alt and effect_allele == alt else "flipped",
                                "vcf_status": "present",
                                "gene_id": gene,
                                "tissue_id": tissue,
                                "outcome_id": outcome,
                            })
                        primary, complete_case, n_complete = correlation_matrices(raw_columns)
                        write_summary(unit_dir / "summary.csv", rows_out)

                    order_digest = ordered_hash(variant_ids)
                    cache_primary = ld_cache / f"{order_digest}.primary.npz"
                    cache_complete = ld_cache / f"{order_digest}.complete_case.npz"
                    if not cache_primary.exists():
                        np.savez_compressed(cache_primary, ld=primary, variant_ids=np.asarray(variant_ids, dtype="U"))
                        np.savez_compressed(cache_complete, ld=complete_case, variant_ids=np.asarray(variant_ids, dtype="U"))
                    else:
                        validate_ld_cache(
                            cache_primary, cache_complete, variant_ids, primary, complete_case, order_digest
                        )
                    primary_mode = hardlink_or_copy(cache_primary, unit_dir / "ld.npz")
                    complete_mode = hardlink_or_copy(cache_complete, unit_dir / "ld_complete_case.npz")
                    overlap = {
                        "n_eqtl_variants_pre": len(eqtl_map),
                        "n_gwas_variants_pre": len(gwas_maps[outcome]),
                        "n_ld_variants_pre": len(reference),
                        "n_three_way_key_overlap_before_LD_missingness": len(set(eqtl_map) & set(gwas_maps[outcome]) & set(reference)),
                        "n_ordered_overlap": len(variant_ids),
                        "ordered_variant_sha256": order_digest,
                        "preparation_code": unit_pre_code if unit_pre_status == "INELIGIBLE" else ("OK" if variant_ids else "E_NO_ORDERED_OVERLAP"),
                        "reference_ambiguous_key_policy_code": "E_REFERENCE_KEY_AMBIGUOUS_EXCLUDED" if ambiguous_matched else "OK_NO_MATCHED_AMBIGUOUS_KEY",
                        "n_reference_ambiguous_keys_in_region": len(ambiguous_reference_keys),
                        "n_matched_discordant_reference_keys": len(ambiguous_matched),
                        "eqtl_strict_liftover_failures": eqtl_liftover_fail,
                        "eqtl_mapped_outside_source_window": eqtl_outside,
                    }
                    write_json(unit_dir / "overlap.json", overlap)
                    write_json(unit_dir / "ld_missingness_trace.json", {
                        "policy": "per_variant_mean_dosage_after_missingness_le_0.02",
                        "no_summary_field_imputation": True,
                        "primary_reference_samples": len(samples),
                        "complete_case_samples": n_complete,
                        "ordered_variant_sha256": order_digest,
                        "reference_ambiguous_keys_excluded": [
                            {"position": key[0], "alleles": sorted(key[1]), "code": "E_REFERENCE_KEY_AMBIGUOUS_EXCLUDED"}
                            for key in sorted(ambiguous_matched, key=lambda item: (item[0], sorted(item[1])))
                        ],
                        "variants": missing_trace,
                    })
                    provenance = [
                        source_record(assets[tissue], tissue),
                        source_record(assets["GCST005195"] if outcome == "CAD" else assets["HERMES2_EUR_outer"], outcome),
                        {
                            "label": "1000G_Phase3_EUR503_remote_region",
                            "url": reference_audit["url"], "version_or_build": "GRCh37_phase3_v5b_EUR503",
                            "license_or_access": "public_reference_panel_original_provider_terms_apply",
                            "retrieval_date": "2026-09-01", "sha256": reference_audit["decoded_content_sha256"],
                            "byte_count": reference_audit["decoded_canonical_byte_count"],
                            "structural_check": reference_audit["transport"],
                            "checksum_scope": reference_audit["fingerprint_scope"],
                        },
                        {
                            "label": "UCSC_hg19ToHg38_and_hg38ToHg19_chains", "url": "https://hgdownload.soe.ucsc.edu/goldenPath/",
                            "version_or_build": "GRCh37_GRCh38_bidirectional", "license_or_access": "public_UCSC_resource",
                            "retrieval_date": "2026-09-01",
                            "sha256": sha256_file(root / "data/liftover/hg19ToHg38.over.chain.gz") + ":" + sha256_file(root / "data/liftover/hg38ToHg19.over.chain.gz"),
                            "byte_count": (root / "data/liftover/hg19ToHg38.over.chain.gz").stat().st_size + (root / "data/liftover/hg38ToHg19.over.chain.gz").stat().st_size,
                            "structural_check": "unique_positive_strand_exact_point_roundtrip_required",
                        },
                    ]
                    if outcome == "HF":
                        member = source_record(assets["HERMES2_EUR_Pheno1"], "HERMES2_FORMAT-METAL_Pheno1_EUR_member")
                        member["parent_archive_sha256"] = assets["HERMES2_EUR_outer"]["sha256"]
                        member["archive_member_path"] = "HERMES2_GWAS_HF_EUR.zip/FORMAT-METAL_Pheno1_EUR.tsv.gz"
                        provenance.append(member)
                    write_json(unit_dir / "provenance.json", provenance)
                    write_json(unit_dir / "preparation_audit.json", {
                        "unit_id": unit_id, "protocol_id": protocol["protocol_id"], "region_id": region_id,
                        "region_registry_sha256": sha256_file(args.regions), "gene_id": gene,
                        "ordered_variant_sha256": order_digest, "ld_primary_materialization": primary_mode,
                        "ld_complete_case_materialization": complete_mode, "summary_sha256": sha256_file(unit_dir / "summary.csv"),
                        "variant_identity_policy": "same_GRCh37_position_and_exact_case-normalized_allele_strings_in_eQTL_GWAS_1000G_after_unique_positive_strand_exact_point_roundtrip;indels_allowed_only_under_same_exact_rule;no_strand_guess_or_normalization",
                    })
                    attempts.append({
                        "unit_id": unit_id, "region_id": region_id, "anchor_outcome": region["anchor_outcome"],
                        "tissue_id": tissue, "outcome_id": outcome, "gene_id": gene,
                        "pre_gate_status": unit_pre_status, "pre_gate_code": unit_pre_code, "n_overlap": len(variant_ids),
                    })
            write_json(out / "regions" / f"{region_id}.json", region_payloads[region_id])
            print(f"[{region_index}/{len(regions)}] {region_id}: gene={gene or 'NONE'} units=8", flush=True)
    finally:
        for reader in eqtl_readers.values():
            reader.close()
        hf_reader.close()

    expected_attempts = len(regions) * len(TISSUES) * len(OUTCOMES)
    if len(attempts) != expected_attempts:
        raise RuntimeError(f"E_ATTEMPT_COUNT_MISMATCH:{len(attempts)}:{expected_attempts}")
    if not smoke and expected_attempts != 192:
        raise RuntimeError(f"E_FULL_FACTORIAL_COUNT:{expected_attempts}")
    expected_unit_ids = {
        f"{region['region_id']}__{tissue}__{outcome}"
        for region in regions for tissue in TISSUES for outcome in OUTCOMES
    }
    observed_unit_ids = {row["unit_id"] for row in attempts}
    prepared_unit_ids = {path.name for path in prepared_root.iterdir() if path.is_dir()}
    if len(genes) != len(regions) or observed_unit_ids != expected_unit_ids or prepared_unit_ids != expected_unit_ids:
        raise RuntimeError("E_FINAL_CARTESIAN_REGISTRY_MISMATCH")
    required_unit_files = {
        "summary.csv", "ld.npz", "ld_complete_case.npz", "overlap.json",
        "provenance.json", "ld_missingness_trace.json", "preparation_audit.json",
    }
    for unit_id in sorted(expected_unit_ids):
        missing = [name for name in required_unit_files if not (prepared_root / unit_id / name).is_file()]
        if missing:
            raise RuntimeError(f"E_FINAL_UNIT_FILES_MISSING:{unit_id}:{','.join(missing)}")
    write_tsv(out / "genes.tsv", genes, list(genes[0]))
    write_tsv(out / "attempts.tsv", attempts, list(attempts[0]))
    write_tsv(out / "duplicate_key_audit.tsv", global_duplicates, list(global_duplicates[0]))
    write_tsv(
        out / "annotation_alias_audit.tsv", annotation_duplicates,
        ["source", "analytic_key", "annotation_field", "annotation_values", "decision"],
    )
    write_json_atomic(out / "run_summary.json", {
        "status": "PASS", "mode": "smoke" if smoke else "full", "protocol_id": protocol["protocol_id"],
        "regions": len(regions), "genes_selected": sum(row["gene_selection_status"] == "SELECTED" for row in genes),
        "attempts": len(attempts), "prepared_units": len(list(prepared_root.iterdir())),
        "cad_scan_count": cad_scan_audit["scan_count"],
        "reference_ambiguous_keys_quarantined": sum(int(row["discordant_duplicates"]) for row in global_duplicates if row["source"] == "1000G"),
        "discordant_duplicate_keys_remaining_in_analysis_view": 0,
        "generated_at_utc": utc_now(), "regions_sha256": sha256_file(args.regions),
        "protocol_sha256": sha256_file(protocol_path),
        "metadata_lock_sha256": sha256_file(metadata_lock_path),
        "source_verification_sha256": sha256_file(out / "source_verification.tsv"),
        "preparation_script_sha256": sha256_file(Path(__file__).resolve()),
        "runtime_indices_freeze_sha256": sha256_file(root / "config/runtime_indices.freeze.json"),
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pysam": pysam.__version__,
            "certifi_CA_sha256": sha256_file(Path(certifi.where())),
        },
        "variant_identity_policy": "exact_GRCh37_position_plus_exact_allele_strings_across_all_three_sources_after_strict_roundtrip;applies_equally_to_SNPs_and_indels",
    })


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(str(error), file=sys.stderr, flush=True)
        raise
