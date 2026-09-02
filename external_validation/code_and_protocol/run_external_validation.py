#!/usr/bin/env python3
"""Frozen, fail-closed cardiovascular eQTL--GWAS external validation runner.

This orchestrator deliberately separates prospective selection, technical
gating and downstream modelling.  It never repairs an input and never invokes
the model for an INELIGIBLE unit.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
_CODE_CANDIDATES = [
    Path(os.environ["HCSMR_CODE_DIR"]) if os.environ.get("HCSMR_CODE_DIR") else None,
    SCRIPT_DIR.parent / "hcsmr_code",
    SCRIPT_DIR.parent,
]
HCSMR_CODE = next(
    (candidate for candidate in _CODE_CANDIDATES if candidate is not None and (candidate / "hcsmr").is_dir()),
    SCRIPT_DIR.parent / "hcsmr_code",
)
if str(HCSMR_CODE) not in sys.path:
    sys.path.insert(0, str(HCSMR_CODE))

try:
    from hcsmr.reanalysis.preflight_gate import run_preflight  # type: ignore[attr-defined]  # noqa: E402
    GATE_IMPLEMENTATION = "hcsmr.reanalysis.preflight_gate"
except ModuleNotFoundError:
    # The historical server snapshot predates this module.  The vendored file is
    # an exact, tested snapshot so deployment remains immutable and auditable.
    from external_gate import run_preflight  # noqa: E402
    GATE_IMPLEMENTATION = "external_gate_vendored_snapshot"


def canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def protocol_hash(protocol: dict[str, Any]) -> str:
    return sha256_bytes(canonical_json(protocol).encode("utf-8"))


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def natural_chromosome(value: Any) -> tuple[int, str]:
    text = str(value).lower().removeprefix("chr")
    mapping = {"x": 23, "y": 24, "xy": 25, "m": 26, "mt": 26}
    try:
        return int(text), text
    except ValueError:
        return mapping.get(text, 10_000), text


def validate_protocol(protocol: dict[str, Any]) -> None:
    errors: list[str] = []
    tissue_ids = [row.get("dataset_id") for row in protocol.get("tissues", [])]
    outcome_ids = [row.get("outcome_id") for row in protocol.get("outcomes", [])]
    selection = protocol.get("region_selection", {})
    gene = protocol.get("gene_selection", {})
    factorial = protocol.get("factorial_design", {})
    if not protocol.get("frozen"):
        errors.append("protocol must be frozen")
    if tissue_ids != ["QTD000131", "QTD000136", "QTD000251", "QTD000256"]:
        errors.append("the four tissue IDs or their prospective order changed")
    if outcome_ids != ["CAD", "HF"]:
        errors.append("outcomes must be CAD then HF")
    if [row.get("dataset_id") for row in protocol.get("outcomes", [])] != ["GCST005195", "HERMES2_EUR_2025"]:
        errors.append("primary external outcomes must remain GCST005195 and HERMES2 EUR")
    protocol_id = protocol.get("protocol_id")
    cad_ancestry = protocol.get("outcomes", [{}])[0].get("ancestry")
    if protocol_id == "cv_external_validation_20260901_v1" and cad_ancestry != "UNRESOLVED_PRE_RUN_STOP":
        errors.append("v1 GCST005195 ancestry must remain unresolved")
    if protocol_id == "cv_external_validation_20260901_v2" and cad_ancestry != "European_UKB_subset":
        errors.append("v2 GCST005195 ancestry must remain the locked European UKB subset")
    if protocol_id not in {"cv_external_validation_20260901_v1", "cv_external_validation_20260901_v2"}:
        errors.append("unknown protocol amendment")
    if selection.get("selection_basis") != "GWAS_chromosome_position_and_p_value_only":
        errors.append("region selection basis changed")
    if int(selection.get("n_leads_per_outcome", -1)) != 12:
        errors.append("exactly 12 lead loci per discovery outcome are required")
    if int(selection.get("window_bp_each_side", -1)) != 500_000:
        errors.append("region window must remain +/-500 kb")
    if int(selection.get("greedy_min_distance_bp", -1)) != 1_000_000:
        errors.append("lead-locus spacing must remain 1 Mb")
    if gene.get("choice_rule") != "lexicographically_first_gene_id_among_shared_complete_candidates":
        errors.append("gene choice rule changed")
    if set(gene.get("forbidden_value_use", [])) < {"beta", "p_value", "LD", "gate_status", "posterior"}:
        errors.append("gene selection does not explicitly forbid result-based fields")
    expected = len(tissue_ids) * len(outcome_ids) * int(selection.get("expected_disease_regions", -1))
    if expected != 192 or int(factorial.get("expected_attempts", -1)) != expected:
        errors.append("factorial design must contain 192 attempts")
    if not protocol.get("gate", {}).get("no_posterior_for_non_READY"):
        errors.append("non-READY units must be forbidden from posterior generation")
    if protocol.get("gate", {}).get("ld_rank_policy") != "diagnostic":
        errors.append("LD rank must remain diagnostic for the 503-sample reference")
    if protocol.get("downstream", {}).get("ld_model_policy", {}).get("conditioning_or_jitter") != "forbidden":
        errors.append("unrecorded LD conditioning or jitter must remain forbidden")
    externality = protocol.get("externality_statement", {})
    if externality.get("permitted_label") != "prospectively_frozen_multi_tissue_external_technical_validation":
        errors.append("externality label must not imply fully independent replication")
    if errors:
        raise ValueError("invalid frozen protocol: " + "; ".join(errors))


def validate_frozen_protocol_path(path: Path) -> str:
    protocol = read_json(path)
    validate_protocol(protocol)
    freeze_path = path.with_name(path.stem + ".freeze.json")
    if not freeze_path.is_file():
        raise FileNotFoundError("frozen protocol requires protocol.freeze.json")
    freeze = read_json(freeze_path)
    observed_file_hash = sha256_file(path)
    observed_canonical_hash = protocol_hash(protocol)
    observed_bytes = path.stat().st_size
    checks = {
        "protocol_id": freeze.get("protocol_id") == protocol.get("protocol_id"),
        "frozen_at_utc": freeze.get("frozen_at_utc") == protocol.get("frozen_at_utc"),
        "protocol_byte_count": int(freeze.get("protocol_byte_count", -1)) == observed_bytes,
        "protocol_file_sha256": freeze.get("protocol_file_sha256") == observed_file_hash,
        "protocol_canonical_json_sha256": freeze.get("protocol_canonical_json_sha256") == observed_canonical_hash,
        "state": freeze.get("state") == "READ_ONLY_FROZEN",
    }
    failed = [key for key, passed in checks.items() if not passed]
    if failed:
        raise ValueError("frozen protocol manifest mismatch: " + ",".join(failed))
    return observed_canonical_hash


def validate_metadata_lock(path: Path, protocol: dict[str, Any]) -> dict[str, Any]:
    lock = read_json(path)
    if lock.get("state") != "LOCKED":
        raise ValueError("E_METADATA_NOT_LOCKED: metadata lock state is not LOCKED")
    if lock.get("protocol_canonical_json_sha256") != protocol_hash(protocol):
        raise ValueError("E_METADATA_PROTOCOL_MISMATCH")
    expected = {row["outcome_id"]: row["dataset_id"] for row in protocol["outcomes"]}
    for outcome, dataset_id in expected.items():
        record = lock.get("outcomes", {}).get(outcome, {})
        required = ["ancestry", "genome_build", "cases", "controls", "sample_size", "case_fraction", "metadata_url", "source_citation"]
        missing = [field for field in required if record.get(field) in (None, "", "UNRESOLVED")]
        if record.get("dataset_id") != dataset_id or missing:
            raise ValueError(f"E_METADATA_NOT_LOCKED: {outcome} unresolved fields {missing}")
        cases, controls, sample_size = int(record["cases"]), int(record["controls"]), int(record["sample_size"])
        fraction = float(record["case_fraction"])
        if cases <= 0 or controls <= 0 or sample_size <= 0 or cases + controls != sample_size:
            raise ValueError(f"E_METADATA_INVALID: {outcome} sample counts")
        if abs(fraction - cases / sample_size) > 1e-12:
            raise ValueError(f"E_METADATA_INVALID: {outcome} case fraction")
    return lock


def read_table(path: Path, separator: str | None = None, usecols: list[str] | None = None) -> pd.DataFrame:
    if separator is None:
        separator = "\t" if not path.name.endswith(".csv") and ".tsv" in path.name else ","
    if separator in {r"\\s+", "whitespace"}:
        separator = r"\s+"
    if separator in {r"\t", "tab"}:
        separator = "\t"
    return pd.read_csv(path, sep=separator, compression="infer", low_memory=False, usecols=usecols)


def select_leads(
    frame: pd.DataFrame,
    *,
    outcome: str,
    chromosome_col: str,
    position_col: str,
    p_col: str,
    variant_col: str,
    n_leads: int,
    p_threshold: float,
    min_distance_bp: int,
    window_bp: int,
) -> list[dict[str, Any]]:
    """Select leads using GWAS coordinates and P values only.

    No eQTL, overlap, LD, gate or downstream quantity is accepted by this
    function.  Greedy distance pruning follows a fully deterministic sort.
    """
    required = [chromosome_col, position_col, p_col, variant_col]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{outcome} GWAS is missing selection columns: {missing}")
    work = frame[required].copy()
    work.columns = ["chromosome", "position", "p_value", "variant_id"]
    work["position"] = pd.to_numeric(work["position"], errors="coerce")
    work["p_value"] = pd.to_numeric(work["p_value"], errors="coerce")
    work["variant_id"] = work["variant_id"].astype(str)
    work = work.loc[
        work["position"].notna()
        & work["p_value"].notna()
        & (work["position"] > 0)
        & (work["p_value"] > 0)
        & (work["p_value"] <= p_threshold)
    ].copy()
    work["position"] = work["position"].astype(int)
    work["chrom_key"] = work["chromosome"].map(natural_chromosome)
    work = work.sort_values(
        ["p_value", "chrom_key", "position", "variant_id"], kind="mergesort"
    )
    chosen: list[dict[str, Any]] = []
    positions: dict[str, list[int]] = defaultdict(list)
    for row in work.itertuples(index=False):
        chrom = str(row.chromosome).removeprefix("chr")
        pos = int(row.position)
        if any(abs(pos - existing) < min_distance_bp for existing in positions[chrom]):
            continue
        positions[chrom].append(pos)
        chosen.append({
            "region_id": f"{outcome}_L{len(chosen) + 1:02d}",
            "anchor_outcome": outcome,
            "lead_variant": str(row.variant_id),
            "chromosome": chrom,
            "lead_position": pos,
            "lead_p_value": float(row.p_value),
            "region_start": max(1, pos - window_bp),
            "region_end": pos + window_bp,
            "selection_rank": len(chosen) + 1,
            "selection_basis": "GWAS_chromosome_position_and_p_value_only",
        })
        if len(chosen) == n_leads:
            break
    if len(chosen) != n_leads:
        raise RuntimeError(
            f"E_LEAD_LOCUS_SHORTFALL: {outcome} yielded {len(chosen)} of {n_leads} required loci"
        )
    return chosen


def cmd_select_loci(args: argparse.Namespace) -> None:
    validate_frozen_protocol_path(args.protocol)
    protocol = read_json(args.protocol)
    validate_metadata_lock(args.metadata_lock, protocol)
    cfg = protocol["region_selection"]
    rows: list[dict[str, Any]] = []
    for outcome, path, columns, separator in (
        ("CAD", args.cad, args.cad_columns, args.cad_separator),
        ("HF", args.hf, args.hf_columns, args.hf_separator),
    ):
        mapping = json.loads(columns)
        frame = read_table(path, separator, usecols=list(mapping.values()))
        rows.extend(select_leads(
            frame,
            outcome=outcome,
            chromosome_col=mapping["chromosome"],
            position_col=mapping["position"],
            p_col=mapping["p_value"],
            variant_col=mapping["variant_id"],
            n_leads=int(cfg["n_leads_per_outcome"]),
            p_threshold=float(cfg["p_value_threshold"]),
            min_distance_bp=int(cfg["greedy_min_distance_bp"]),
            window_bp=int(cfg["window_bp_each_side"]),
        ))
    if len(rows) != int(cfg["expected_disease_regions"]):
        raise RuntimeError("E_DISEASE_REGION_COUNT_MISMATCH")
    write_tsv(args.out, rows, list(rows[0]))
    write_json(args.out.with_suffix(".audit.json"), {
        "protocol_id": protocol["protocol_id"],
        "protocol_sha256": protocol_hash(protocol),
        "generated_at_utc": utc_now(),
        "selection_inputs": {
            "CAD": {"sha256": sha256_file(args.cad), "columns": json.loads(args.cad_columns), "separator": args.cad_separator},
            "HF": {"sha256": sha256_file(args.hf), "columns": json.loads(args.hf_columns), "separator": args.hf_separator},
        },
        "n_regions": len(rows),
        "selection_basis": cfg["selection_basis"],
        "downstream_results_consulted": False,
    })


def complete_gene_candidates(
    frame: pd.DataFrame, region: pd.Series, required: list[str], minimum: int
) -> set[str]:
    missing = [field for field in required if field not in frame.columns]
    if missing:
        raise ValueError(f"eQTL candidate table missing fields: {missing}")
    chromosome = frame["chromosome"].astype(str).str.removeprefix("chr")
    position = pd.to_numeric(frame["position"], errors="coerce")
    subset = frame.loc[
        (chromosome == str(region["chromosome"]).removeprefix("chr"))
        & (position >= int(region["region_start"]))
        & (position <= int(region["region_end"]))
    ].copy()
    # Completeness only; magnitudes/significance are not examined.
    complete = subset.loc[subset[required].notna().all(axis=1)]
    counts = complete.groupby(complete["gene_id"].astype(str), sort=True).size()
    return set(counts.loc[counts >= minimum].index)


def cmd_select_genes(args: argparse.Namespace) -> None:
    validate_frozen_protocol_path(args.protocol)
    protocol = read_json(args.protocol)
    regions = pd.read_csv(args.regions, sep="\t", dtype={"chromosome": str})
    tissue_paths = json.loads(args.tissue_tables)
    expected_tissues = [row["dataset_id"] for row in protocol["tissues"]]
    if sorted(tissue_paths) != sorted(expected_tissues):
        raise ValueError("exactly the four frozen tissue tables are required")
    tissue_frames = {key: read_table(Path(value), args.separator) for key, value in tissue_paths.items()}
    gene_cfg = protocol["gene_selection"]
    required = list(gene_cfg["required_nonmissing_fields"])
    minimum = int(gene_cfg["min_complete_rows_per_tissue"])
    output: list[dict[str, Any]] = []
    for _, region in regions.iterrows():
        candidates_by_tissue = {
            tissue: complete_gene_candidates(frame, region, required, minimum)
            for tissue, frame in tissue_frames.items()
        }
        shared = set.intersection(*candidates_by_tissue.values()) if candidates_by_tissue else set()
        chosen = sorted(shared)[0] if shared else ""
        output.append({
            "region_id": region["region_id"],
            "gene_id": chosen,
            "gene_selection_status": "SELECTED" if chosen else "INELIGIBLE",
            "gene_selection_code": "OK" if chosen else gene_cfg["no_candidate_code"],
            "n_shared_complete_candidates": len(shared),
            **{f"n_candidates_{tissue}": len(candidates_by_tissue[tissue]) for tissue in expected_tissues},
        })
    fields = list(output[0])
    write_tsv(args.out, output, fields)
    write_json(args.out.with_suffix(".audit.json"), {
        "protocol_sha256": protocol_hash(protocol),
        "region_registry_sha256": sha256_file(args.regions),
        "tissue_input_sha256": {key: sha256_file(Path(value)) for key, value in tissue_paths.items()},
        "rule": gene_cfg["choice_rule"],
        "forbidden_values_used": [],
        "generated_at_utc": utc_now(),
    })


def cmd_plan(args: argparse.Namespace) -> None:
    validate_frozen_protocol_path(args.protocol)
    protocol = read_json(args.protocol)
    regions = pd.read_csv(args.regions, sep="\t", dtype=str)
    genes = pd.read_csv(args.genes, sep="\t", dtype=str).set_index("region_id")
    attempts: list[dict[str, Any]] = []
    for region in regions.to_dict(orient="records"):
        gene = genes.loc[region["region_id"]]
        for tissue in protocol["tissues"]:
            for outcome in protocol["outcomes"]:
                unit_id = f"{region['region_id']}__{tissue['dataset_id']}__{outcome['outcome_id']}"
                eligible = gene["gene_selection_status"] == "SELECTED"
                attempts.append({
                    "unit_id": unit_id,
                    "region_id": region["region_id"],
                    "anchor_outcome": region["anchor_outcome"],
                    "tissue_id": tissue["dataset_id"],
                    "outcome_id": outcome["outcome_id"],
                    "gene_id": gene["gene_id"] if eligible else "",
                    "pre_gate_status": "PENDING" if eligible else "INELIGIBLE",
                    "pre_gate_code": "PENDING" if eligible else gene["gene_selection_code"],
                })
    if len(attempts) != int(protocol["factorial_design"]["expected_attempts"]):
        raise RuntimeError("E_ATTEMPT_COUNT_MISMATCH")
    write_tsv(args.out, attempts, list(attempts[0]))
    write_json(args.out.with_suffix(".audit.json"), {
        "protocol_sha256": protocol_hash(protocol),
        "region_registry_sha256": sha256_file(args.regions),
        "gene_registry_sha256": sha256_file(args.genes),
        "n_attempts": len(attempts),
        "generated_at_utc": utc_now(),
    })


def load_gate_inputs(unit_dir: Path) -> tuple[pd.DataFrame, np.ndarray, list[str], list[dict[str, Any]], dict[str, Any]]:
    summary_path = unit_dir / "summary.csv"
    ld_path = unit_dir / "ld.npz"
    provenance_path = unit_dir / "provenance.json"
    overlap_path = unit_dir / "overlap.json"
    missing = [str(path.name) for path in (summary_path, ld_path, provenance_path, overlap_path) if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing prepared unit inputs: " + ",".join(missing))
    summary = pd.read_csv(summary_path, dtype=str).where(lambda x: x.notna(), None)
    archive = np.load(ld_path, allow_pickle=False)
    if "ld" not in archive.files or "variant_ids" not in archive.files:
        raise ValueError("ld.npz must contain ld and variant_ids")
    variant_ids = [str(value) for value in archive["variant_ids"]]
    overlap = read_json(overlap_path)
    required_overlap = {
        "n_eqtl_variants_pre", "n_gwas_variants_pre", "n_ld_variants_pre",
        "n_ordered_overlap", "ordered_variant_sha256",
    }
    if required_overlap.difference(overlap):
        raise ValueError("E_OVERLAP_AUDIT_INVALID: required overlap fields are missing")
    if int(overlap["n_ordered_overlap"]) != len(summary) or len(summary) != len(variant_ids):
        raise ValueError("E_OVERLAP_AUDIT_INVALID: overlap, summary and LD identifier counts differ")
    observed_order_hash = sha256_bytes(("\n".join(variant_ids) + "\n").encode("utf-8"))
    if overlap["ordered_variant_sha256"] != observed_order_hash:
        raise ValueError("E_OVERLAP_AUDIT_INVALID: ordered variant hash differs from LD archive")
    return summary, archive["ld"], variant_ids, read_json(provenance_path), overlap


def case_fraction_for(outcome: str, case_fractions: dict[str, Any]) -> float:
    if outcome not in case_fractions:
        raise ValueError(f"case fraction was not prospectively supplied for {outcome}")
    value = float(case_fractions[outcome])
    if not 0 < value < 1:
        raise ValueError(f"invalid case fraction for {outcome}")
    return value


def cmd_gate(args: argparse.Namespace) -> None:
    validate_frozen_protocol_path(args.protocol)
    protocol = read_json(args.protocol)
    metadata_lock = validate_metadata_lock(args.metadata_lock, protocol)
    attempts = pd.read_csv(args.attempts, sep="\t", dtype=str)
    case_fractions = {key: value["case_fraction"] for key, value in metadata_lock["outcomes"].items()}
    output_root = args.out
    trace: list[dict[str, Any]] = []
    for attempt in attempts.to_dict(orient="records"):
        unit_id = attempt["unit_id"]
        unit_out = output_root / "units" / unit_id
        gate_path = unit_out / "gate.json"
        if attempt["pre_gate_status"] == "INELIGIBLE":
            result = {
                "schema_version": "external_validation_1.0",
                "status": "INELIGIBLE",
                "status_codes": [attempt["pre_gate_code"]],
                "failures": [{"code": attempt["pre_gate_code"], "message": "Prospective pre-gate preparation stop"}],
                "gate_ready_for_dispatch": False,
                "decision_scope": "declared_gate_policy_conformance_only",
            }
            write_json(gate_path, result)
        else:
            prepared = args.prepared / unit_id
            try:
                summary, ld, variant_ids, provenance, overlap = load_gate_inputs(prepared)
                if int(overlap["n_ordered_overlap"]) == 0:
                    result = {
                        "schema_version": "external_validation_1.0",
                        "status": "INELIGIBLE",
                        "status_codes": ["E_NO_ORDERED_OVERLAP"],
                        "failures": [{"code": "E_NO_ORDERED_OVERLAP", "message": "No exactly ordered eQTL-GWAS-LD variant overlap"}],
                        "n_summary_rows": 0,
                        "n_ld_variant_ids": 0,
                        "ld_diagnostics": {"n_variants": 0, "numerical_rank": 0, "condition_number_2": None},
                    }
                else:
                    result = run_preflight(
                        rows=summary.to_dict(orient="records"),
                        ld=ld,
                        ld_variant_ids=variant_ids,
                        provenance=provenance,
                        eqtl_trait_type="quant",
                        gwas_trait_type="cc",
                        case_fraction=case_fraction_for(attempt["outcome_id"], case_fractions),
                        case_fraction_policy=protocol["gate"]["case_fraction_policy"],
                        ld_rank_policy=protocol["gate"]["ld_rank_policy"],
                        max_condition_number=float(protocol["gate"]["max_condition_number"]),
                    )
                result["gate_ready_for_dispatch"] = result["status"] == "READY"
                result["overlap_diagnostics"] = overlap
                result["protocol_sha256"] = protocol_hash(protocol)
                result["prepared_input_sha256"] = {
                    name: sha256_file(prepared / name) for name in ("summary.csv", "ld.npz", "provenance.json", "overlap.json")
                }
                write_json(gate_path, result)
            except Exception as error:  # fail closed; preserve audit instead of aborting other units
                message = str(error)
                error_code = message.split(":", 1)[0] if message.startswith("E_") else "E_INPUT_PREPARATION_OR_GATE_EXCEPTION"
                result = {
                    "schema_version": "external_validation_1.0",
                    "status": "INELIGIBLE",
                    "status_codes": [error_code],
                    "failures": [{"code": error_code, "message": message}],
                    "gate_ready_for_dispatch": False,
                    "protocol_sha256": protocol_hash(protocol),
                }
                write_json(gate_path, result)
        trace.append({
            "unit_id": unit_id,
            "stage": "gate",
            "status": result["status"],
            "status_codes": ";".join(result["status_codes"]),
            "gate_sha256": sha256_file(gate_path),
        })
    write_tsv(output_root / "trace.tsv", trace, list(trace[0]))
    write_run_audit(output_root, protocol, args)


def write_run_audit(out: Path, protocol: dict[str, Any], args: argparse.Namespace) -> None:
    config_rows = [
        {"key": "protocol_id", "value": protocol["protocol_id"]},
        {"key": "protocol_sha256", "value": protocol_hash(protocol)},
        {"key": "seed", "value": protocol["seed"]},
        {"key": "generated_at_utc", "value": utc_now()},
    ]
    write_tsv(out / "run_configuration.tsv", config_rows, ["key", "value"])
    (out / "environment.txt").write_text(
        f"python={sys.version}\nplatform={platform.platform()}\nnumpy={np.__version__}\npandas={pd.__version__}\ngate_implementation={GATE_IMPLEMENTATION}\n",
        encoding="utf-8",
    )


def cmd_dispatch(args: argparse.Namespace) -> None:
    validate_frozen_protocol_path(args.protocol)
    protocol = read_json(args.protocol)
    metadata_lock = validate_metadata_lock(args.metadata_lock, protocol)
    attempts = pd.read_csv(args.attempts, sep="\t", dtype=str)
    case_fractions = {key: value["case_fraction"] for key, value in metadata_lock["outcomes"].items()}
    for attempt in attempts.to_dict(orient="records"):
        unit_id = attempt["unit_id"]
        gate_path = args.gates / "units" / unit_id / "gate.json"
        model_out = args.out / "units" / unit_id
        model_out.mkdir(parents=True, exist_ok=True)
        if not gate_path.is_file():
            write_json(model_out / "model_status.json", {"status": "NOT_RUN", "code": "E_GATE_MISSING"})
            continue
        gate = read_json(gate_path)
        posterior_files = list(model_out.glob("*posterior*")) + list(model_out.glob("*.rds"))
        if gate.get("status") != "READY":
            if posterior_files:
                raise RuntimeError(f"E_STALE_POSTERIOR_FOR_INELIGIBLE: {unit_id}")
            write_json(model_out / "model_status.json", {
                "status": "NOT_RUN",
                "code": "BLOCKED_BY_FAIL_CLOSED_GATE",
                "gate_status": gate.get("status"),
                "gate_sha256": sha256_file(gate_path),
            })
            continue
        prepared = args.prepared / unit_id
        summary, ld, variant_ids, _, _ = load_gate_inputs(prepared)
        if summary["target_variant"].astype(str).tolist() != variant_ids:
            raise RuntimeError(f"E_LD_VARIANT_ORDER_MISMATCH before model dispatch: {unit_id}")
        locked = model_out / "locked_input"
        locked.mkdir(exist_ok=True)
        summary.to_csv(locked / "summary.csv", index=False)
        pd.DataFrame(ld, index=variant_ids, columns=variant_ids).to_csv(locked / "ld.csv")
        write_json(locked / "binding.json", {
            "protocol_sha256": protocol_hash(protocol),
            "gate_sha256": sha256_file(gate_path),
            "summary_sha256": sha256_file(locked / "summary.csv"),
            "ld_sha256": sha256_file(locked / "ld.csv"),
        })
        command = [
            args.rscript,
            str(args.r_runner),
            str(locked),
            str(case_fraction_for(attempt["outcome_id"], case_fractions)),
            str(model_out),
        ]
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        (model_out / "stdout.log").write_text(completed.stdout, encoding="utf-8")
        (model_out / "stderr.log").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            write_json(model_out / "model_status.json", {"status": "FAILED", "code": "E_MODEL_RUNNER", "returncode": completed.returncode})
        else:
            write_json(model_out / "model_status.json", {
                "status": "COMPLETED",
                "code": "OK",
                "gate_sha256": sha256_file(gate_path),
                "runner_sha256": sha256_file(args.r_runner),
            })


def cmd_summarize(args: argparse.Namespace) -> None:
    attempts = pd.read_csv(args.attempts, sep="\t", dtype=str)
    rows: list[dict[str, Any]] = []
    for attempt in attempts.to_dict(orient="records"):
        gate_path = args.gates / "units" / attempt["unit_id"] / "gate.json"
        gate = read_json(gate_path) if gate_path.is_file() else {"status": "MISSING", "status_codes": ["E_GATE_MISSING"]}
        diagnostics = gate.get("ld_diagnostics", {})
        rows.append({
            **attempt,
            "gate_status": gate.get("status"),
            "error_codes": ";".join(gate.get("status_codes", [])),
            "variant_overlap_n": gate.get("overlap_diagnostics", {}).get("n_ordered_overlap", gate.get("n_summary_rows")),
            "variant_overlap_fraction_eqtl": (
                gate.get("overlap_diagnostics", {}).get("n_ordered_overlap", 0)
                / gate.get("overlap_diagnostics", {}).get("n_eqtl_variants_pre", 1)
                if gate.get("overlap_diagnostics", {}).get("n_eqtl_variants_pre", 0) else None
            ),
            "variant_overlap_fraction_gwas": (
                gate.get("overlap_diagnostics", {}).get("n_ordered_overlap", 0)
                / gate.get("overlap_diagnostics", {}).get("n_gwas_variants_pre", 1)
                if gate.get("overlap_diagnostics", {}).get("n_gwas_variants_pre", 0) else None
            ),
            "ld_dimension": diagnostics.get("n_variants"),
            "ld_rank": diagnostics.get("numerical_rank"),
            "ld_condition_number": diagnostics.get("condition_number_2"),
            "posterior_generated": any((args.models / "units" / attempt["unit_id"]).rglob("*.rds")) if args.models else False,
        })
    write_tsv(args.out, rows, list(rows[0]))
    bad = [row for row in rows if row["gate_status"] != "READY" and row["posterior_generated"]]
    status_counts = Counter(row["gate_status"] for row in rows)
    write_json(args.out.with_suffix(".json"), {
        "n_attempts": len(rows),
        "gate_status_counts": status_counts,
        "non_READY_with_posterior": len(bad),
        "integrity_status": "PASS" if not bad else "FAIL",
    })
    if bad:
        raise RuntimeError("E_POSTERIOR_GENERATED_FOR_NON_READY")


def cmd_compare_reruns(args: argparse.Namespace) -> None:
    attempts = pd.read_csv(args.attempts, sep="\t", dtype=str)
    rows: list[dict[str, Any]] = []
    for unit_id in attempts["unit_id"].tolist():
        left_path = args.left / "units" / unit_id / "gate.json"
        right_path = args.right / "units" / unit_id / "gate.json"
        if not left_path.is_file() or not right_path.is_file():
            rows.append({
                "unit_id": unit_id, "status": "FAIL", "byte_identical": False,
                "semantic_identical": False, "reason": "missing_gate_record",
            })
            continue
        left, right = read_json(left_path), read_json(right_path)
        semantic_fields = ["status", "status_codes", "n_summary_rows", "n_ld_variant_ids", "ld_diagnostics", "overlap_diagnostics"]
        semantic_identical = all(left.get(key) == right.get(key) for key in semantic_fields)
        byte_identical = sha256_file(left_path) == sha256_file(right_path)
        rows.append({
            "unit_id": unit_id,
            "status": "PASS" if byte_identical and semantic_identical else "FAIL",
            "byte_identical": byte_identical,
            "semantic_identical": semantic_identical,
            "reason": "" if byte_identical and semantic_identical else "gate_record_mismatch",
        })
    write_tsv(args.out, rows, list(rows[0]))
    failures = sum(row["status"] == "FAIL" for row in rows)
    write_json(args.out.with_suffix(".json"), {
        "n_units": len(rows), "n_failures": failures,
        "exact_rerun_status": "PASS" if failures == 0 else "FAIL",
    })
    if failures:
        raise RuntimeError(f"E_RERUN_MISMATCH: {failures} unit(s)")


def cmd_verify_sources(args: argparse.Namespace) -> None:
    sources = read_json(args.sources)
    local_map = read_json(args.local_files)
    locked = []
    for source in sources:
        label = source["label"]
        path = Path(local_map[label]) if label in local_map else None
        if path is None or not path.is_file():
            raise FileNotFoundError(f"source {label} has no local file")
        record = dict(source)
        record.update(retrieval_date=args.retrieval_date, sha256=sha256_file(path), byte_count=path.stat().st_size)
        locked.append(record)
    write_json(args.out, locked)


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    sub = root.add_subparsers(dest="command", required=True)
    p = sub.add_parser("validate-protocol")
    p.add_argument("--protocol", type=Path, required=True)
    p.set_defaults(func=lambda args: print(validate_frozen_protocol_path(args.protocol)))

    p = sub.add_parser("select-loci")
    p.add_argument("--protocol", type=Path, required=True)
    p.add_argument("--metadata-lock", type=Path, required=True)
    p.add_argument("--cad", type=Path, required=True)
    p.add_argument("--hf", type=Path, required=True)
    p.add_argument("--cad-columns", required=True, help='JSON mapping with chromosome,position,p_value,variant_id')
    p.add_argument("--hf-columns", required=True)
    p.add_argument("--cad-separator", default="\\s+")
    p.add_argument("--hf-separator", default="\t")
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_select_loci)

    p = sub.add_parser("select-genes")
    p.add_argument("--protocol", type=Path, required=True)
    p.add_argument("--regions", type=Path, required=True)
    p.add_argument("--tissue-tables", required=True, help="JSON mapping of QTD ID to regional-candidate table")
    p.add_argument("--separator", default="\t")
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_select_genes)

    p = sub.add_parser("plan")
    p.add_argument("--protocol", type=Path, required=True)
    p.add_argument("--regions", type=Path, required=True)
    p.add_argument("--genes", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_plan)

    p = sub.add_parser("gate")
    p.add_argument("--protocol", type=Path, required=True)
    p.add_argument("--metadata-lock", type=Path, required=True)
    p.add_argument("--attempts", type=Path, required=True)
    p.add_argument("--prepared", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_gate)

    p = sub.add_parser("dispatch")
    p.add_argument("--protocol", type=Path, required=True)
    p.add_argument("--metadata-lock", type=Path, required=True)
    p.add_argument("--attempts", type=Path, required=True)
    p.add_argument("--prepared", type=Path, required=True)
    p.add_argument("--gates", type=Path, required=True)
    p.add_argument("--rscript", default="Rscript")
    p.add_argument("--r-runner", type=Path, default=SCRIPT_DIR / "run_coloc_susie_rank_diagnostic.R")
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_dispatch)

    p = sub.add_parser("summarize")
    p.add_argument("--attempts", type=Path, required=True)
    p.add_argument("--gates", type=Path, required=True)
    p.add_argument("--models", type=Path)
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_summarize)

    p = sub.add_parser("compare-reruns")
    p.add_argument("--attempts", type=Path, required=True)
    p.add_argument("--left", type=Path, required=True)
    p.add_argument("--right", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_compare_reruns)

    p = sub.add_parser("verify-sources")
    p.add_argument("--sources", type=Path, required=True)
    p.add_argument("--local-files", type=Path, required=True)
    p.add_argument("--retrieval-date", required=True)
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_verify_sources)
    return root


def main() -> None:
    args = parser().parse_args()
    np.random.seed(20260901)
    args.func(args)


if __name__ == "__main__":
    main()
