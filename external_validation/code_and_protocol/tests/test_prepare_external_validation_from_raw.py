import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


SCRIPT = Path(__file__).parents[1] / "prepare_external_validation_from_raw.py"
SPEC = importlib.util.spec_from_file_location("prepare_external_validation_from_raw", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_register_unique_collapses_only_exact_duplicates():
    registry = {}
    assert MODULE.register_unique(registry, (1, "A", "G"), {"beta": "1"}, label="fixture") == "new"
    assert MODULE.register_unique(registry, (1, "A", "G"), {"beta": "1"}, label="fixture") == "exact_duplicate"
    with pytest.raises(RuntimeError, match="E_DISCORDANT_DUPLICATE_KEY"):
        MODULE.register_unique(registry, (1, "A", "G"), {"beta": "2"}, label="fixture")


def test_correlation_is_on_final_order_and_has_complete_case_sensitivity():
    first = np.array([0.0, 1.0, 2.0, np.nan])
    second = np.array([0.0, 1.0, 2.0, 1.0])
    primary, complete, n_complete = MODULE.correlation_matrices([first, second])
    assert primary.shape == (2, 2)
    assert complete.shape == (2, 2)
    assert n_complete == 3
    assert np.allclose(np.diag(primary), 1.0)
    assert np.allclose(complete, np.ones((2, 2)))


def test_empty_overlap_has_stable_zero_dimensional_ld():
    primary, complete, n_complete = MODULE.correlation_matrices([])
    assert primary.shape == complete.shape == (0, 0)
    assert n_complete == 0
    assert MODULE.ordered_hash([]) == MODULE.ordered_hash([])


def test_rsid_only_duplicate_is_analytic_exact_and_audited():
    registry, audit = {}, []
    base = {"gene_id": "ENSG1", "position": "10", "beta": "0.1", "se": "0.02", "rsid": "rs1"}
    assert MODULE.register_eqtl_analytic(registry, ("ENSG1", 10), base, label="qtl", annotation_audit=audit) == "new"
    alias = {**base, "rsid": "rs2"}
    assert MODULE.register_eqtl_analytic(registry, ("ENSG1", 10), alias, label="qtl", annotation_audit=audit) == "analytic_exact_duplicate"
    assert audit[0]["annotation_values"] == "rs1;rs2"
    with pytest.raises(RuntimeError, match="E_DISCORDANT_DUPLICATE_KEY"):
        MODULE.register_eqtl_analytic(registry, ("ENSG1", 10), {**alias, "beta": "0.2"}, label="qtl", annotation_audit=audit)


def test_protocol_freeze_detects_post_freeze_tamper(tmp_path):
    protocol_path = tmp_path / "protocol.v2.json"
    freeze_path = tmp_path / "protocol.v2.freeze.json"
    protocol = {"protocol_id": "p", "frozen_at_utc": "t", "frozen": True}
    protocol_path.write_text(json.dumps(protocol) + "\n", encoding="utf-8")
    freeze = {
        "state": "READ_ONLY_FROZEN", "protocol_id": "p", "frozen_at_utc": "t",
        "protocol_byte_count": protocol_path.stat().st_size,
        "protocol_file_sha256": MODULE.sha256_file(protocol_path),
        "protocol_canonical_json_sha256": MODULE.canonical_json_sha256(protocol),
    }
    freeze_path.write_text(json.dumps(freeze), encoding="utf-8")
    assert MODULE.validate_protocol_freeze(protocol_path, freeze_path) == protocol
    protocol_path.write_text(json.dumps({**protocol, "seed": 2}) + "\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="E_PROTOCOL_FREEZE_MISMATCH"):
        MODULE.validate_protocol_freeze(protocol_path, freeze_path)


def test_regions_design_tamper_is_blocked(tmp_path):
    rows = []
    for outcome in MODULE.OUTCOMES:
        for rank in range(1, 13):
            lead = 2_000_000 + rank * 2_000_000
            rows.append({
                "region_id": f"{outcome}_L{rank:02d}", "anchor_outcome": outcome,
                "selection_rank": str(rank), "lead_position": str(lead),
                "region_start": str(lead - 500_000), "region_end": str(lead + 500_000),
                "lead_p_value": "1e-9", "selection_basis": "GWAS_chromosome_position_and_p_value_only",
            })
    rows[0]["region_end"] = str(int(rows[0]["region_end"]) + 1)
    audit = tmp_path / "regions.audit.json"
    audit.write_text(json.dumps({"n_regions": 24, "selection_basis": "GWAS_chromosome_position_and_p_value_only", "CAD_sha256": "c", "HF_sha256": "h"}), encoding="utf-8")
    regions_path = tmp_path / "regions.tsv"
    regions_path.write_text("fixture\n", encoding="utf-8")
    freeze_path = tmp_path / "regions.freeze.json"
    protocol = {"outcomes": []}
    freeze_path.write_text(json.dumps({"state": "READ_ONLY_FROZEN", "regions_sha256": MODULE.sha256_file(regions_path), "audit_sha256": MODULE.sha256_file(audit), "protocol_canonical_json_sha256": MODULE.canonical_json_sha256(protocol), "n_regions": 24}), encoding="utf-8")
    with pytest.raises(RuntimeError, match="E_REGION_WINDOW_OR_THRESHOLD"):
        MODULE.validate_regions(rows, regions_path, audit, freeze_path, {"GCST005195": "c", "HERMES2_EUR_Pheno1": "h"}, protocol)


def test_complete_case_cache_tamper_is_blocked(tmp_path):
    variants = ["1:1:A:G", "1:2:C:T"]
    primary = np.eye(2)
    complete = np.ones((2, 2))
    p = tmp_path / "primary.npz"
    c = tmp_path / "complete.npz"
    np.savez_compressed(p, ld=primary, variant_ids=np.asarray(variants))
    np.savez_compressed(c, ld=np.eye(2), variant_ids=np.asarray(variants))
    with pytest.raises(RuntimeError, match="E_LD_CACHE_HASH_COLLISION"):
        MODULE.validate_ld_cache(p, c, variants, primary, complete, "digest")


def test_manifest_tbi_hash_tamper_is_blocked(tmp_path):
    shared = tmp_path / "shared.bin"
    shared.write_bytes(b"source")
    hf = tmp_path / "data/hf/primary/FORMAT-METAL_Pheno1_EUR.tsv.gz"
    hf.parent.mkdir(parents=True)
    hf.write_bytes(b"hf")
    (hf.with_name(hf.name + ".md5")).write_text(MODULE.md5_file(hf) + "  member\n", encoding="utf-8")
    hf_tbi = hf.with_name(hf.name + ".tbi")
    hf_tbi.write_bytes(b"index")
    config = tmp_path / "config"
    config.mkdir()
    (config / "runtime_indices.freeze.json").write_text(json.dumps({"state": "READ_ONLY_FROZEN", "indices": {"HERMES2_EUR_Pheno1_tbi": {"project_relative_path": "data/hf/primary/FORMAT-METAL_Pheno1_EUR.tsv.gz.tbi", "byte_count": hf_tbi.stat().st_size, "sha256": MODULE.sha256_file(hf_tbi)}}}), encoding="utf-8")
    ids = set(MODULE.TISSUES) | {f"{x}_tbi" for x in MODULE.TISSUES} | {
        "GCST005195", "HERMES2_EUR_outer", "UCSC_hg19ToHg38_chain",
        "UCSC_hg38ToHg19_chain", "1000G_panel", "1000G_chr1_tbi_cache",
    }
    assets = {
        asset_id: {"project_relative_path": str(shared.relative_to(tmp_path)), "byte_size": shared.stat().st_size, "sha256": MODULE.sha256_file(shared)}
        for asset_id in ids
    }
    assets["HERMES2_EUR_Pheno1"] = {"project_relative_path": str(hf.relative_to(tmp_path)), "byte_size": hf.stat().st_size, "sha256": ""}
    assets["QTD000131_tbi"]["sha256"] = "0" * 64
    with pytest.raises(RuntimeError, match="E_MANIFEST_SHA256_MISMATCH:QTD000131_tbi"):
        MODULE.verify_manifest_assets(tmp_path, assets, [{"chromosome": "1"}])


def test_reference_overlap_record_is_countable_but_other_outside_records_stop():
    assert MODULE.reference_record_position_decision("9", 90, 110, "9", 100, 200) == "exclude_interval_overlap_outside_POS_window"
    assert MODULE.reference_record_position_decision("9", 100, 100, "9", 100, 200) == "include"
    with pytest.raises(RuntimeError, match="E_REMOTE_REFERENCE_RECORD_OUTSIDE_QUERY"):
        MODULE.reference_record_position_decision("9", 80, 90, "9", 100, 200)
    with pytest.raises(RuntimeError, match="E_REMOTE_REFERENCE_CONTIG_MISMATCH"):
        MODULE.reference_record_position_decision("8", 100, 100, "9", 100, 200)


def test_discordant_reference_key_is_permanently_quarantined():
    registry, ambiguous, details = {}, set(), {}
    key = (10, frozenset(("A", "AAAAAC")))
    first = {"ref": "A", "alt": "AAAAAC", "dosage": np.array([0.0, 1.0]), "audit_descriptor": {"ref": "A", "alt": "AAAAAC"}}
    reciprocal = {"ref": "AAAAAC", "alt": "A", "dosage": np.array([0.0, 2.0]), "audit_descriptor": {"ref": "AAAAAC", "alt": "A"}}
    assert MODULE.register_reference_key(registry, ambiguous, details, key, first) == "new"
    assert MODULE.register_reference_key(registry, ambiguous, details, key, reciprocal) == "discordant_key_permanently_quarantined"
    assert key not in registry and key in ambiguous and len(details[key]) == 2
    assert MODULE.register_reference_key(registry, ambiguous, details, key, first) == "ambiguous_additional_record_quarantined"
    assert key not in registry and len(details[key]) == 3
