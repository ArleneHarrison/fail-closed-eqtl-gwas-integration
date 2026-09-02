from __future__ import annotations

import argparse
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("external_validation", ROOT / "run_external_validation.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)
PROTOCOL = ROOT / "config" / "protocol.json"
PROTOCOL_V2 = ROOT / "config" / "protocol.v2.json"


class FrozenProtocolTests(unittest.TestCase):
    def test_protocol_is_valid_and_has_expected_factorial(self):
        protocol = MODULE.read_json(PROTOCOL)
        MODULE.validate_protocol(protocol)
        self.assertEqual(MODULE.validate_frozen_protocol_path(PROTOCOL), MODULE.protocol_hash(protocol))
        self.assertEqual(protocol["factorial_design"]["expected_attempts"], 192)
        self.assertEqual(protocol["externality_statement"]["forbidden_label"], "fully_independent_replication")

    def test_protocol_rejects_development_cad_as_primary(self):
        protocol = MODULE.read_json(PROTOCOL)
        protocol["outcomes"][0]["dataset_id"] = "GCST90132314"
        with self.assertRaisesRegex(ValueError, "GCST005195"):
            MODULE.validate_protocol(protocol)

    def test_amended_v2_protocol_and_metadata_are_locked(self):
        protocol = MODULE.read_json(PROTOCOL_V2)
        self.assertEqual(MODULE.validate_frozen_protocol_path(PROTOCOL_V2), MODULE.protocol_hash(protocol))
        lock = MODULE.validate_metadata_lock(ROOT / "config" / "metadata.locked.v2.json", protocol)
        self.assertEqual(lock["outcomes"]["CAD"]["cases"], 34541)
        self.assertEqual(lock["outcomes"]["HF"]["cases"], 139533)


class ProspectiveSelectionTests(unittest.TestCase):
    def test_lead_selection_uses_declared_columns_and_fixed_spacing(self):
        rows = []
        for index in range(24):
            rows.append({
                "chr": "1",
                "bp": 1_000_000 + index * 1_100_000,
                "p": (index + 1) * 1e-10,
                "snp": f"rs{index + 1}",
                "downstream_posterior": 1 - index / 100,
            })
        leads = MODULE.select_leads(
            pd.DataFrame(rows), outcome="CAD", chromosome_col="chr", position_col="bp",
            p_col="p", variant_col="snp", n_leads=12, p_threshold=5e-8,
            min_distance_bp=1_000_000, window_bp=500_000,
        )
        self.assertEqual(len(leads), 12)
        self.assertEqual([row["lead_variant"] for row in leads], [f"rs{i}" for i in range(1, 13)])
        self.assertTrue(all(b["lead_position"] - a["lead_position"] >= 1_000_000 for a, b in zip(leads, leads[1:])))

    def test_lead_shortfall_is_fail_closed(self):
        frame = pd.DataFrame({"chr": ["1"], "bp": [100], "p": [1e-9], "snp": ["rs1"]})
        with self.assertRaisesRegex(RuntimeError, "E_LEAD_LOCUS_SHORTFALL"):
            MODULE.select_leads(
                frame, outcome="HF", chromosome_col="chr", position_col="bp", p_col="p",
                variant_col="snp", n_leads=12, p_threshold=5e-8,
                min_distance_bp=1_000_000, window_bp=500_000,
            )

    def test_gene_choice_is_lexicographic_and_not_effect_ranked(self):
        frame = pd.DataFrame({
            "gene_id": ["ENSG000002", "ENSG000001"] * 20,
            "chromosome": ["1"] * 40,
            "position": list(range(110, 150)),
            "variant": [f"v{i}" for i in range(40)],
            "ref": ["A"] * 40,
            "alt": ["G"] * 40,
            "beta": [1000, -0.0001] * 20,
            "se": [1.0] * 40,
            "n": [100] * 40,
        })
        region = pd.Series({"chromosome": "1", "region_start": 100, "region_end": 200})
        required = ["gene_id", "chromosome", "position", "variant", "ref", "alt", "beta", "se", "n"]
        candidates = MODULE.complete_gene_candidates(frame, region, required, 20)
        self.assertEqual(sorted(candidates)[0], "ENSG000001")


class FailClosedGateTests(unittest.TestCase):
    def test_summary_detects_hash_bound_nested_model_result(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            attempts = root / "attempts.tsv"
            gates = root / "gates"
            models = root / "models"
            out = root / "summary.tsv"
            pd.DataFrame([{
                "unit_id": "u_ready", "region_id": "CAD_L01", "anchor_outcome": "CAD",
                "tissue_id": "QTD000131", "outcome_id": "CAD",
                "pre_gate_status": "PENDING", "pre_gate_code": "PENDING",
            }]).to_csv(attempts, sep="\t", index=False)
            gate_dir = gates / "units" / "u_ready"
            gate_dir.mkdir(parents=True)
            (gate_dir / "gate.json").write_text(json.dumps({
                "status": "READY", "status_codes": ["OK"], "n_summary_rows": 2,
                "ld_diagnostics": {"n_variants": 2, "numerical_rank": 2, "condition_number_2": 1.0},
                "overlap_diagnostics": {
                    "n_ordered_overlap": 2, "n_eqtl_variants_pre": 4, "n_gwas_variants_pre": 5,
                },
            }), encoding="utf-8")
            result_dir = models / "units" / "u_ready" / "run_bindinghash"
            result_dir.mkdir(parents=True)
            (result_dir / "coloc_susie_result.rds").write_bytes(b"test-result")
            MODULE.cmd_summarize(argparse.Namespace(
                attempts=attempts, gates=gates, models=models, out=out,
            ))
            summary = pd.read_csv(out, sep="\t")
            self.assertTrue(bool(summary.loc[0, "posterior_generated"]))

    def test_r_runner_uses_trait_specific_maf_columns(self):
        runner = (ROOT / "run_coloc_susie_rank_diagnostic.R").read_text(encoding="utf-8")
        self.assertIn('"eqtl_maf" %in% names(summary_data)', runner)
        self.assertIn('"gwas_maf" %in% names(summary_data)', runner)

    def test_r_runner_has_explicit_version_and_susie_implementation_lock(self):
        runner = (ROOT / "run_coloc_susie_rank_diagnostic.R").read_text(encoding="utf-8")
        for literal in (
            'EXPECTED_COLOC_VERSION <- "5.2.3"',
            'EXPECTED_SUSIER_VERSION <- "0.14.2"',
            "SUSIE_L <- 10L",
            "SUSIE_MAX_ITER <- 1000L",
            "SUSIE_REPEAT_UNTIL_CONVERGENCE <- FALSE",
            "SUSIE_ESTIMATE_RESIDUAL_VARIANCE <- FALSE",
            "SUSIE_TOL <- 0.001",
            "SUSIE_COVERAGE <- 0.95",
            "SUSIE_MIN_ABS_CORR <- 0.5",
            "SUSIE_SCALED_PRIOR_VARIANCE <- 0.2",
            "SUSIE_CHECK_PRIOR <- TRUE",
            "SUSIE_Z_LD_WEIGHT <- 0",
            "COLOC_P1 <- 1e-4",
            "COLOC_P2 <- 1e-4",
            "COLOC_P12 <- 5e-6",
            "COLOC_OVERLAP_MIN <- 0.5",
            "COLOC_TRIM_BY_POSTERIOR <- TRUE",
            "COLOC_BACK_CALCULATE_LBF <- FALSE",
            "repeat_until_convergence = SUSIE_REPEAT_UNTIL_CONVERGENCE",
            "scaled_prior_variance = SUSIE_SCALED_PRIOR_VARIANCE",
            "check_prior = SUSIE_CHECK_PRIOR",
            "z_ld_weight = SUSIE_Z_LD_WEIGHT",
            "back_calculate_lbf = COLOC_BACK_CALCULATE_LBF",
            "overlap.min = COLOC_OVERLAP_MIN",
            "trim_by_posterior = COLOC_TRIM_BY_POSTERIOR",
            "p12 = COLOC_P12",
        ):
            self.assertIn(literal, runner)

    def test_gate_dispatch_marker_is_distinct_and_bound_by_dispatcher(self):
        gate_script = (ROOT / "run_external_validation.py").read_text(encoding="utf-8")
        dispatcher = (ROOT / "dispatch_ready_models.py").read_text(encoding="utf-8")
        self.assertNotIn("model_eligible", gate_script)
        self.assertNotIn("model_eligible", dispatcher)
        self.assertIn('result["gate_ready_for_dispatch"] = result["status"] == "READY"', gate_script)
        self.assertIn('gate.get("gate_ready_for_dispatch") is not True', dispatcher)
        self.assertIn('"gate_ready_for_dispatch": True', dispatcher)
        self.assertIn('"scaled_prior_variance": lock["scaled_prior_variance"]', dispatcher)

    def test_ready_and_ineligible_units_are_recorded_without_model_dispatch(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            prepared = root / "prepared"
            out = root / "gates"
            unit = prepared / "u_ready"
            unit.mkdir(parents=True)
            pd.DataFrame([
                {"target_variant": "v1", "eqtl_se": 0.1, "gwas_se": 0.2, "maf": 0.2,
                 "eqtl_n": 100, "gwas_n": 1000, "alignment_status": "aligned", "vcf_status": "present"},
                {"target_variant": "v2", "eqtl_se": 0.1, "gwas_se": 0.2, "maf": 0.3,
                 "eqtl_n": 100, "gwas_n": 1000, "alignment_status": "flipped", "vcf_status": "present"},
            ]).to_csv(unit / "summary.csv", index=False)
            np.savez_compressed(unit / "ld.npz", ld=np.eye(2), variant_ids=np.array(["v1", "v2"]))
            ordered_hash = MODULE.sha256_bytes(b"v1\nv2\n")
            (unit / "overlap.json").write_text(json.dumps({
                "n_eqtl_variants_pre": 3, "n_gwas_variants_pre": 4,
                "n_ld_variants_pre": 2, "n_ordered_overlap": 2,
                "ordered_variant_sha256": ordered_hash,
            }), encoding="utf-8")
            provenance = [{
                "label": "synthetic", "url": "local:test", "version_or_build": "test",
                "license_or_access": "test", "retrieval_date": "2026-09-01",
                "sha256": "0" * 64, "byte_count": 1, "structural_check": "passed",
            }]
            (unit / "provenance.json").write_text(json.dumps(provenance), encoding="utf-8")
            empty = prepared / "u_empty"
            empty.mkdir()
            pd.DataFrame(columns=["target_variant"]).to_csv(empty / "summary.csv", index=False)
            np.savez_compressed(empty / "ld.npz", ld=np.empty((0, 0)), variant_ids=np.array([], dtype=str))
            (empty / "provenance.json").write_text(json.dumps(provenance), encoding="utf-8")
            (empty / "overlap.json").write_text(json.dumps({
                "n_eqtl_variants_pre": 3, "n_gwas_variants_pre": 4,
                "n_ld_variants_pre": 2, "n_ordered_overlap": 0,
                "ordered_variant_sha256": MODULE.sha256_bytes(b"\n"),
            }), encoding="utf-8")
            attempts = root / "attempts.tsv"
            pd.DataFrame([
                {"unit_id": "u_ready", "outcome_id": "CAD", "pre_gate_status": "PENDING", "pre_gate_code": "PENDING"},
                {"unit_id": "u_stop", "outcome_id": "HF", "pre_gate_status": "INELIGIBLE", "pre_gate_code": "E_NO_ELIGIBLE_GENE"},
                {"unit_id": "u_empty", "outcome_id": "CAD", "pre_gate_status": "PENDING", "pre_gate_code": "PENDING"},
            ]).to_csv(attempts, sep="\t", index=False)
            metadata_lock = root / "metadata.lock.json"
            metadata_lock.write_text(json.dumps({
                "state": "LOCKED",
                "protocol_canonical_json_sha256": MODULE.protocol_hash(MODULE.read_json(PROTOCOL)),
                "outcomes": {
                    "CAD": {"dataset_id": "GCST005195", "ancestry": "test", "genome_build": "GRCh37",
                            "cases": 400, "controls": 600, "sample_size": 1000, "case_fraction": 0.4,
                            "metadata_url": "local:test", "source_citation": "test"},
                    "HF": {"dataset_id": "HERMES2_EUR_2025", "ancestry": "European", "genome_build": "GRCh37",
                           "cases": 100, "controls": 900, "sample_size": 1000, "case_fraction": 0.1,
                           "metadata_url": "local:test", "source_citation": "test"},
                },
            }), encoding="utf-8")
            args = argparse.Namespace(
                protocol=PROTOCOL, attempts=attempts, prepared=prepared,
                metadata_lock=metadata_lock, out=out,
            )
            MODULE.cmd_gate(args)
            ready = MODULE.read_json(out / "units/u_ready/gate.json")
            self.assertEqual(ready["status"], "READY")
            self.assertTrue(ready["gate_ready_for_dispatch"])
            self.assertIsNone(ready["model_eligibility"])
            stop = MODULE.read_json(out / "units/u_stop/gate.json")
            self.assertEqual(stop["status"], "INELIGIBLE")
            self.assertEqual(stop["status_codes"], ["E_NO_ELIGIBLE_GENE"])
            self.assertFalse(stop["gate_ready_for_dispatch"])
            self.assertFalse(any((out / "units/u_stop").glob("*.rds")))
            self.assertEqual(MODULE.read_json(out / "units/u_empty/gate.json")["status_codes"], ["E_NO_ORDERED_OVERLAP"])


if __name__ == "__main__":
    unittest.main()
