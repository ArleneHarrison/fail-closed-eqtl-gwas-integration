import unittest

import numpy as np

from hcsmr.reanalysis.preflight_gate import prepare_analysis_rows, run_preflight


PROVENANCE = [{
    "label": "synthetic_input", "url": "https://example.org/synthetic", "version_or_build": "test-v1",
    "license_or_access": "CC0", "retrieval_date": "2026-08-26", "sha256": "0" * 64,
    "byte_count": 1, "structural_check": "passed",
}]
ROW = {
    "target_variant": "chr1_1_A_C", "eqtl_se": "0.1", "gwas_se": "0.2", "maf": "0.2",
    "an": "100", "gwas_n": "1000", "alignment_status": "aligned", "vcf_status": "exact_ref_alt_match",
}


class PreflightGateTests(unittest.TestCase):
    def call(self, rows=None, ld=None, provenance=PROVENANCE):
        active_rows = rows or [ROW]
        return run_preflight(
            rows=active_rows, ld=np.eye(1) if ld is None else ld, provenance=provenance,
            eqtl_trait_type="quant", gwas_trait_type="cc", case_fraction=0.2,
            max_condition_number=1e12, ld_variant_ids=[row["target_variant"] for row in active_rows],
        )

    def test_ready_for_valid_fixture(self):
        result = self.call()
        self.assertEqual(result["status"], "READY")
        self.assertEqual(result["status_codes"], ["OK"])

    def test_maf_upper_boundary_is_inclusive_and_values_above_it_fail(self):
        boundary = self.call(rows=[dict(ROW, maf="0.5")])
        self.assertEqual(boundary["status"], "READY")
        above_boundary = self.call(rows=[dict(ROW, maf="0.5001")])
        self.assertIn("E_MAF_INVALID", above_boundary["status_codes"])

    def test_unmapped_one_to_many_strand_and_indel_are_fail_closed(self):
        rows = []
        for index, status in enumerate(("liftover_unmapped", "liftover_one_to_many", "strand_unresolved", "indel_unresolved")):
            row = dict(ROW, target_variant=f"chr1_{index + 1}_A_C", alignment_status=status)
            rows.append(row)
        result = self.call(rows=rows)
        self.assertEqual(result["status"], "INELIGIBLE")
        self.assertIn("E_ALLELE_UNRESOLVED", result["status_codes"])
        self.assertEqual(result["failures"][0]["affected_rows"], 4)

    def test_duplicate_and_vcf_absence_are_fail_closed(self):
        absent = dict(ROW, vcf_status="position_absent")
        result = self.call(rows=[ROW, absent])
        self.assertIn("E_VARIANT_DUPLICATE", result["status_codes"])
        self.assertIn("E_VCF_ABSENT", result["status_codes"])

    def test_ill_conditioned_ld_is_fail_closed_without_repair(self):
        matrix = np.array([[1.0, 1.0 - 1e-14], [1.0 - 1e-14, 1.0]])
        rows = [ROW, dict(ROW, target_variant="chr1_2_A_C")]
        result = self.call(rows=rows, ld=matrix)
        self.assertIn("E_LD_ILL_CONDITIONED", result["status_codes"])
        self.assertGreater(result["ld_diagnostics"]["condition_number_2"], 1e12)

    def test_rank_deficient_ld_is_reported_separately(self):
        matrix = np.ones((3, 3))
        rows = [dict(ROW, target_variant=f"chr1_{index + 1}_A_C") for index in range(3)]
        result = self.call(rows=rows, ld=matrix)
        self.assertIn("E_LD_RANK_DEFICIENT", result["status_codes"])
        self.assertTrue(result["ld_diagnostics"]["rank_deficient"])
        self.assertIsNone(result["ld_diagnostics"]["condition_number_2"])

    def test_rank_deficiency_can_be_diagnostic_under_a_named_policy(self):
        matrix = np.ones((3, 3))
        rows = [dict(ROW, target_variant=f"chr1_{index + 1}_A_C") for index in range(3)]
        result = run_preflight(
            rows=rows, ld=matrix, provenance=PROVENANCE,
            eqtl_trait_type="quant", gwas_trait_type="cc", case_fraction=0.2,
            max_condition_number=1e12, ld_rank_policy="diagnostic",
            ld_variant_ids=[row["target_variant"] for row in rows],
        )
        self.assertEqual(result["status"], "READY")
        self.assertIsNone(result["model_eligibility"])
        self.assertTrue(result["ld_diagnostics"]["rank_deficient"])
        self.assertFalse(result["ld_diagnostics"]["condition_number_threshold_applied"])

    def test_indefinite_correlation_like_matrix_is_rejected(self):
        matrix = np.array([
            [1.0, -0.6, -0.6],
            [-0.6, 1.0, -0.6],
            [-0.6, -0.6, 1.0],
        ])
        rows = [dict(ROW, target_variant=f"chr1_{index + 1}_A_C") for index in range(3)]
        result = self.call(rows=rows, ld=matrix)
        self.assertIn("E_LD_NOT_PSD", result["status_codes"])
        self.assertFalse(result["ld_diagnostics"]["positive_semidefinite"])

    def test_requires_complete_provenance_and_valid_case_fraction(self):
        result = run_preflight(
            rows=[ROW], ld=np.eye(1), provenance=[{"label": "missing"}],
            eqtl_trait_type="quant", gwas_trait_type="cc", case_fraction=1.0,
            max_condition_number=1e12, ld_variant_ids=[ROW["target_variant"]],
        )
        self.assertIn("E_PROVENANCE_INCOMPLETE", result["status_codes"])
        self.assertIn("E_CASE_FRACTION_INVALID", result["status_codes"])

    def test_rejects_heterogeneous_row_case_fractions(self):
        rows = [
            dict(ROW, gwas_cases="100", gwas_n="1000"),
            dict(ROW, target_variant="chr1_2_A_C", gwas_cases="200", gwas_n="1000"),
        ]
        result = self.call(rows=rows, ld=np.eye(2))
        self.assertIn("E_CASE_FRACTION_HETEROGENEOUS", result["status_codes"])

    def test_study_level_case_fraction_treats_row_heterogeneity_as_diagnostic(self):
        rows = [
            dict(ROW, gwas_cases="100", gwas_n="1000"),
            dict(ROW, target_variant="chr1_2_A_C", gwas_cases="200", gwas_n="1000"),
        ]
        result = run_preflight(
            rows=rows, ld=np.eye(2), provenance=PROVENANCE,
            eqtl_trait_type="quant", gwas_trait_type="cc", case_fraction=0.2,
            case_fraction_policy="study_level", max_condition_number=1e12,
            ld_variant_ids=[row["target_variant"] for row in rows],
        )
        self.assertNotIn("E_CASE_FRACTION_HETEROGENEOUS", result["status_codes"])
        self.assertEqual(result["status"], "READY")

    def test_rejects_ld_variant_set_mismatch(self):
        result = run_preflight(
            rows=[ROW], ld=np.eye(1), provenance=PROVENANCE,
            eqtl_trait_type="quant", gwas_trait_type="cc", case_fraction=0.2,
            max_condition_number=1e12, ld_variant_ids=["chr1_2_A_C"],
        )
        self.assertIn("E_LD_VARIANT_SET_MISMATCH", result["status_codes"])

    def test_rejects_ld_identifier_duplicates_count_and_order_mismatch(self):
        rows = [ROW, dict(ROW, target_variant="chr1_2_A_C")]
        common = dict(
            rows=rows, ld=np.eye(2), provenance=PROVENANCE,
            eqtl_trait_type="quant", gwas_trait_type="cc", case_fraction=0.2,
            max_condition_number=1e12,
        )
        duplicated = run_preflight(**common, ld_variant_ids=[ROW["target_variant"], ROW["target_variant"]])
        self.assertIn("E_LD_VARIANT_DUPLICATE", duplicated["status_codes"])
        shortened = run_preflight(**common, ld_variant_ids=[ROW["target_variant"]])
        self.assertIn("E_LD_VARIANT_COUNT_MISMATCH", shortened["status_codes"])
        reordered = run_preflight(**common, ld_variant_ids=[rows[1]["target_variant"], rows[0]["target_variant"]])
        self.assertIn("E_LD_VARIANT_ORDER_MISMATCH", reordered["status_codes"])

    def test_analysis_view_collapses_only_exact_duplicates_and_audits_exclusions(self):
        duplicate = dict(ROW)
        absent = dict(ROW, target_variant="chr1_9_A_C")
        selected, audit = prepare_analysis_rows(
            [ROW, duplicate, absent], ld_variant_ids=[ROW["target_variant"]],
        )
        self.assertEqual(selected, [ROW])
        self.assertEqual(
            [item["decision"] for item in audit],
            ["retained_first_occurrence", "excluded_exact_duplicate", "excluded_not_in_locked_ld"],
        )

    def test_analysis_view_does_not_choose_between_discordant_duplicates(self):
        discordant = dict(ROW, eqtl_se="0.3")
        selected, audit = prepare_analysis_rows(
            [ROW, discordant], ld_variant_ids=[ROW["target_variant"]],
        )
        self.assertEqual(len(selected), 2)
        self.assertEqual(audit[-1]["decision"], "retained_discordant_duplicate_for_gate_failure")
