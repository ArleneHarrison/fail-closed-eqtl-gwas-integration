import unittest

import pandas as pd

from hcsmr.reanalysis.candidate_confirmation import build_confirmation_table
from hcsmr.reanalysis.context_registry import ContextRecord


class CandidateConfirmationTableTests(unittest.TestCase):
    def setUp(self):
        self.contexts = [
            ContextRecord(
                name="OneK1K_NK",
                evidence_class="immune_population",
                biological_source="PBMC",
                covariance_available=False,
                source_study="OneK1K",
            )
        ]

    def test_coloc_failure_prevents_primary_candidate(self):
        diagnostics = pd.DataFrame(
            {"candidate_id": ["CAD|ENSG1|OneK1K_NK"], "instrument_pass": [True]}
        )
        coloc = pd.DataFrame(
            {
                "region_id": ["CAD|ENSG1|OneK1K_NK"],
                "status": ["ineligible"],
                "ineligibility_reason": ["missing_ld_reference"],
            }
        )

        result = build_confirmation_table(diagnostics, self.contexts, coloc)

        self.assertFalse(result.iloc[0]["primary_candidate"])
        self.assertIn("formal_colocalization", result.iloc[0]["failed_criteria"])
        self.assertEqual(result.iloc[0]["context_evidence_class"], "immune_population")

    def test_weak_instrument_is_explicit_failure(self):
        diagnostics = pd.DataFrame(
            {"candidate_id": ["CAD|ENSG1|OneK1K_NK"], "instrument_pass": [False]}
        )
        coloc = pd.DataFrame(columns=["region_id", "status", "ineligibility_reason"])

        result = build_confirmation_table(diagnostics, self.contexts, coloc)

        self.assertFalse(result.iloc[0]["primary_candidate"])
        self.assertIn("instrument_strength", result.iloc[0]["failed_criteria"])
        self.assertEqual(result.iloc[0]["coloc_status"], "missing_coloc_record")

    def test_uses_coordinate_audit_when_supplied(self):
        diagnostics = pd.DataFrame(
            {"candidate_id": ["CAD|ENSG1|OneK1K_NK"], "instrument_pass": [True]}
        )
        coordinates = pd.DataFrame(
            {"candidate_id": ["CAD|ENSG1|OneK1K_NK"], "coordinate_pass": [True]}
        )
        coloc = pd.DataFrame(columns=["region_id", "status", "ineligibility_reason"])

        result = build_confirmation_table(diagnostics, self.contexts, coloc, coordinates)

        self.assertTrue(result.iloc[0]["coordinate_pass"])
