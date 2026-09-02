import unittest

import pandas as pd

from hcsmr.reanalysis.coordinate_contract import (
    SourceCoordinateRecord,
    build_coordinate_audit,
    validate_source_pair,
)


class CoordinateContractTests(unittest.TestCase):
    def setUp(self):
        self.registry = [
            SourceCoordinateRecord("GTEx v8", "GRCh38", "eQTL"),
            SourceCoordinateRecord("CAD Aragam 2022", "GRCh37", "GWAS"),
            SourceCoordinateRecord("HERMES HF", "unknown", "GWAS"),
        ]

    def test_requires_liftover_for_declared_build_mismatch(self):
        result = validate_source_pair(self.registry, "GTEx v8", "CAD Aragam 2022")

        self.assertEqual(result.status, "needs_liftover")
        self.assertFalse(result.coordinate_pass)

    def test_retains_unknown_build_as_an_explicit_failure(self):
        result = validate_source_pair(self.registry, "GTEx v8", "HERMES HF")

        self.assertEqual(result.status, "unknown_build")
        self.assertFalse(result.coordinate_pass)

    def test_audit_assigns_sources_from_context_and_outcome(self):
        candidates = pd.DataFrame(
            {"candidate_id": ["CAD|ENSG1|GTEx_Fibroblast"]}
        )

        audit = build_coordinate_audit(candidates, self.registry)

        self.assertEqual(audit.iloc[0]["status"], "needs_liftover")
        self.assertIn("liftover", audit.iloc[0]["required_action"])
