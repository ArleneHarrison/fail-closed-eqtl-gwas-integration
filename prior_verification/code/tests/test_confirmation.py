import unittest

from hcsmr.reanalysis.confirmation import evaluate_candidate, validate_manuscript_values
from hcsmr.reanalysis.manifest import validate_manifest


class ConfirmationTests(unittest.TestCase):
    def setUp(self):
        self.valid_candidate = {
            "evidence_class_pass": True,
            "coordinate_pass": True,
            "instrument_pass": True,
            "convergence_pass": True,
            "sensitivity_pass": True,
            "multiplicity_pass": True,
            "coloc_pass": True,
        }

    def test_fails_candidate_without_formal_colocalization(self):
        candidate = {**self.valid_candidate, "coloc_pass": False}

        decision = evaluate_candidate(candidate)

        self.assertFalse(decision.primary_candidate)
        self.assertIn("formal_colocalization", decision.failed_criteria)

    def test_fails_candidate_without_coordinate_harmonization(self):
        candidate = {**self.valid_candidate, "coordinate_pass": False}

        decision = evaluate_candidate(candidate)

        self.assertFalse(decision.primary_candidate)
        self.assertIn("coordinate_harmonization", decision.failed_criteria)

    def test_rejects_manifest_without_input_hashes(self):
        manifest = {"locked_run_id": "run-001", "input_hashes": {}, "analysis_version": "v2"}

        decision = validate_manifest(manifest)

        self.assertFalse(decision.valid)
        self.assertIn("input_hashes", decision.errors)

    def test_rejects_manuscript_values_without_locked_provenance(self):
        decision = validate_manuscript_values({"headline_hit_count": 0})

        self.assertFalse(decision.valid)
