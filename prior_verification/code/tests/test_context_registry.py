import unittest

from hcsmr.reanalysis.context_registry import ContextRecord, validate_context_pair


class ContextRegistryTests(unittest.TestCase):
    def test_rejects_cardiac_localization_for_pbmc_and_bulk_contexts(self):
        pbmc = ContextRecord(
            "OneK1K_NK",
            "immune_population",
            "PBMC",
            False,
            "OneK1K",
        )
        bulk = ContextRecord(
            "GTEx_Heart_LV",
            "bulk_tissue",
            "heart_left_ventricle",
            False,
            "GTEx v8",
        )

        decision = validate_context_pair([pbmc, bulk], pbmc.name, bulk.name)

        self.assertFalse(decision.can_claim_cardiac_cell_state_localization)
        self.assertEqual(decision.analysis_mode, "separate")

    def test_permits_joint_analysis_only_for_related_contexts_with_covariance(self):
        left = ContextRecord("Fairfax_LPS2", "immune_stimulation", "monocyte", True, "Fairfax")
        right = ContextRecord("Fairfax_LPS24", "immune_stimulation", "monocyte", True, "Fairfax")

        decision = validate_context_pair([left, right], left.name, right.name)

        self.assertEqual(decision.analysis_mode, "joint")
        self.assertFalse(decision.can_claim_cardiac_cell_state_localization)
