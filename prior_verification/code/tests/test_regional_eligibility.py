import unittest

import pandas as pd

from hcsmr.reanalysis.regional_eligibility import assess_regional_eligibility


class RegionalEligibilityTests(unittest.TestCase):
    def test_marks_disjoint_contexts_as_not_localizable(self):
        variants = pd.DataFrame(
            {
                "context": ["A", "A", "A", "B", "B", "B"],
                "variant_id": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6"],
                "f_stat": [20.0] * 6,
            }
        )

        decision = assess_regional_eligibility(variants, min_shared_variants=3)

        self.assertFalse(decision.joint_fit_allowed)
        self.assertEqual(decision.reason, "insufficient_shared_variants")

    def test_marks_low_f_stat_context_as_ineligible(self):
        variants = pd.DataFrame(
            {
                "context": ["A", "A", "A", "B", "B", "B"],
                "variant_id": ["rs1", "rs2", "rs3", "rs1", "rs2", "rs3"],
                "f_stat": [20.0, 20.0, 20.0, 20.0, 9.9, 20.0],
            }
        )

        decision = assess_regional_eligibility(variants)

        self.assertFalse(decision.joint_fit_allowed)
        self.assertIn("weak_instruments", decision.reason)
