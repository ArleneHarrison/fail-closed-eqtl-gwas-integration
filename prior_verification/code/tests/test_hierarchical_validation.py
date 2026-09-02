from __future__ import annotations

import unittest

from hcsmr.reanalysis.hierarchical_validation import required_hierarchical_validation_scenarios


class HierarchicalValidationProtocolTests(unittest.TestCase):
    def test_grid_covers_null_sparse_dense_and_pleiotropic_conditions(self):
        scenarios = required_hierarchical_validation_scenarios()
        names = {scenario.name for scenario in scenarios}
        self.assertIn("null_sparse_clean", names)
        self.assertIn("signal_dense_pleiotropy", names)
        self.assertTrue(any(scenario.pleiotropy_fraction > 0 for scenario in scenarios))
        self.assertTrue(any(scenario.ivs_per_state == 3 for scenario in scenarios))
        self.assertTrue(any(scenario.ivs_per_state == 10 for scenario in scenarios))


if __name__ == "__main__":
    unittest.main()
