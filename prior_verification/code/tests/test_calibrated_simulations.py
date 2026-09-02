import unittest

from hcsmr.reanalysis.calibrated_simulations import (
    run_calibrated_replicates,
    summarize_method_performance,
)
from hcsmr.reanalysis.simulation_protocol import SimulationScenarioV2


class CalibratedSimulationTests(unittest.TestCase):
    def test_preserves_adverse_scenario_labels_and_method_eligibility(self):
        scenario = SimulationScenarioV2(3, "directional", "residual_ld", 0.25)

        results = run_calibrated_replicates(scenario, replicates=4, true_theta=0.0)

        self.assertEqual(len(results), 12)
        self.assertSetEqual(set(results["pleiotropy"]), {"directional"})
        self.assertSetEqual(set(results["ld"]), {"residual_ld"})
        self.assertTrue((results.loc[results["method"] == "egger", "eligible"] == False).all())

    def test_summary_reports_wilson_intervals_and_failures(self):
        scenario = SimulationScenarioV2(5, "balanced", "independent", 0.0)
        results = run_calibrated_replicates(scenario, replicates=8, true_theta=0.0)

        summary = summarize_method_performance(results, true_theta=0.0)

        row = summary.loc[summary["method"] == "ivw"].iloc[0]
        self.assertEqual(row["n_replicates"], 8)
        self.assertIn("type1_ci_low", summary.columns)
        self.assertIn("ineligible_replicates", summary.columns)
