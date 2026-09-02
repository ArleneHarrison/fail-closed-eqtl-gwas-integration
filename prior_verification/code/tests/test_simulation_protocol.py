import unittest

from hcsmr.reanalysis.simulation_protocol import binomial_interval, required_scenarios


class SimulationProtocolTests(unittest.TestCase):
    def test_required_scenarios_include_sparse_and_dense_instruments(self):
        iv_counts = {scenario.n_instruments for scenario in required_scenarios()}

        self.assertTrue({3, 4, 5, 10, 20}.issubset(iv_counts))

    def test_required_scenarios_include_directional_and_inside_violating_pleiotropy(self):
        mechanisms = {scenario.pleiotropy for scenario in required_scenarios()}

        self.assertIn("directional", mechanisms)
        self.assertIn("inside_violation", mechanisms)

    def test_binomial_interval_is_nonzero_width_for_zero_of_1000(self):
        low, high = binomial_interval(0, 1000)

        self.assertEqual(low, 0.0)
        self.assertGreater(high, 0.0)
        self.assertLess(high, 0.01)

    def test_binomial_interval_rejects_invalid_counts(self):
        with self.assertRaises(ValueError):
            binomial_interval(3, 2)
