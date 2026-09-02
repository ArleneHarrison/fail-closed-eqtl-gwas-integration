import unittest
from dataclasses import dataclass

import numpy as np

from hcsmr.reanalysis.benchmarks import run_comparator_suite


@dataclass
class ComparatorData:
    beta_hat: np.ndarray
    gamma_hat: np.ndarray
    se_x: np.ndarray
    se_y: np.ndarray

    @property
    def n_instruments(self) -> int:
        return len(self.beta_hat)


class BenchmarkTests(unittest.TestCase):
    def test_records_egger_ineligibility_for_three_instruments(self):
        data = ComparatorData(
            beta_hat=np.array([0.2, 0.3, 0.4]),
            gamma_hat=np.array([0.02, 0.03, 0.04]),
            se_x=np.array([0.02, 0.02, 0.02]),
            se_y=np.array([0.01, 0.01, 0.01]),
        )

        result = run_comparator_suite(data, methods=["ivw", "egger", "weighted_median"])
        egger = result.loc[result["method"].eq("egger")].iloc[0]

        self.assertFalse(egger["eligible"])
        self.assertEqual(egger["ineligibility_reason"], "fewer_than_four_instruments")

    def test_returns_standard_columns_for_supported_methods(self):
        data = ComparatorData(
            beta_hat=np.array([0.2, 0.3, 0.4, 0.5]),
            gamma_hat=np.array([0.02, 0.03, 0.04, 0.05]),
            se_x=np.array([0.02] * 4),
            se_y=np.array([0.01] * 4),
        )

        result = run_comparator_suite(data, methods=["ivw", "weighted_median"])

        self.assertEqual(
            set(result.columns),
            {"method", "estimate", "standard_error", "p_value", "eligible", "ineligibility_reason"},
        )
        self.assertTrue(result["eligible"].all())

    def test_weighted_median_is_reproducible(self):
        data = ComparatorData(
            beta_hat=np.array([0.1, 0.15, 0.2, 0.3]),
            gamma_hat=np.array([0.02, 0.03, 0.04, 0.06]),
            se_x=np.array([0.01, 0.01, 0.01, 0.01]),
            se_y=np.array([0.02, 0.02, 0.02, 0.02]),
        )

        first = run_comparator_suite(data, methods=["weighted_median"])
        second = run_comparator_suite(data, methods=["weighted_median"])

        self.assertEqual(first["estimate"].item(), second["estimate"].item())
        self.assertEqual(first["standard_error"].item(), second["standard_error"].item())
