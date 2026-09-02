from __future__ import annotations

import unittest

import pandas as pd

from hcsmr.reanalysis.candidate_sensitivity import candidate_sensitivity, ivw_summary


def _variants(n: int) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "candidate_id": ["CAD|GENE|context"] * n,
            "rsid": [f"rs{i}" for i in range(n)],
            "beta": [0.20, 0.17, 0.13, 0.11][:n],
            "se": [0.02] * n,
            "gwas_beta": [0.03, 0.026, 0.018, 0.017][:n],
            "gwas_se": [0.01] * n,
            "f_stat": [100.0, 72.25, 42.25, 30.25][:n],
        }
    )


class CandidateSensitivityTests(unittest.TestCase):
    def test_ivw_reports_heterogeneity_fields(self):
        result = ivw_summary(_variants(4))
        self.assertEqual(result["n_iv"], 4)
        self.assertIn("q_statistic", result)
        self.assertIn("q_p_value", result)

    def test_four_signals_produce_leave_one_out_rows(self):
        summary, leave_one_out = candidate_sensitivity(_variants(4))
        self.assertEqual(len(leave_one_out), 4)
        self.assertIn(summary.iloc[0]["leave_one_out_status"], {"stable", "unstable"})
        self.assertEqual(summary.iloc[0]["instrument_strength_status"], "pass")

    def test_three_signals_do_not_claim_leave_one_out_assessment(self):
        summary, leave_one_out = candidate_sensitivity(_variants(3))
        self.assertTrue(leave_one_out.empty)
        self.assertEqual(
            summary.iloc[0]["leave_one_out_status"],
            "not_assessed_fewer_than_four_signals",
        )


if __name__ == "__main__":
    unittest.main()
