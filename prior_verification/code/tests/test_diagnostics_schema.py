import unittest

import pandas as pd

from hcsmr.reanalysis.diagnostics import diagnostic_row, with_f_statistics


class DiagnosticsSchemaTests(unittest.TestCase):
    def test_diagnostic_row_has_every_required_field(self):
        variants = pd.DataFrame(
            {
                "variant_id": ["rs1", "rs2", "rs3"],
                "f_stat": [11.0, 20.0, 30.0],
            }
        )

        row = diagnostic_row("candidate-1", variants)

        required = {
            "candidate_id",
            "n_iv",
            "min_f_stat",
            "mean_f_stat",
            "leave_one_out_pass",
            "strongest_iv_pass",
            "convergence_pass",
            "overlap_status",
        }
        self.assertTrue(required.issubset(row))
        self.assertEqual(row["candidate_id"], "candidate-1")
        self.assertEqual(row["overlap_status"], "not_assessed")
        self.assertTrue(row["instrument_pass"])

    def test_marks_candidate_as_instrument_failure_when_any_iv_has_f_below_ten(self):
        variants = pd.DataFrame(
            {
                "variant_id": ["rs1", "rs2", "rs3"],
                "f_stat": [9.99, 20.0, 30.0],
            }
        )

        row = diagnostic_row("candidate-weak", variants)

        self.assertFalse(row["instrument_pass"])

    def test_diagnostic_row_rejects_missing_f_stat(self):
        with self.assertRaises(ValueError):
            diagnostic_row("candidate-1", pd.DataFrame({"variant_id": ["rs1"]}))

    def test_adds_f_statistics_from_eqtl_effect_and_standard_error(self):
        variants = pd.DataFrame({"beta": [0.2, 0.3], "se": [0.1, 0.1]})

        result = with_f_statistics(variants)

        self.assertAlmostEqual(result["f_stat"].iloc[0], 4.0)
        self.assertAlmostEqual(result["f_stat"].iloc[1], 9.0)
