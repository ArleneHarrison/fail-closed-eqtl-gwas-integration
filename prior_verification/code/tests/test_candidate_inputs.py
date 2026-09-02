import unittest

import pandas as pd

from hcsmr.reanalysis.candidate_inputs import build_candidate_variants


class CandidateInputTests(unittest.TestCase):
    def test_retains_empty_candidate_with_explicit_status(self):
        hits = pd.DataFrame(
            {
                "outcome": ["HF"],
                "gene_id": ["ENSG1"],
                "cell_state": ["state_a"],
                "is_robust_corrected": [1],
                "passes_bonferroni_corrected": [1],
            }
        )
        merged = {"HF": pd.DataFrame(columns=["gene_id", "state", "beta", "se"])}

        variants, status = build_candidate_variants(hits, merged)

        self.assertTrue(variants.empty)
        self.assertEqual(status.iloc[0]["status"], "missing_variant_rows")

    def test_adds_candidate_identifier_and_f_stat(self):
        hits = pd.DataFrame(
            {
                "outcome": ["HF"],
                "gene_id": ["ENSG1"],
                "cell_state": ["state_a"],
                "is_robust_corrected": [1],
                "passes_bonferroni_corrected": [1],
            }
        )
        merged = {
            "HF": pd.DataFrame(
                {"gene_id": ["ENSG1"], "state": ["state_a"], "beta": [0.2], "se": [0.1]}
            )
        }

        variants, status = build_candidate_variants(hits, merged)

        self.assertEqual(variants.iloc[0]["candidate_id"], "HF|ENSG1|state_a")
        self.assertAlmostEqual(variants.iloc[0]["f_stat"], 4.0)
        self.assertEqual(status.iloc[0]["status"], "ready_for_diagnostics")
