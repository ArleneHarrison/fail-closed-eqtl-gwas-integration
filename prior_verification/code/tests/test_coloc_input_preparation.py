import unittest

import numpy as np
import pandas as pd

from hcsmr.reanalysis.coloc_input_preparation import align_harmonized_to_ld


class ColocInputPreparationTests(unittest.TestCase):
    def test_reorders_summary_statistics_to_ld_variant_order(self):
        table = pd.DataFrame(
            {"target_variant": ["v2", "v1", "v3"], "beta": [2.0, 1.0, 3.0]}
        )
        ld = np.eye(2)

        aligned, matrix, dropped = align_harmonized_to_ld(table, ["v1", "v2"], ld)

        self.assertEqual(aligned["target_variant"].tolist(), ["v1", "v2"])
        self.assertEqual(matrix.shape, (2, 2))
        self.assertEqual(dropped, ["v3"])

    def test_selects_one_molecular_trait_before_deduplicating_variants(self):
        table = pd.DataFrame(
            {
                "molecular_trait_id": ["probe_a", "probe_b", "probe_a", "probe_b"],
                "target_variant": ["v1", "v1", "v2", "v2"],
                "beta": [1.0, 2.0, 3.0, 4.0],
            }
        )

        aligned, matrix, dropped = align_harmonized_to_ld(
            table, ["v1", "v2"], np.eye(2), molecular_trait_id="probe_b"
        )

        self.assertEqual(aligned["beta"].tolist(), [2.0, 4.0])
        self.assertEqual(matrix.shape, (2, 2))
        self.assertEqual(dropped, [])
