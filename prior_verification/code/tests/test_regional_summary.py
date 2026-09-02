import unittest

import pandas as pd

from hcsmr.reanalysis.regional_summary import filter_eqtl_region, filter_gwas_region


class RegionalSummaryTests(unittest.TestCase):
    def test_filters_eqtl_by_gene_and_inclusive_coordinate_window(self):
        frame = pd.DataFrame(
            {
                "gene_id": ["ENSG1", "ENSG1", "ENSG2"],
                "chromosome": [4, 4, 4],
                "position": [100, 201, 100],
            }
        )
        result = filter_eqtl_region(frame, "ENSG1", "4", 100, 200)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["position"], 100)

    def test_filters_gwas_hg19_coordinates(self):
        frame = pd.DataFrame({"CHR": [4, 4, 5], "BP": [99, 100, 100]})
        result = filter_gwas_region(frame, "4", 100, 100, "CHR", "BP")
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["BP"], 100)
