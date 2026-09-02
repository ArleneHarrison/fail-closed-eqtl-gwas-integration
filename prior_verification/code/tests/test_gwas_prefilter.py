import unittest

from hcsmr.reanalysis.gwas_prefilter import filter_rows_for_regions, resolve_delimiter


class GwasPrefilterTests(unittest.TestCase):
    def test_retains_header_and_rows_in_any_requested_region(self):
        lines = [
            "chromosome\tbase_pair_location\tbeta\n",
            "11\t308015\t0.1\n",
            "11\t500000\t0.2\n",
            "11\t700000\t0.3\n",
        ]
        result = list(filter_rows_for_regions(lines, [("11", 300000, 400000), ("11", 650000, 750000)]))
        self.assertEqual(result, [lines[0], lines[1], lines[3]])

    def test_resolves_named_tab_delimiter(self):
        self.assertEqual(resolve_delimiter("tab"), "\t")
        self.assertEqual(resolve_delimiter(","), ",")
