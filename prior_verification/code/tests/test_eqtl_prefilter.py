import unittest

from hcsmr.reanalysis.eqtl_prefilter import filter_rows_for_gene


class EqtlPrefilterTests(unittest.TestCase):
    def test_retains_header_and_exact_gene_id_rows(self):
        lines = [
            "variant\tgene_id\tbeta\n",
            "v1\tENSG1\t0.1\n",
            "v2\tENSG10\t0.2\n",
            "v3\tENSG1\t0.3\n",
        ]
        result = list(filter_rows_for_gene(lines, "ENSG1"))
        self.assertEqual(result, [lines[0], lines[1], lines[3]])
