import unittest

from hcsmr.reanalysis.coloc_contract import validate_coloc_inputs


class ColocalizationContractTests(unittest.TestCase):
    def test_marks_region_ineligible_without_ancestry_matched_ld(self):
        region = {
            "region_id": "CYP4V1_HF",
            "eqtl_path": "eqtl.tsv.gz",
            "gwas_path": "gwas.tsv.gz",
            "ld_path": None,
            "ld_ancestry": "EUR",
        }

        result = validate_coloc_inputs(region)

        self.assertEqual(result.status, "ineligible")
        self.assertEqual(result.ineligibility_reason, "missing_ld_reference")

    def test_marks_region_ineligible_when_ld_file_is_absent(self):
        region = {
            "region_id": "CYP4V1_HF",
            "eqtl_path": "eqtl.tsv.gz",
            "gwas_path": "gwas.tsv.gz",
            "ld_path": "not_a_real_file.ld",
            "ld_ancestry": "EUR",
        }

        result = validate_coloc_inputs(region)

        self.assertEqual(result.status, "ineligible")
        self.assertEqual(result.ineligibility_reason, "ld_file_not_found")
