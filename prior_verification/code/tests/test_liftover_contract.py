import unittest

from hcsmr.reanalysis.liftover_contract import complement_allele, transformed_alleles


class LiftoverContractTests(unittest.TestCase):
    def test_complements_dna_alleles(self):
        self.assertEqual(complement_allele("AGT"), "TCA")

    def test_flips_effect_alleles_on_negative_strand(self):
        self.assertEqual(transformed_alleles("A", "G", "+"), ("A", "G"))
        self.assertEqual(transformed_alleles("A", "G", "-"), ("T", "C"))
