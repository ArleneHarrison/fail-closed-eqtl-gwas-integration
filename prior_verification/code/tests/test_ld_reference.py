import unittest

import numpy as np

from hcsmr.reanalysis.ld_reference import (
    build_variant_key,
    compute_ld,
    parse_region,
    select_requested_variants,
)


class LdReferenceTests(unittest.TestCase):
    def test_parses_hg19_region(self):
        self.assertEqual(parse_region("chr4:185191567-187191567"), ("4", 185191567, 187191567))

    def test_builds_allele_aware_variant_key(self):
        self.assertEqual(build_variant_key("4", 186202630, "T", "C"), "chr4_186202630_T_C")

    def test_computes_symmetric_ld_and_excludes_missing_genotypes(self):
        dosages = np.array(
            [
                [0.0, 1.0, 2.0, np.nan],
                [0.0, 1.0, 2.0, 1.0],
                [2.0, 1.0, 0.0, 1.0],
            ]
        )

        ld = compute_ld(dosages)

        self.assertTrue(np.allclose(np.diag(ld), 1.0))
        self.assertTrue(np.allclose(ld, ld.T))
        self.assertAlmostEqual(ld[0, 2], -1.0)

    def test_selects_only_requested_allele_aware_variants(self):
        selected, missing = select_requested_variants(
            ["chr4_10_A_G", "chr4_11_C_T"], {"chr4_11_C_T", "chr4_12_G_A"}
        )

        self.assertEqual(selected, ["chr4_11_C_T"])
        self.assertEqual(missing, ["chr4_12_G_A"])
