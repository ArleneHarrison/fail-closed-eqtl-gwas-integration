import unittest

from hcsmr.reanalysis.tabix_bgzf import reg2bins


class TabixBinningTests(unittest.TestCase):
    def test_bins_include_root_and_are_deterministic(self):
        bins = reg2bins(1_999_999, 2_100_000)
        self.assertEqual(bins[0], 0)
        self.assertEqual(bins, reg2bins(1_999_999, 2_100_000))
        self.assertGreater(len(bins), 1)

    def test_rejects_invalid_interval(self):
        with self.assertRaises(ValueError):
            reg2bins(5, 5)
