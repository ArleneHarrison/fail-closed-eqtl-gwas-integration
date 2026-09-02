import unittest

from audit_1000g_eur_reference import classify_record


class Audit1000GEurReferenceTests(unittest.TestCase):
    def test_retains_multiallelic_and_high_missing_records_with_reasons(self):
        fields = ["11", "100", ".", "A", "C,G", ".", "PASS", ".", "GT", "0|1", "./."]
        record = classify_record(fields, [0, 1], 0.05)
        self.assertFalse(record["include_for_ld"])
        self.assertIn("multiallelic", record["exclusion_reason"])
        self.assertIn("missingness_exceeds_threshold", record["exclusion_reason"])

    def test_accepts_biallelic_pass_record_without_missing_calls(self):
        fields = ["11", "100", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0|1:8", "1|1:7"]
        record = classify_record(fields, [0, 1], 0.05)
        self.assertTrue(record["include_for_ld"])
        self.assertEqual(record["eur_missing_n"], 0)
