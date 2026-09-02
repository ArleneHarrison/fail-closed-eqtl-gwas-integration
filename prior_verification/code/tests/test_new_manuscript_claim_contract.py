import json
from pathlib import Path
import unittest


CODE_ROOT = Path(__file__).parents[1]
CONTRACT = CODE_ROOT / "manuscript_claim_contract_2026-09-01.json"


class NewManuscriptClaimContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract = json.loads(CONTRACT.read_text(encoding="utf-8"))

    def test_manuscript_retains_the_nonbiological_scope_boundary(self):
        boundaries = self.contract["scope_boundaries"]
        self.assertFalse(boundaries["biological_association"])
        self.assertFalse(boundaries["causal_inference"])
        self.assertFalse(boundaries["clinical_utility"])
        self.assertFalse(boundaries["comparative_superiority"])
        self.assertEqual(boundaries["real_data_role"], "workflow-behavior demonstration")

    def test_verified_reference_identifiers_are_present(self):
        dois = set(self.contract["verified_dois"])
        for doi in (
            "10.1038/s41588-021-00924-w",
            "10.1371/journal.pgen.1004383",
            "10.1371/journal.pgen.1009440",
            "10.1038/s41588-022-01233-6",
            "10.1038/nature15393",
            "10.1093/nar/gkae1070",
            "10.1093/bioinformatics/btab665",
            "10.1038/s41467-021-27438-7",
            "10.1371/journal.pgen.1010299",
        ):
            self.assertIn(doi, dois)

    def test_strengthened_validation_counts_are_locked(self):
        validation = self.contract["validation_counts"]
        self.assertEqual(validation["automated_tests"], 93)
        self.assertEqual(validation["fault_injection_expected_codes"], "2500/2500")
        self.assertEqual(validation["numerical_matrices"], 80)
        self.assertEqual(validation["full_rank_condition_numbers"], 64)
        self.assertEqual(validation["permutation_matrices"], 75)
        self.assertEqual(validation["downstream_replicate_rows"], 500)
        self.assertEqual(validation["ancestry_panels"], 5)
        self.assertEqual(validation["archived_cases"], 8)

    def test_ld_policy_documents_no_silent_repair(self):
        policy = self.contract["ld_policy"]
        self.assertEqual(policy["ill_conditioned_code"], "E_LD_ILL_CONDITIONED")
        self.assertEqual(policy["matrix_action"], "do not repair the matrix")

    def test_submission_metadata_remains_author_owned(self):
        status = self.contract["submission_status"]
        self.assertEqual(status["state"], "NOT READY FOR FINAL SUBMISSION")
        self.assertIn("corresponding-author confirmation", status["authorization"])
        self.assertIn("archival DOI", status["repository_release"])
        self.assertIn("open-source license", status["repository_release"])

