import json
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

from hcsmr.reanalysis.harmonization import build_coordinate_coloc_audit
from hcsmr.reanalysis.preflight_gate import run_preflight


FIXTURE_PATH = Path(__file__).parent / "fixtures" / "fail_closed_synthetic_cases.json"
PROVENANCE = [{
    "label": "synthetic", "url": "https://example.org/synthetic", "version_or_build": "v1",
    "license_or_access": "CC0", "retrieval_date": "2026-08-26", "sha256": "0" * 64,
    "byte_count": 1, "structural_check": "passed",
}]


class FailClosedIntegrationFixtureTests(unittest.TestCase):
    def test_fixture_contract_lists_all_required_failure_modes(self):
        cases = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))["cases"]
        self.assertEqual({case["case_id"] for case in cases}, {
            "unmapped", "one_to_many", "strand", "indel", "duplicate", "vcf_absent", "ill_conditioned_ld"
        })

    def test_coordinate_pipeline_preserves_mapping_and_representation_failures(self):
        eqtl = pd.DataFrame({
            "variant": ["unmapped", "one_to_many", "strand", "indel"],
            "ref": ["A", "A", "C", "AT"], "alt": ["C", "C", "T", "A"],
            "beta": [0.1] * 4, "se": [0.1] * 4,
        })
        lifted = pd.DataFrame({
            "source_variant": ["unmapped", "one_to_many", "strand", "indel"],
            "target_variant": [None, None, "chr1_3_C_T", "chr1_4_AT_A"],
            "status": ["unmapped", "one_to_many", "mapped", "mapped"],
        })
        gwas = pd.DataFrame({
            "chromosome": [1, 1], "base_pair_location": [3, 4],
            "effect_allele": ["A", "AT"], "other_allele": ["G", "A"],
            "beta": [0.2, 0.2], "standard_error": [0.1, 0.1], "n": [100, 100],
        })
        audit = build_coordinate_coloc_audit(
            eqtl, lifted, gwas,
            gwas_chrom_column="chromosome", gwas_position_column="base_pair_location",
            gwas_effect_column="effect_allele", gwas_other_column="other_allele",
            gwas_beta_column="beta", gwas_se_column="standard_error", gwas_n_column="n",
        )
        self.assertEqual(set(audit["alignment_status"]), {
            "liftover_unmapped", "liftover_one_to_many", "strand_unresolved", "indel_unresolved"
        })

    def test_gate_retains_duplicate_absent_vcf_and_ill_conditioned_ld_as_distinct_codes(self):
        row = {
            "target_variant": "chr1_1_A_C", "eqtl_se": "0.1", "gwas_se": "0.1", "maf": "0.2",
            "an": "100", "gwas_n": "100", "alignment_status": "aligned", "vcf_status": "position_absent",
        }
        result = run_preflight(
            rows=[row, dict(row)],
            ld=np.array([[1.0, 1.0 - 1e-14], [1.0 - 1e-14, 1.0]]),
            provenance=PROVENANCE, eqtl_trait_type="quant", gwas_trait_type="cc",
            case_fraction=0.2, max_condition_number=1e12, ld_variant_ids=["chr1_1_A_C", "chr1_1_A_C"],
        )
        self.assertEqual(result["status"], "INELIGIBLE")
        self.assertTrue({"E_VARIANT_DUPLICATE", "E_VCF_ABSENT", "E_LD_ILL_CONDITIONED"}.issubset(result["status_codes"]))
