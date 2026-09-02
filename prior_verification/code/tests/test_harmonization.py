import unittest

import pandas as pd

from hcsmr.reanalysis.harmonization import (
    align_outcome_to_exposure,
    build_coordinate_coloc_audit,
    build_coordinate_coloc_table,
    join_by_lifted_position,
)


class HarmonizationTests(unittest.TestCase):
    def test_preserves_effect_when_outcome_effect_is_exposure_alt(self):
        result = align_outcome_to_exposure("C", "T", "T", "C", 0.2)
        self.assertEqual(result.status, "aligned")
        self.assertAlmostEqual(result.beta, 0.2)

    def test_flips_effect_when_outcome_effect_is_exposure_ref(self):
        result = align_outcome_to_exposure("C", "T", "C", "T", 0.2)
        self.assertEqual(result.status, "flipped")
        self.assertAlmostEqual(result.beta, -0.2)

    def test_rejects_ambiguous_alleles(self):
        result = align_outcome_to_exposure("A", "T", "A", "T", 0.2)
        self.assertEqual(result.status, "palindromic_ambiguous")
        self.assertIsNone(result.beta)

    def test_joins_eqtl_liftover_rows_to_coordinate_gwas(self):
        eqtl = pd.DataFrame({"variant": ["v1"], "ref": ["C"], "alt": ["T"]})
        lifted = pd.DataFrame(
            {"source_variant": ["v1"], "target_variant": ["chr4_100_C_T"], "status": ["mapped"]}
        )
        gwas = pd.DataFrame({"chromosome": [4], "base_pair_location": [100], "beta": [0.2]})

        result = join_by_lifted_position(eqtl, lifted, gwas, "chromosome", "base_pair_location")

        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["target_variant"], "chr4_100_C_T")

    def test_coordinate_coloc_table_has_required_effect_and_precision_fields(self):
        eqtl = pd.DataFrame(
            {"variant": ["v1"], "ref": ["C"], "alt": ["T"], "beta": [0.01], "se": [0.05]}
        )
        lifted = pd.DataFrame(
            {"source_variant": ["v1"], "target_variant": ["chr4_100_C_T"], "status": ["mapped"]}
        )
        gwas = pd.DataFrame(
            {
                "chromosome": [4], "base_pair_location": [100],
                "effect_allele": ["C"], "other_allele": ["T"],
                "beta": [0.2], "standard_error": [0.1], "n": [500],
            }
        )

        result = build_coordinate_coloc_table(
            eqtl, lifted, gwas,
            gwas_chrom_column="chromosome", gwas_position_column="base_pair_location",
            gwas_effect_column="effect_allele", gwas_other_column="other_allele",
            gwas_beta_column="beta", gwas_se_column="standard_error", gwas_n_column="n",
        )

        self.assertEqual(result.iloc[0]["alignment_status"], "flipped")
        self.assertAlmostEqual(result.iloc[0]["gwas_beta_aligned"], -0.2)
        self.assertEqual(result.iloc[0]["BP"], 100)
        self.assertAlmostEqual(result.iloc[0]["se_eqtl"], 0.05)
        self.assertAlmostEqual(result.iloc[0]["se_gwas"], 0.1)

    def test_audit_retains_unmapped_and_ambiguous_rows(self):
        eqtl = pd.DataFrame({
            "variant": ["v_mapped", "v_unmapped", "v_pal"],
            "ref": ["C", "A", "A"], "alt": ["T", "G", "T"],
            "beta": [0.01, 0.02, 0.03], "se": [0.05, 0.06, 0.07],
        })
        lifted = pd.DataFrame({
            "source_variant": ["v_mapped", "v_unmapped", "v_pal"],
            "target_variant": ["chr4_100_C_T", None, "chr4_101_A_T"],
            "status": ["mapped", "unmapped", "mapped"],
        })
        gwas = pd.DataFrame({
            "chromosome": [4, 4], "base_pair_location": [100, 101],
            "effect_allele": ["T", "A"], "other_allele": ["C", "T"],
            "beta": [0.2, 0.3], "standard_error": [0.1, 0.2], "n": [500, 500],
        })
        result = build_coordinate_coloc_audit(
            eqtl, lifted, gwas,
            gwas_chrom_column="chromosome", gwas_position_column="base_pair_location",
            gwas_effect_column="effect_allele", gwas_other_column="other_allele",
            gwas_beta_column="beta", gwas_se_column="standard_error", gwas_n_column="n",
        )
        self.assertEqual(set(result["alignment_status"]), {"aligned", "liftover_unmapped", "palindromic_ambiguous"})
