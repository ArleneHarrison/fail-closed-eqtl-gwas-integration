import pandas as pd

from hcsmr.reanalysis.raw_gwas_verification import (
    candidate_query_keys,
    filter_raw_rows_for_candidates,
    summarize_candidate_qc,
    verify_candidate_rows,
)


def test_verifies_matching_candidate_against_raw_gwas_row():
    candidates = pd.DataFrame(
        {
            "candidate_id": ["CAD|ENSG00000185201|GTEx_Fibroblast"],
            "variant": ["chr11_308065_T_C"],
            "gwas_beta": [0.00619087],
            "gwas_se": [0.011217038],
            "gwas_p": [0.581],
        }
    )
    raw = pd.DataFrame(
        {
            "chromosome": [11],
            "base_pair_location": [308065],
            "effect_allele": ["C"],
            "other_allele": ["T"],
            "beta": [0.00619087],
            "standard_error": [0.011217038],
            "p_value": [0.581],
            "markername": ["11:308065_T_C"],
        }
    )

    result = verify_candidate_rows(candidates, raw)

    assert len(result) == 1
    assert result.loc[0, "match_status"] == "verified"
    assert bool(result.loc[0, "effect_allele_match"]) is True
    assert bool(result.loc[0, "beta_match"]) is True
    assert bool(result.loc[0, "se_match"]) is True
    assert bool(result.loc[0, "p_value_match"]) is True
    assert result.loc[0, "raw_effect_alleles"] == "C"
    assert result.loc[0, "raw_other_alleles"] == "T"
    assert result.loc[0, "allele_relation"] == "effect_matches_candidate_alt"


def test_records_missing_raw_row_without_implying_a_match():
    candidates = pd.DataFrame(
        {
            "candidate_id": ["CAD|ENSG00000185201|GTEx_Fibroblast"],
            "variant": ["chr11_308065_T_C"],
            "gwas_beta": [0.1],
            "gwas_se": [0.1],
            "gwas_p": [0.1],
        }
    )
    raw = pd.DataFrame(
        columns=[
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "p_value", "markername",
        ]
    )

    result = verify_candidate_rows(candidates, raw)

    assert result.loc[0, "match_status"] == "raw_row_not_found"
    assert bool(result.loc[0, "beta_match"]) is False


def test_records_alternative_raw_alleles_at_same_coordinate():
    candidates = pd.DataFrame(
        {
            "candidate_id": ["CAD|ENSG00000185201|GTEx_Fibroblast"],
            "variant": ["chr11_308065_T_C"],
            "gwas_beta": [0.1],
            "gwas_se": [0.1],
            "gwas_p": [0.1],
        }
    )
    raw = pd.DataFrame(
        {
            "chromosome": [11], "base_pair_location": [308065],
            "effect_allele": ["A"], "other_allele": ["T"],
            "beta": [0.2], "standard_error": [0.1], "p_value": [0.1],
            "markername": ["11:308065_T_A"],
        }
    )

    result = verify_candidate_rows(candidates, raw)

    assert result.loc[0, "match_status"] == "coordinate_found_marker_allele_mismatch"
    assert result.loc[0, "raw_markername"] == "11:308065_T_A"
    assert result.loc[0, "raw_effect_alleles"] == "A"


def test_candidate_gate_fails_when_any_raw_row_is_not_verified():
    qc = pd.DataFrame(
        {
            "candidate_id": ["CAD|gene|context", "CAD|gene|context"],
            "match_status": ["verified", "raw_row_not_found"],
        }
    )

    summary = summarize_candidate_qc(qc)

    assert bool(summary.loc[0, "raw_gwas_qc_pass"]) is False
    assert summary.loc[0, "raw_gwas_qc_status"] == "blocked_raw_row_not_found"


def test_raw_filter_keeps_same_coordinate_alternative_alleles():
    candidates = pd.DataFrame(
        {
            "candidate_id": ["CAD|gene|context"],
            "variant": ["chr11_308065_T_C"],
        }
    )
    raw = pd.DataFrame(
        {
            "chromosome": [11, 11],
            "base_pair_location": [308065, 308066],
            "markername": ["11:308065_T_A", "11:308066_T_A"],
        }
    )

    markers, coordinates = candidate_query_keys(candidates)
    kept = filter_raw_rows_for_candidates(raw, markers, coordinates)

    assert kept["markername"].tolist() == ["11:308065_T_A"]


def test_marks_reversed_effect_orientation_as_pending_reference_validation():
    candidates = pd.DataFrame(
        {
            "candidate_id": ["CAD|gene|context"],
            "variant": ["chr11_308065_T_C"],
            "gwas_beta": [0.1], "gwas_se": [0.1], "gwas_p": [0.1],
        }
    )
    raw = pd.DataFrame(
        {
            "chromosome": [11], "base_pair_location": [308065],
            "effect_allele": ["T"], "other_allele": ["C"],
            "beta": [0.1], "standard_error": [0.1], "p_value": [0.1],
            "markername": ["11:308065_C_T"],
        }
    )

    result = verify_candidate_rows(candidates, raw)

    assert result.loc[0, "allele_relation"] == "effect_matches_candidate_ref"
    assert result.loc[0, "required_resolution"] == "validate_build_strand_and_indel_normalization_before_effect_flip"
