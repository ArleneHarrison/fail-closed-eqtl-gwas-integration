from pathlib import Path

import numpy as np
import pandas as pd

from build_cross_ancestry_ld_benchmark import matrix_diagnostics, variant_key, write_matrix


def test_variant_key_normalizes_chromosome_and_alleles() -> None:
    assert variant_key("chr11", 123, "a", "g") == "chr11_123_A_G"


def test_matrix_diagnostics_separates_rank_from_positive_spectrum_condition() -> None:
    matrix = np.array([[1.0, 1.0], [1.0, 1.0]])
    result = matrix_diagnostics(matrix)
    assert result["rank"] == 1
    assert result["condition_number_positive_spectrum"] == 1.0


def test_write_matrix_preserves_identifier_order(tmp_path: Path) -> None:
    variants = ["chr11_1_A_G", "chr11_2_C_T"]
    matrix = np.array([[1.0, 0.25], [0.25, 1.0]])
    path = tmp_path / "ld.tsv"
    write_matrix(path, variants, matrix)
    observed = pd.read_csv(path, sep="\t")
    assert observed["variant_id"].tolist() == variants
    np.testing.assert_allclose(observed.iloc[:, 1:].to_numpy(), matrix)
