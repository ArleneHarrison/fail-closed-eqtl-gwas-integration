"""Cross-check LD numerical diagnostics with an R/SVD oracle and metamorphic tests."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from hcsmr.reanalysis.preflight_gate import validate_ld


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def correlation_from_spectrum(size: int, log10_span: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(size, size)))
    eigenvalues = np.logspace(0.0, -log10_span, size)
    covariance = q @ np.diag(eigenvalues) @ q.T
    scale = np.sqrt(np.diag(covariance))
    correlation = covariance / np.outer(scale, scale)
    return (correlation + correlation.T) / 2.0


def rank_deficient_correlation(size: int, latent_rank: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    loadings = rng.normal(size=(latent_rank, size))
    covariance = loadings.T @ loadings
    scale = np.sqrt(np.diag(covariance))
    correlation = covariance / np.outer(scale, scale)
    return (correlation + correlation.T) / 2.0


def equicorrelation(size: int, rho: float = 0.2) -> np.ndarray:
    """Construct a deterministic full-rank correlation matrix at realistic scale."""
    matrix = np.full((size, size), rho, dtype=float)
    np.fill_diagonal(matrix, 1.0)
    return matrix


def gate_record(case_id: str, matrix: np.ndarray, case_type: str, seed: int) -> dict:
    failures, diagnostics = validate_ld(
        matrix,
        max_condition_number=1e14,
        rank_policy="diagnostic",
    )
    return {
        "case_id": case_id,
        "case_type": case_type,
        "seed": seed,
        "n_variants": matrix.shape[0],
        "gate_codes": ";".join(failure.code for failure in failures) or "OK",
        **diagnostics,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--rscript", default=shutil.which("Rscript"))
    parser.add_argument(
        "--real-ld-npz",
        type=Path,
        help="Optional stored LD archive containing an 'ld' matrix for a construction-level cross-check.",
    )
    args = parser.parse_args()
    if not args.rscript:
        raise SystemExit("Rscript is required for the independent SVD oracle")

    out_dir = args.out_dir.resolve()
    matrix_dir = out_dir / "oracle_matrices"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    gate_rows: list[dict] = []
    permutation_rows: list[dict] = []
    manifest_rows: list[dict] = []

    generated: list[tuple[str, str, int, np.ndarray, bool]] = []
    for size in (8, 16, 32):
        for span in (0.0, 2.0, 6.0, 10.0):
            for seed in range(5):
                case_id = f"full_n{size}_span{span:g}_seed{seed}"
                generated.append((case_id, "full_rank", seed, correlation_from_spectrum(size, span, seed + size * 100), True))
    for size in (12, 24, 36):
        for seed in range(5):
            latent_rank = max(3, size // 3)
            case_id = f"deficient_n{size}_rank{latent_rank}_seed{seed}"
            generated.append((case_id, "rank_deficient", seed, rank_deficient_correlation(size, latent_rank, seed + size * 100), True))

    for size in (225, 620, 1000, 1776):
        case_id = f"large_equicorrelation_n{size}_rho0.2"
        generated.append((case_id, "large_full_rank", 0, equicorrelation(size), False))

    if args.real_ld_npz:
        real_archive = np.load(args.real_ld_npz.resolve(), allow_pickle=False)
        real_matrix = np.asarray(real_archive["ld"], dtype=float)
        generated.append(("stored_coordinate_ld", "stored_real_ld", 0, real_matrix, False))

    for case_id, case_type, seed, matrix, permutation_tested in generated:
        record = gate_record(case_id, matrix, case_type, seed)
        gate_rows.append(record)
        matrix_path = matrix_dir / f"{case_id}.tsv"
        np.savetxt(matrix_path, matrix, delimiter="\t", fmt="%.17g")
        manifest_rows.append({
            "case_id": case_id,
            "matrix_path": str(matrix_path),
            "rank_tolerance": record["rank_tolerance"],
        })
        if permutation_tested:
            rng = np.random.default_rng(seed + matrix.shape[0] * 1000)
            permutation = rng.permutation(matrix.shape[0])
            permuted = matrix[np.ix_(permutation, permutation)]
            permuted_record = gate_record(f"{case_id}_permuted", permuted, case_type, seed)
            original_condition = record["condition_number_2"]
            permuted_condition = permuted_record["condition_number_2"]
            relative_condition_difference = (
                None if original_condition is None else abs(permuted_condition - original_condition) / original_condition
            )
            permutation_rows.append({
                "case_id": case_id,
                "rank_equal": record["numerical_rank"] == permuted_record["numerical_rank"],
                "psd_equal": record["positive_semidefinite"] == permuted_record["positive_semidefinite"],
                "relative_condition_difference": relative_condition_difference,
                "permutation_invariant": (
                    record["numerical_rank"] == permuted_record["numerical_rank"]
                    and record["positive_semidefinite"] == permuted_record["positive_semidefinite"]
                    and (relative_condition_difference is None or relative_condition_difference <= 1e-6)
                ),
            })

    gate_frame = pd.DataFrame(gate_rows)
    manifest = pd.DataFrame(manifest_rows)
    permutation_frame = pd.DataFrame(permutation_rows)
    gate_path = out_dir / "gate_diagnostics.tsv"
    manifest_path = out_dir / "oracle_manifest.tsv"
    oracle_path = out_dir / "r_svd_oracle.tsv"
    gate_frame.to_csv(gate_path, sep="\t", index=False)
    manifest.to_csv(manifest_path, sep="\t", index=False)
    permutation_frame.to_csv(out_dir / "permutation_invariance.tsv", sep="\t", index=False)

    oracle_script = Path(__file__).with_name("run_svd_oracle.R")
    # R on some Windows installations cannot open Unicode project paths. Stage
    # only the oracle inputs in an ASCII temporary directory; keep the canonical
    # matrices and all reported outputs in the project tree.
    with tempfile.TemporaryDirectory(prefix="hcsmr_oracle_") as temporary:
        stage = Path(temporary)
        staged_manifest_rows = []
        for item in manifest_rows:
            source = Path(item["matrix_path"])
            staged = stage / source.name
            shutil.copyfile(source, staged)
            staged_manifest_rows.append({**item, "matrix_path": str(staged)})
        staged_manifest = stage / "oracle_manifest.tsv"
        staged_output = stage / "r_svd_oracle.tsv"
        pd.DataFrame(staged_manifest_rows).to_csv(staged_manifest, sep="\t", index=False)
        subprocess.run(
            [str(args.rscript), str(oracle_script), str(staged_manifest), str(staged_output)],
            check=True,
        )
        oracle = pd.read_csv(staged_output, sep="\t")
    oracle.to_csv(oracle_path, sep="\t", index=False)
    comparison = gate_frame.merge(oracle, on="case_id", validate="one_to_one")
    comparison["rank_agreement"] = comparison["numerical_rank"] == comparison["r_svd_rank"]
    comparison["condition_relative_error"] = np.where(
        comparison["condition_number_2"].notna(),
        np.abs(comparison["condition_number_2"] - comparison["r_svd_condition_number_2"])
        / comparison["r_svd_condition_number_2"],
        np.nan,
    )
    full_rank_mask = comparison["case_type"].isin(["full_rank", "large_full_rank"])
    comparison["condition_agreement"] = np.where(
        full_rank_mask,
        comparison["condition_relative_error"] <= 1e-6,
        comparison["condition_number_2"].isna(),
    )
    comparison.to_csv(out_dir / "oracle_comparison.tsv", sep="\t", index=False)

    indefinite_cases = []
    for rho in (-0.51, -0.60, -0.75):
        matrix = np.full((3, 3), rho)
        np.fill_diagonal(matrix, 1.0)
        failures, diagnostics = validate_ld(matrix, max_condition_number=1e14, rank_policy="diagnostic")
        codes = [failure.code for failure in failures]
        indefinite_cases.append({
            "rho": rho,
            "minimum_eigenvalue": diagnostics["min_eigenvalue"],
            "positive_semidefinite": diagnostics["positive_semidefinite"],
            "status_codes": ";".join(codes),
            "expected_code_present": "E_LD_NOT_PSD" in codes,
        })
    indefinite_frame = pd.DataFrame(indefinite_cases)
    indefinite_frame.to_csv(out_dir / "indefinite_matrix_rejection.tsv", sep="\t", index=False)

    full_rank = comparison.loc[full_rank_mask]
    summary = {
        "schema_version": "1.1",
        "scope": "independent R SVD oracle and permutation-based metamorphic validation of LD numerical diagnostics",
        "matrices_evaluated": int(len(comparison)),
        "generated_matrices": int((comparison["case_type"] != "stored_real_ld").sum()),
        "stored_real_ld_matrices": int((comparison["case_type"] == "stored_real_ld").sum()),
        "large_scale_generated_matrices": int((comparison["case_type"] == "large_full_rank").sum()),
        "maximum_matrix_dimension": int(comparison["n_variants"].max()),
        "full_rank_condition_agreements": int(full_rank["condition_agreement"].sum()),
        "full_rank_condition_cases": int(len(full_rank)),
        "rank_agreements": int(comparison["rank_agreement"].sum()),
        "rank_cases": int(len(comparison)),
        "permutation_invariance_agreements": int(permutation_frame["permutation_invariant"].sum()),
        "permutation_cases": int(len(permutation_frame)),
        "indefinite_matrix_rejections": int(indefinite_frame["expected_code_present"].sum()),
        "indefinite_matrix_cases": int(len(indefinite_frame)),
        "maximum_condition_relative_error": float(full_rank["condition_relative_error"].max()),
        "interpretation_boundary": "This cross-checks numerical implementation and invariance only; it is not external biological or downstream-model validation.",
    }
    (out_dir / "independent_validation_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    files = [path for path in out_dir.rglob("*") if path.is_file()]
    artifact_manifest = {
        "files": [
            {"path": str(path.relative_to(out_dir)), "byte_count": path.stat().st_size, "sha256": sha256_file(path)}
            for path in sorted(files) if path.name != "artifact_manifest.json"
        ]
    }
    (out_dir / "artifact_manifest.json").write_text(
        json.dumps(artifact_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
