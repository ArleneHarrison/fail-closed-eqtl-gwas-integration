from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "complete_case_sensitivity", ROOT / "analyze_complete_case_ld_sensitivity.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)
PROTOCOL = ROOT / "config/protocol.v2.json"


def write_unit(root: Path, unit_id: str, primary: np.ndarray, complete: np.ndarray,
               ids: list[str], complete_samples: int, gate_status: str = "READY") -> None:
    unit = root / "prepared" / unit_id
    unit.mkdir(parents=True)
    np.savez_compressed(unit / "ld.npz", ld=primary, variant_ids=np.asarray(ids, dtype="U"))
    np.savez_compressed(unit / "ld_complete_case.npz", ld=complete, variant_ids=np.asarray(ids, dtype="U"))
    (unit / "ld_missingness_trace.json").write_text(json.dumps({
        "primary_reference_samples": 503,
        "complete_case_samples": complete_samples,
        "ordered_variant_sha256": MODULE.ordered_hash(ids),
    }), encoding="utf-8")
    gate_dir = root / "gates" / "units" / unit_id
    gate_dir.mkdir(parents=True)
    if ids:
        failures, diagnostics = MODULE.validate_ld(
            primary, max_condition_number=1e12, rank_policy="diagnostic"
        )
        codes = [failure.code for failure in failures]
    else:
        diagnostics = {"n_variants": 0, "numerical_rank": 0, "condition_number_2": None}
        codes = ["E_NO_ORDERED_OVERLAP"]
    (gate_dir / "gate.json").write_text(json.dumps({
        "status": gate_status if not codes else "INELIGIBLE",
        "status_codes": codes or ["OK"],
        "ld_diagnostics": diagnostics,
        "prepared_input_sha256": {
            "ld.npz": MODULE.sha256_file(unit / "ld.npz"),
        },
    }), encoding="utf-8")


class CompleteCaseSensitivityTests(unittest.TestCase):
    def make_inputs(self, root: Path) -> Path:
        primary = np.array([[1.0, 0.2, -0.1], [0.2, 1.0, 0.3], [-0.1, 0.3, 1.0]])
        complete = np.array([[1.0, 0.25, -0.05], [0.25, 1.0, 0.2], [-0.05, 0.2, 1.0]])
        write_unit(root, "u_matrix", primary, complete, ["v1", "v2", "v3"], 490)
        write_unit(root, "u_zero", np.empty((0, 0)), np.empty((0, 0)), [], 0)
        attempts = root / "attempts.tsv"
        pd.DataFrame([
            {"unit_id": "u_matrix", "region_id": "CAD_L01", "anchor_outcome": "CAD",
             "tissue_id": "QTD000131", "outcome_id": "CAD", "gene_id": "ENSG1"},
            {"unit_id": "u_zero", "region_id": "HF_L01", "anchor_outcome": "HF",
             "tissue_id": "QTD000136", "outcome_id": "HF", "gene_id": ""},
        ]).to_csv(attempts, sep="\t", index=False)
        return attempts

    def test_metrics_zero_overlap_and_determinism(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            attempts = self.make_inputs(root)
            first, second = root / "out1", root / "out2"
            MODULE.analyze(attempts, root / "prepared", root / "gates", PROTOCOL, first, 2)
            MODULE.analyze(attempts, root / "prepared", root / "gates", PROTOCOL, second, 2)
            self.assertEqual(
                (first / "complete_case_ld_sensitivity.tsv").read_bytes(),
                (second / "complete_case_ld_sensitivity.tsv").read_bytes(),
            )
            result = pd.read_csv(first / "complete_case_ld_sensitivity.tsv", sep="\t", dtype=str, keep_default_na=False)
            matrix = result.loc[result.unit_id == "u_matrix"].iloc[0]
            self.assertEqual(matrix.complete_case_ld_integrity_status, "PASS")
            self.assertEqual(matrix.primary_vs_complete_case_ld_status_concordant, "TRUE")
            self.assertEqual(matrix.primary_vs_complete_case_rank_concordant, "TRUE")
            delta = np.array([[0, .05, .05], [.05, 0, -.1], [.05, -.1, 0]])
            expected_relative = np.linalg.norm(delta, "fro") / np.linalg.norm(
                np.array([[1.0, .2, -.1], [.2, 1.0, .3], [-.1, .3, 1.0]]), "fro"
            )
            self.assertAlmostEqual(float(matrix.relative_frobenius_distance), expected_relative)
            self.assertAlmostEqual(float(matrix.max_absolute_off_diagonal_delta_r), 0.1)
            self.assertAlmostEqual(float(matrix.median_absolute_off_diagonal_delta_r), 0.05)
            self.assertAlmostEqual(float(matrix.p95_absolute_off_diagonal_delta_r), 0.095)
            zero = result.loc[result.unit_id == "u_zero"].iloc[0]
            self.assertEqual(zero.analysis_status, "N_A_ZERO_OVERLAP")
            self.assertEqual(zero.relative_frobenius_distance, "NA")
            self.assertEqual(zero.primary_vs_complete_case_rank_concordant, "NA")
            summary = json.loads((first / "run_summary.json").read_text(encoding="utf-8"))
            self.assertFalse(summary["primary_gate_or_dispatch_mutation"])
            self.assertFalse(summary["difference_threshold_selected"])

    def test_nonfinite_complete_case_is_reported_not_repaired(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            primary = np.eye(2)
            complete = np.full((2, 2), np.nan)
            write_unit(root, "u", primary, complete, ["v1", "v2"], 1)
            attempt = {"unit_id": "u"}
            row = MODULE.analyze_unit(attempt, root / "prepared", root / "gates", 1e12)
            self.assertEqual(row["complete_case_ld_integrity_status"], "FAIL")
            self.assertEqual(row["complete_case_ld_failure_codes"], "E_LD_NONFINITE")
            self.assertEqual(row["relative_frobenius_distance"], "NA")

    def test_variant_order_mismatch_stops_analysis(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            write_unit(root, "u", np.eye(2), np.eye(2), ["v1", "v2"], 500)
            np.savez_compressed(
                root / "prepared/u/ld_complete_case.npz", ld=np.eye(2),
                variant_ids=np.asarray(["v2", "v1"]),
            )
            with self.assertRaisesRegex(RuntimeError, "E_COMPLETE_CASE_VARIANT_ORDER_MISMATCH"):
                MODULE.analyze_unit({"unit_id": "u"}, root / "prepared", root / "gates", 1e12)

    def test_primary_gate_input_binding_mismatch_stops_analysis(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            write_unit(root, "u", np.eye(2), np.eye(2), ["v1", "v2"], 500)
            np.savez_compressed(
                root / "prepared/u/ld.npz", ld=np.array([[1.0, 0.1], [0.1, 1.0]]),
                variant_ids=np.asarray(["v1", "v2"]),
            )
            with self.assertRaisesRegex(RuntimeError, "E_PRIMARY_GATE_INPUT_BINDING_MISMATCH"):
                MODULE.analyze_unit({"unit_id": "u"}, root / "prepared", root / "gates", 1e12)


if __name__ == "__main__":
    unittest.main()
