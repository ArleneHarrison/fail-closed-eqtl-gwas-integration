"""Portable checks of the review-corrected evidence; no raw data required."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
V = ROOT / "revision_validation"

def read(path):
    return pd.read_csv(path, sep="\t")

def test_multisignal_denominator_ledger():
    d = read(ROOT/"multisignal_revised/multisignal_replicates.tsv")
    assert len(d) == 1800
    groups = d.groupby(["sample_size", "architecture", "replicate"])
    assert len(groups) == 600
    assert (groups.size() == 3).all()
    assert (~d.completed).sum() == 48
    assert (d.completed & d.max_pp_h4_susie.isna()).sum() == 66
    assert d.max_pp_h4_susie.notna().sum() == 1686

def test_abf_independent_and_equal_across_ld_conditions():
    d = read(ROOT/"multisignal_revised/multisignal_replicates.tsv")
    assert d.abf_completed_independently.all()
    assert np.isfinite(d.pp_h4_abf).all()
    assert (d.groupby(["sample_size","architecture","replicate"]).pp_h4_abf.nunique() == 1).all()

def test_original_susie_outputs_preserved():
    old = read(ROOT/"multisignal_coloc/multisignal_replicates.tsv")
    new = read(ROOT/"multisignal_revised/multisignal_replicates.tsv")
    columns = [c for c in old if c != "pp_h4_abf"]
    pd.testing.assert_frame_equal(old[columns], new[columns])
    mask = old.pp_h4_abf.notna()
    assert np.max(abs(old.loc[mask,"pp_h4_abf"]-new.loc[mask,"pp_h4_abf"])) < 1e-12

def test_challenge_pairing_and_regions():
    d = read(V/"challenge_metadata.tsv")
    assert len(d) == 3520
    assert d.region.nunique() == 11
    assert d.base_id.nunique() == 440
    assert (d.groupby("base_id").size() == 8).all()
    assert (d.groupby(["region","fault"]).size() == 40).all()
    manifest = read(V/"challenge_region_manifest.tsv")
    assert manifest.min_position.min() > 5_000_000
    assert (manifest.variants == 100).all()

@pytest.mark.parametrize("model",["signed_graph","unsigned_graph","pooled_mlp","logistic","susie_kriging"])
def test_comparison_scores_and_operating_points(model):
    d = read(V/"all_challenge_predictions.tsv")
    summary = read(V/"all_model_challenge_summary.tsv").set_index("model").loc[model]
    p = d[model+"_anomaly_score"].to_numpy()
    y = (d.fault != "valid").to_numpy()
    assert np.isfinite(p).all()
    assert abs(roc_auc_score(y,p)-summary.anomaly_auc) < 1e-12
    assert int((p[~y] > summary.threshold).sum()) == summary.valid_false_positives
    assert summary.failed_scores == 0
    assert summary.region_bootstrap_low <= summary.anomaly_auc <= summary.region_bootstrap_high

def test_source64_content_and_identifier_gate_are_distinct():
    d = read(V/"source_float64_gate_challenge.tsv")
    d = d[d.representation == "source64"]
    content = d[d.fault != "explicit_ID_order"]
    structural = d[d.fault == "explicit_ID_order"]
    assert len(content) == 3520 and (content.status == "READY").all()
    assert len(structural) == 440 and (structural.status == "INELIGIBLE").all()
    assert structural.codes.str.contains("ORDER").all()

def test_float32_failures_retained():
    d = read(V/"source_float64_gate_challenge.tsv")
    d = d[(d.representation == "float32") & (d.fault != "explicit_ID_order")]
    assert (d.status == "READY").sum() == 640
    assert (d.status == "INELIGIBLE").sum() == 2880
    bad = d[d.codes.str.contains("E_LD_NOT_PSD")]
    assert bad.region.nunique() == 9

def test_outer_replicate_cluster_structure():
    d = read(ROOT/"ld_guardnet/ood_predictions.tsv")
    assert len(d) == 3200
    assert d.replicate.nunique() == 100
    assert (d.groupby("replicate").size() == 32).all()
    audit = json.loads((V/"audit.json").read_text())
    assert audit["bootstrap_repeats"] == 2000
    assert audit["frozen_checkpoint_max_abs_difference"] < 2e-6

def test_ablations_have_three_seeds_each():
    d = read(V/"ablation_training.tsv")
    assert len(d) == 6
    for _, group in d.groupby("mode"):
        assert set(group.seed) == {11,29,47}
    for name in ("signed_graph","pooled_mlp","unsigned_graph"):
        for seed in (11,29,47):
            p = ROOT/f"ld_guardnet/ld_guardnet_seed_{seed}.pt" if name == "signed_graph" else V/f"{name}_{seed}.pt"
            assert p.exists()
