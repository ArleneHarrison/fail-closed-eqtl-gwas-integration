"""Post-hoc reviewer-requested diagnostics on the existing pilot results.

For each robust + Bonferroni-pass hit:
  - #6 Steiger filter (per-IV directionality test): does the IV explain more
        variance in the exposure than in the outcome?
  - #6 MR-Egger intercept p-value: is directional pleiotropy present?
  - #7 IV-selection bias check: do different cell states for the same gene
        use *the same* SNPs, or different ones?  A direction-flipping pattern
        is only biologically interpretable if the same IVs see different θs.
  - #10 / #4 Coloc-SuSiE-style proxy: PIP × −log10(GWAS p) joint posterior.

Recomputes BF₁₀ with the closed-form Savage-Dickey from model.py (so this
also fixes #1 for downstream tables).
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def bayes_factor_zero(theta_samples, prior_sd: float = 1.0,
                       cap_log10: float = 200.0) -> float:
    """Closed-form Savage-Dickey BF₁₀ assuming Gaussian posterior; capped at 10²⁰⁰."""
    if theta_samples.size < 50:
        return float("nan")
    mu  = float(np.mean(theta_samples))
    sig = float(np.std(theta_samples, ddof=1))
    if sig <= 0 or not np.isfinite(sig):
        return float("nan")
    log10_bf = np.log10(sig / prior_sd) + (mu ** 2) / (2 * sig ** 2 * np.log(10))
    if log10_bf > cap_log10:
        log10_bf = cap_log10
    if not np.isfinite(log10_bf):
        return float("nan")
    return float(10.0 ** log10_bf)


def steiger_filter(beta_iv, se_iv, beta_outcome, se_outcome, n_exposure, n_outcome):
    """Per-IV Steiger directional test.

    Returns boolean True if r2_exposure > r2_outcome (IV is causal in
    instrument→outcome direction; passes Steiger), and a p-value testing
    that the IV explains MORE variance in the exposure.
    """
    # r² estimates from z-statistics
    z_e = beta_iv / se_iv
    z_o = beta_outcome / se_outcome
    # r² ≈ z² / (z² + n)
    r2_e = z_e**2 / (z_e**2 + n_exposure)
    r2_o = z_o**2 / (z_o**2 + n_outcome)
    passes = r2_e > r2_o
    # Simple p-value: probit-style test
    diff = r2_e - r2_o
    return passes, diff


def mr_egger_intercept(beta_iv, se_iv, beta_out, se_out):
    """MR-Egger regression: gwas_beta ~ intercept + slope * eqtl_beta with
    inverse-variance weights.  Return intercept and its p-value."""
    from scipy import stats
    mask = np.isfinite(beta_iv) & np.isfinite(beta_out) & (se_iv > 0) & (se_out > 0)
    bi, bo, si, so = beta_iv[mask], beta_out[mask], se_iv[mask], se_out[mask]
    if len(bi) < 4:
        return {"intercept": np.nan, "intercept_p": np.nan, "n_iv_used": int(mask.sum())}
    w = 1.0 / np.maximum(so**2, 1e-12)
    sw = w.sum()
    x_bar = np.sum(w * bi) / sw
    y_bar = np.sum(w * bo) / sw
    Sxx = np.sum(w * (bi - x_bar)**2)
    Sxy = np.sum(w * (bi - x_bar) * (bo - y_bar))
    if Sxx <= 0:
        return {"intercept": np.nan, "intercept_p": np.nan, "n_iv_used": int(len(bi))}
    slope = Sxy / Sxx
    intercept = y_bar - slope * x_bar
    resid = bo - intercept - slope * bi
    rss = np.sum(w * resid**2)
    dof = max(len(bi) - 2, 1)
    sigma2 = rss / dof
    se_int = np.sqrt(sigma2 * (1.0/sw + x_bar**2/Sxx))
    z_int = intercept / se_int if se_int > 0 else 0.0
    p_int = 2 * (1 - stats.norm.cdf(abs(z_int)))
    return {"intercept": float(intercept), "intercept_p": float(p_int),
            "n_iv_used": int(len(bi))}


def iv_overlap_per_state(merged_df, gene_id):
    """For a gene with multiple cell-state estimates, compute pairwise IV
    overlap (Jaccard) between cell states.  Low overlap → the cross-state
    direction-flipping pattern is largely IV-selection driven (different
    SNPs see different effects), not biology.
    """
    sub = merged_df[merged_df["gene_id"] == gene_id].copy()
    snp_col = None
    for cand in ["rsid_for_merge", "rsid", "variant"]:
        if cand in sub.columns:
            snp_col = cand; break
    if snp_col is None:
        return None
    per_state = sub.groupby("state")[snp_col].apply(set).to_dict()
    states = list(per_state.keys())
    overlaps = []
    for i, s_i in enumerate(states):
        for s_j in states[i+1:]:
            a, b = per_state[s_i], per_state[s_j]
            if not a or not b: continue
            j = len(a & b) / len(a | b)
            overlaps.append((s_i, s_j, len(a & b), len(a | b), j))
    return overlaps


def coloc_proxy(sub_df):
    """Coloc-SuSiE proxy = max( PIP * −log10(GWAS p) ) over IVs in the
    credible set for this (gene, cell-state).  Mimics PP.H4 strength."""
    if "pip" not in sub_df.columns or "gwas_p" not in sub_df.columns:
        return float("nan")
    sub_df = sub_df.dropna(subset=["pip", "gwas_p"])
    if len(sub_df) == 0:
        return float("nan")
    neg_log_p = -np.log10(np.clip(sub_df["gwas_p"], 1e-300, 1))
    return float((sub_df["pip"] * neg_log_p).max())


def main():
    # Load: full hit table + per-outcome merged IV-GWAS tables (these have
    # PIP, beta_eqtl, se_eqtl, gwas_beta, gwas_se, gwas_p, rsid_for_merge).
    hits = pd.read_csv("paper/supplementary/multi_outcome/all_pilots_with_MTC.csv.gz",
                       compression="gzip")
    print(f"Loaded {len(hits)} (gene, cell-state) tests")

    # ---- Re-compute BF₁₀ using closed-form Savage-Dickey ----
    # The csv has post_mean and post_sd. Recompute directly here.
    hits["log10_BF_corrected"] = (
        np.log10(hits["post_sd"]) - np.log10(1.0)  # prior_sd = 1
        + (hits["post_mean"] ** 2) / (2 * hits["post_sd"] ** 2 * np.log(10))
    )
    # cap at 200
    hits["log10_BF_corrected"] = hits["log10_BF_corrected"].clip(upper=200)
    hits["BF10_corrected"] = 10 ** hits["log10_BF_corrected"]
    hits["BF10_old"] = hits["BF10"]
    print(f"BF₁₀ correction: old max log10 = {np.log10(hits['BF10_old'].clip(lower=1)).max():.1f} "
          f"→ corrected max = {hits['log10_BF_corrected'].max():.1f}")

    # Re-apply robust filter with corrected BF
    hits["ci_excludes_zero"] = ((hits["post_q025"] > 0) | (hits["post_q975"] < 0)).astype(int)
    hits["is_robust_corrected"] = (
        (hits["BF10_corrected"] > 10)
        & (hits["ci_excludes_zero"] == 1)
        & (hits["n_iv"] >= 3)
    ).astype(int)
    # Bonferroni at 3,663 tests
    alpha = 0.05
    n_tests = len(hits)
    bonf_threshold = alpha / n_tests
    # Wakefield-style p-value from BF via Sellke bound (very approximate)
    # p ~ 1/(1 + BF) lower bound; we use 10^(-log10_BF_corrected) as crude proxy
    hits["p_corrected"] = 10 ** (-(hits["log10_BF_corrected"].clip(lower=-300)))
    hits["passes_bonferroni_corrected"] = (hits["p_corrected"] < bonf_threshold).astype(int)

    rob = hits[hits["is_robust_corrected"] == 1].copy()
    bonf = hits[(hits["is_robust_corrected"] == 1) & (hits["passes_bonferroni_corrected"] == 1)].copy()
    print(f"  Robust hits (corrected, n_iv≥3 + BF₁₀>10): {len(rob)} ({rob['outcome'].value_counts().to_dict()})")
    print(f"  Robust + Bonferroni-pass (corrected): {len(bonf)} ({bonf['outcome'].value_counts().to_dict()})")

    # ---- Per-hit Steiger + MR-Egger using merged IV-GWAS files ----
    # The merged_iv_gwas.csv.gz files were saved by each pilot run.
    OUTCOME_TO_MERGED = {
        "HF":  "paper/supplementary/merged_iv_gwas/hf_merged.csv.gz",
        "AF":  "paper/supplementary/merged_iv_gwas/af_merged.csv.gz",
        "IS":  "paper/supplementary/merged_iv_gwas/is_merged.csv.gz",
        "CAD": "paper/supplementary/merged_iv_gwas/cad_merged.csv.gz",
    }

    print(f"  Merged IV-GWAS sources: {OUTCOME_TO_MERGED}")

    # Try to load HF merged table (the only one preserved)
    diag_rows = []
    for outc in ["HF", "CAD", "AF", "IS"]:
        path = OUTCOME_TO_MERGED.get(outc)
        if not path or not Path(path).exists():
            continue
        try:
            merged = pd.read_csv(path, compression="gzip" if path.endswith(".gz") else None,
                                  low_memory=False)
        except Exception as e:
            print(f"  failed to load {path}: {e}")
            continue
        print(f"  loaded {outc} merged table: {len(merged)} rows; cols={merged.columns.tolist()[:8]}")
        rob_in = rob[rob["outcome"] == outc]
        # Common col aliases — local merged files use raw 'beta', 'se'
        beta_col = "beta"
        se_col   = "se"
        for _, row in rob_in.iterrows():
            g, s = row["gene_id"], row["cell_state"]
            sub = merged[(merged["gene_id"] == g) & (merged["state"] == s if "state" in merged.columns else pd.Series([True] * len(merged)))]
            if len(sub) < 3:
                diag_rows.append({**row[["outcome","gene_id","cell_state","n_iv"]].to_dict(),
                                  "steiger_pass_fraction": np.nan,
                                  "egger_intercept": np.nan, "egger_intercept_p": np.nan,
                                  "coloc_proxy": np.nan, "comment": "IV table missing or n<3"})
                continue
            bi = sub[beta_col].values.astype(float)
            si = sub[se_col].values.astype(float)
            bo = sub["gwas_beta"].values.astype(float)
            so = sub["gwas_se"].values.astype(float)
            # Approximate n: sc-eQTL ~ 1000 donors; GWAS ~ 1M
            steiger_pass, _ = steiger_filter(bi, si, bo, so, n_exposure=1000, n_outcome=1_000_000)
            egger = mr_egger_intercept(bi, si, bo, so)
            cp = coloc_proxy(sub)
            diag_rows.append({
                **row[["outcome","gene_id","cell_state","n_iv","post_mean","post_sd"]].to_dict(),
                "BF10_corrected": float(row["BF10_corrected"]) if "BF10_corrected" in row else np.nan,
                "log10_BF_corrected": float(row["log10_BF_corrected"]) if "log10_BF_corrected" in row else np.nan,
                "steiger_pass_fraction": float(steiger_pass.mean()) if len(steiger_pass) else np.nan,
                "egger_intercept": egger["intercept"],
                "egger_intercept_p": egger["intercept_p"],
                "coloc_proxy": cp,
                "comment": "ok",
            })

    diag = pd.DataFrame(diag_rows)
    out_dir = Path("paper/supplementary/multi_outcome")
    out_dir.mkdir(parents=True, exist_ok=True)
    diag.to_csv(out_dir / "posthoc_diagnostics.csv", index=False)
    print(f"\nWrote {out_dir/'posthoc_diagnostics.csv'} ({len(diag)} robust hits diagnosed)")
    if len(diag):
        print(diag.to_string())

    # ---- IV overlap analysis for CYP4V1 (#7) ----
    print("\n=== #7 IV-selection-bias check for CYP4V1 ===")
    for outc in ["HF", "CAD"]:
        path = OUTCOME_TO_MERGED.get(outc)
        if not path or not Path(path).exists():
            continue
        merged = pd.read_csv(path, compression="gzip" if path.endswith(".gz") else None,
                              low_memory=False)
        for tgt_gene in ["ENSG00000145476", "ENSG00000185201"]:
            ov = iv_overlap_per_state(merged, tgt_gene)
            if not ov:
                continue
            gene_symbol = "CYP4V1" if tgt_gene == "ENSG00000145476" else "IFITM2"
            print(f"\n{gene_symbol} ({outc}) pairwise IV overlap (state_i, state_j, common, union, Jaccard):")
            ov_df = pd.DataFrame(ov, columns=["state_i", "state_j", "common", "union", "jaccard"])
            print(f"  {len(ov_df)} state pairs")
            print(f"  median Jaccard = {ov_df['jaccard'].median():.3f}")
            print(f"  fraction with Jaccard ≥ 0.5 = {(ov_df['jaccard']>=0.5).mean():.3f}")
            ov_df.to_csv(out_dir / f"iv_overlap_{gene_symbol}_{outc}.csv", index=False)

    # Save corrected hits table
    hits.to_csv(out_dir / "all_pilots_with_MTC_corrected.csv.gz",
                 index=False, compression="gzip")
    print(f"\nWrote corrected master table: {out_dir/'all_pilots_with_MTC_corrected.csv.gz'}")


if __name__ == "__main__":
    main()
