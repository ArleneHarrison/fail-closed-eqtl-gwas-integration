"""Apply Bonferroni MTC + BH-FDR to the pilot results.

Adds adjusted thresholds:
  - bonf_alpha = 0.05 / n_tests (BF threshold equivalent ~ log(1/bonf))
  - p-value approx from BF via Wakefield approx; FDR via BH
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def bf_to_pval(bf, prior_odds=1.0):
    """Approximate two-sided p-value from a BF₁₀.

    Using the Sellke et al. bound: p_min >= 1 / (1 + e * BF * log(BF))
    For BF >> 1, p ~ e * BF * log(BF))^-1 (very approximate).
    Simpler: Bayes-factor → posterior odds → posterior probability,
    treat 1 - PP as a frequentist analogue p.
    """
    bf = np.asarray(bf, dtype=float)
    bf = np.where(np.isfinite(bf), bf, np.nan)
    post_odds = bf * prior_odds
    pp = post_odds / (1 + post_odds)
    # frequentist analogue
    p = 1 - pp
    # Avoid log(0)
    p = np.clip(p, 1e-300, 1.0)
    return p


def bh_fdr(p):
    """Standard BH 1995 FDR."""
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    q = p * n / ranks
    # enforce monotonicity
    q_sorted = q[order]
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    out = np.empty(n)
    out[order] = q_sorted
    return np.clip(out, 0, 1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--inputs", nargs="+", required=True,
                   help="Pairs of label,path  e.g. HF,paper/.../hf.csv AF,paper/.../af.csv")
    p.add_argument("--out", required=True)
    p.add_argument("--alpha", type=float, default=0.05)
    args = p.parse_args()

    dfs = []
    for spec in args.inputs:
        label, path = spec.split(",", 1)
        d = pd.read_csv(path)
        d["outcome"] = label
        d["p_approx"] = bf_to_pval(d["BF10"])
        d["ci_excludes_zero"] = ((d["post_q025"] > 0) | (d["post_q975"] < 0)).astype(int)
        dfs.append(d)
    all_df = pd.concat(dfs, ignore_index=True)
    n = len(all_df)

    # Bonferroni
    all_df["bonf_threshold"] = args.alpha / n
    all_df["passes_bonferroni"] = (all_df["p_approx"] < all_df["bonf_threshold"]).astype(int)
    # BH FDR
    all_df["q_bh"] = bh_fdr(all_df["p_approx"].fillna(1.0).values)
    all_df["passes_BH_q05"] = (all_df["q_bh"] < args.alpha).astype(int)
    all_df["is_robust_hit"] = (
        (all_df["BF10"] > 10) & (all_df["ci_excludes_zero"] == 1) & (all_df["n_iv"] >= 3)
    ).astype(int)

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    all_df.to_csv(out / "all_pilots_with_MTC.csv.gz", index=False, compression="gzip")

    print(f"Total tests: {n}")
    print(f"Bonferroni threshold (alpha={args.alpha}): p < {args.alpha/n:.3e}")

    # Summary
    summary = all_df.groupby("outcome").agg(
        n_tests=("BF10", "count"),
        n_robust_hits=("is_robust_hit", "sum"),
        n_bonf_pass=("passes_bonferroni", "sum"),
        n_BH_pass=("passes_BH_q05", "sum"),
        min_p=("p_approx", "min"),
        min_q=("q_bh", "min"),
    ).reset_index()
    summary.to_csv(out / "mtc_summary.csv", index=False)
    print("\n=== Per outcome (with MTC) ===")
    print(summary.to_string())

    # Top robust+Bonf hits
    top = all_df[(all_df["is_robust_hit"] == 1) & (all_df["passes_bonferroni"] == 1)].sort_values(
        ["outcome", "BF10"], ascending=[True, False]
    )
    top.to_csv(out / "robust_bonferroni_hits.csv", index=False)
    print(f"\n=== Robust + Bonferroni-pass hits ({len(top)} total) ===")
    if len(top):
        print(top[["outcome","gene_id","cell_state","parental_type","n_iv","post_mean","post_q025","post_q975","BF10","p_approx","q_bh"]].head(30).to_string())
    else:
        print("(none)")

    # Top robust+BH hits
    top_bh = all_df[(all_df["is_robust_hit"] == 1) & (all_df["passes_BH_q05"] == 1)].sort_values(
        ["outcome", "BF10"], ascending=[True, False]
    )
    top_bh.to_csv(out / "robust_BH_hits.csv", index=False)
    print(f"\n=== Robust + BH q<0.05 hits ({len(top_bh)} total) ===")
    if len(top_bh):
        print(top_bh[["outcome","gene_id","cell_state","parental_type","n_iv","post_mean","BF10","p_approx","q_bh"]].head(30).to_string())


if __name__ == "__main__":
    main()
