"""Aggregate hcsMR pilot results across multiple CVD outcomes for the cross-outcome
Figure 4 of the manuscript.

For each outcome (HF/AF/IS/CAD), load pilot_results.csv, filter hits, and
combine into:
  - hits_all_outcomes.csv: every hit across outcomes
  - per_outcome_summary.csv: hit counts and top hit per outcome
  - cross_outcome_heatmap.png: gene × outcome × cell-state matrix
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load_outcome(label, path):
    p = Path(path)
    if not p.exists():
        print(f"missing {p}")
        return None
    df = pd.read_csv(p)
    df["outcome"] = label
    df["ci_excludes_zero"] = ((df["post_q025"] > 0) | (df["post_q975"] < 0)).astype(int)
    df["is_hit"] = ((df["BF10"] > 10) & (df["ci_excludes_zero"] == 1)).astype(int)
    # Robust filter: require n_iv >= 3 to mitigate small-IV overfit
    df["is_hit_robust"] = (df["is_hit"] & (df["n_iv"] >= 3)).astype(int)
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--hf", default="paper/supplementary/pilot_partial/pilot_results_final.csv")
    p.add_argument("--af", default=None)
    p.add_argument("--is", dest="is_path", default=None)
    p.add_argument("--cad", default=None)
    p.add_argument("--out", default="paper/supplementary/multi_outcome")
    args = p.parse_args()

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    dfs = []
    for label, path in (("HF", args.hf), ("AF", args.af), ("IS", args.is_path), ("CAD", args.cad)):
        if path:
            d = load_outcome(label, path)
            if d is not None:
                print(f"  {label}: {len(d)} (gene, cell-state) fits, {d['is_hit'].sum()} loose hits, {d['is_hit_robust'].sum()} robust hits")
                dfs.append(d)

    if not dfs:
        print("no outcomes loaded")
        return

    all_df = pd.concat(dfs, ignore_index=True)
    all_df.to_csv(out / "all_outcomes.csv.gz", index=False, compression="gzip")

    # Hits-only table
    hits = all_df[all_df["is_hit_robust"] == 1].sort_values(["outcome", "BF10"], ascending=[True, False])
    hits.to_csv(out / "all_hits_robust.csv", index=False)
    print(f"\nTotal robust hits across outcomes: {len(hits)}")

    # Per outcome
    per_outcome = all_df.groupby("outcome").agg(
        n_genes=("gene_id", "nunique"),
        n_cell_state_fits=("gene_id", "count"),
        n_hits_loose=("is_hit", "sum"),
        n_hits_robust=("is_hit_robust", "sum"),
        max_BF10=("BF10", "max"),
        median_BF10=("BF10", "median"),
    ).reset_index()
    per_outcome.to_csv(out / "per_outcome_summary.csv", index=False)
    print("\n=== Per outcome ===")
    print(per_outcome.to_string())

    # Per outcome × parental
    per_op = all_df.groupby(["outcome", "parental_type"]).agg(
        n_hits=("is_hit_robust", "sum"),
        n_tests=("gene_id", "count"),
    ).reset_index()
    per_op_piv = per_op.pivot(index="parental_type", columns="outcome", values="n_hits").fillna(0).astype(int)
    per_op_piv.to_csv(out / "hits_by_outcome_x_parental.csv")
    print("\n=== Hits by outcome × parental type ===")
    print(per_op_piv.to_string())

    # Cross-outcome top hits
    print("\n=== Top 30 hits by BF10 across all outcomes (robust only) ===")
    print(hits[["outcome","gene_id","cell_state","parental_type","n_iv","post_mean","post_q025","post_q975","BF10"]].head(30).to_string())

    # Genes hit in multiple outcomes
    gene_outcome_hits = hits.groupby("gene_id")["outcome"].nunique()
    multi = gene_outcome_hits[gene_outcome_hits >= 2].index.tolist()
    if multi:
        print(f"\n=== Genes hit in MULTIPLE outcomes ({len(multi)} total) ===")
        multi_hits = hits[hits["gene_id"].isin(multi)]
        print(multi_hits[["outcome","gene_id","cell_state","parental_type","post_mean","BF10"]].to_string())

    # --- Figures ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        return

    plt.figure(figsize=(8, 5))
    sns.heatmap(per_op_piv, annot=True, fmt="d", cmap="YlGnBu",
                cbar_kws={"label": "# (gene × cell-state) hits"})
    plt.title("Fig 4 (cross-outcome) — Robust hits by parental cell type × outcome")
    plt.tight_layout()
    plt.savefig(out / "fig4_cross_outcome_heatmap.pdf")
    plt.savefig(out / "fig4_cross_outcome_heatmap.png", dpi=200)
    plt.close()
    print(f"\nFigure saved to {out}")


if __name__ == "__main__":
    main()
