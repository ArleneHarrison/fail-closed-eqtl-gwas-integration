"""Analyze hcsMR HF pilot results: hits, per-state summary, top genes per state,
gold-standard recapitulation, and Figure 4 panels.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

GOLD_HF = {
    # Cardiomyocyte/structural
    "TTN", "LMNA", "MYH7", "MYH6", "NPPA", "NPPB", "ACTC1", "BAG3",
    "FLNC", "RBM20", "PLN", "JPH2", "TBX5",
    # Lipid/cholesterol (CAD-relevant but HERMES HF picks up too)
    "PCSK9", "LDLR", "APOB", "APOE", "ANGPTL3", "SORT1", "LPA",
    # Inflammation
    "IL6R", "TREM2", "OLR1",
    # Conduction
    "GJA1", "SCN5A", "SCN10A", "PITX2", "ZFHX3", "KCNN3", "HCN4",
    # Adrenergic / RAAS
    "ADRB1", "ADRB2", "AGTR1", "GRK5",
    # Coagulation
    "F11", "F2", "F13B",
    # Endothelial / vascular
    "EDN1", "NOS3", "MIA3",
}


def load_results(results_path: Path) -> pd.DataFrame:
    df = pd.read_csv(results_path)
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--bf-threshold", type=float, default=10.0)
    p.add_argument("--gene-symbol-map", default=None,
                   help="Optional CSV with gene_id, gene_symbol columns")
    args = p.parse_args()

    df = load_results(args.results)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    print(f"Loaded {len(df)} (gene, cell_state) results across "
          f"{df['gene_id'].nunique()} genes × {df['cell_state'].nunique()} cell states")

    # Optionally map gene_id → symbol (if mapping provided)
    if args.gene_symbol_map:
        sym = pd.read_csv(args.gene_symbol_map)[["gene_id", "gene_symbol"]].drop_duplicates()
        df = df.merge(sym, on="gene_id", how="left")
    else:
        df["gene_symbol"] = df["gene_id"]

    # Filter hits
    df["ci_excludes_zero"] = ((df["post_q025"] > 0) | (df["post_q975"] < 0)).astype(int)
    df["is_hit"] = ((df["BF10"] > args.bf_threshold) & (df["ci_excludes_zero"] == 1)).astype(int)
    df["is_gold"] = df["gene_symbol"].isin(GOLD_HF).astype(int)

    n_hits = df["is_hit"].sum()
    print(f"\n=== HITS (BF₁₀ > {args.bf_threshold} and CI excludes 0): {n_hits} ===")

    if n_hits > 0:
        hits = df[df["is_hit"] == 1].sort_values("BF10", ascending=False)
        hits.to_csv(out / "pilot_hits.csv", index=False)
        print(hits.head(40).to_string())

    # Per-cell-state summary
    per_state = df.groupby(["cell_state", "parental_type"]).agg(
        n_genes_tested=("gene_id", "nunique"),
        n_hits=("is_hit", "sum"),
        mean_BF=("BF10", "mean"),
        median_BF=("BF10", "median"),
    ).reset_index().sort_values("n_hits", ascending=False)
    per_state.to_csv(out / "per_state_summary.csv", index=False)

    print("\n=== PER CELL STATE (top 15 by n_hits) ===")
    print(per_state.head(15).to_string())

    # Per-parental-type aggregation
    per_type = per_state.groupby("parental_type").agg(
        n_states=("cell_state", "count"),
        total_n_hits=("n_hits", "sum"),
        mean_mean_BF=("mean_BF", "mean"),
    ).reset_index().sort_values("total_n_hits", ascending=False)
    per_type.to_csv(out / "per_parental_type_summary.csv", index=False)

    print("\n=== PER PARENTAL TYPE ===")
    print(per_type.to_string())

    # Gold standard recapitulation
    if df["is_gold"].sum() > 0:
        gold_results = df[df["is_gold"] == 1].copy()
        gold_results.to_csv(out / "gold_standard_results.csv", index=False)
        gold_hit_rate = gold_results[gold_results["is_hit"] == 1]
        print(f"\n=== GOLD STANDARD GENES ===")
        print(f"  Tested: {gold_results['gene_symbol'].nunique()} of {len(GOLD_HF)}")
        print(f"  Hits: {gold_hit_rate['gene_symbol'].nunique()}")
        if len(gold_hit_rate):
            print(gold_hit_rate[["gene_symbol","cell_state","post_mean","BF10"]].to_string())

    # ----- Figures -----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        print("\n(matplotlib not available — skipping figures)")
        return

    # Figure 4a: BF10 distribution per cell state
    plt.figure(figsize=(12, 6))
    state_order = per_state["cell_state"].tolist()
    sns.boxplot(data=df, x="cell_state", y="BF10", order=state_order, showfliers=False)
    plt.yscale("log")
    plt.xticks(rotation=80, fontsize=7)
    plt.title(f"Fig 4a — Posterior BF₁₀ per cell state (HF / HERMES; n={df['gene_id'].nunique()} genes)")
    plt.axhline(args.bf_threshold, color="red", linestyle="--", label=f"BF₁₀ = {args.bf_threshold}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out / "fig4a_bf10_per_state.pdf")
    plt.savefig(out / "fig4a_bf10_per_state.png", dpi=200)
    plt.close()

    # Figure 4b: per-parental-type n_hits
    plt.figure(figsize=(8, 5))
    sns.barplot(data=per_type, x="parental_type", y="total_n_hits",
                order=per_type["parental_type"].tolist())
    plt.xticks(rotation=45, ha="right")
    plt.title(f"Fig 4b — Total hits per parental cell type (HF / HERMES)")
    plt.ylabel("# (gene × cell-state) hits")
    plt.tight_layout()
    plt.savefig(out / "fig4b_hits_per_parental_type.pdf")
    plt.savefig(out / "fig4b_hits_per_parental_type.png", dpi=200)
    plt.close()

    # Figure 4c: gene × cell-state heatmap of top genes
    if n_hits >= 5:
        top_hit_genes = df.loc[df["is_hit"] == 1, "gene_id"].value_counts().head(30).index.tolist()
        heat = df[df["gene_id"].isin(top_hit_genes)].pivot_table(
            index="gene_id", columns="cell_state", values="post_mean", aggfunc="mean"
        )
        plt.figure(figsize=(18, max(6, len(top_hit_genes) * 0.3)))
        sns.heatmap(heat, cmap="RdBu_r", center=0, cbar_kws={"label": "Posterior mean θ"})
        plt.xticks(rotation=80, fontsize=7)
        plt.yticks(fontsize=8)
        plt.title("Fig 4c — Posterior θ (gene × cell state) for top hit genes (HF)")
        plt.tight_layout()
        plt.savefig(out / "fig4c_gene_state_heatmap.pdf")
        plt.savefig(out / "fig4c_gene_state_heatmap.png", dpi=200)
        plt.close()

    print(f"\nFigures written to {out}")


if __name__ == "__main__":
    main()
