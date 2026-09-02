"""Aggregate the per-scenario JSONL output of run_simulation.py into the
metrics needed for Figure 3 of the manuscript:

  - Type-I error per scenario  (causal=False scenarios)
  - Power per scenario          (causal=True scenarios)
  - hcsMR vs IVW comparison
  - Bayes-factor calibration at threshold 10

Outputs both CSV summary tables and PDF panels (Fig 3a/3b/3c/3f).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def load_jsonl(path: Path):
    out = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return out


def compute_decisions(records, *, causal_state_idx: int = 3):
    """For each scenario × replicate, derive per-state decisions for hcsMR and IVW.

    Returns a long DataFrame: one row per (scenario × replicate × cell-state).
    """
    rows = []
    for r in records:
        sc = r["scenario"]
        theta_true = r["theta_true"]
        bf10 = r["bf10"]
        q025 = r["posterior_q025"]
        q975 = r["posterior_q975"]
        hcsmr_pos = r["hcsmr_pos"]
        ivw_pos = r["ivw_pos"]
        for c, (tt, bf, lo, hi, hp, ip) in enumerate(zip(theta_true, bf10, q025, q975, hcsmr_pos, ivw_pos)):
            rows.append({
                **sc,
                "cell_state": c,
                "is_causal_state": int(tt != 0),
                "theta_true": tt,
                "bf10": bf,
                "ci_excludes_zero": int((lo > 0) or (hi < 0)),
                "hcsmr_positive": int(hp),
                "ivw_positive": int(ip),
            })
    return pd.DataFrame(rows)


def summarise_by_scenario(df: pd.DataFrame):
    grp_cols = ["cell_state_freq", "pleiotropy_fraction", "ivs_per_state", "F_stat_mean", "is_causal"]
    out = []
    for keys, g in df.groupby(grp_cols, dropna=False):
        # Type-I error: among "is_causal_state==0" cell states, fraction that hcsMR classifies as positive
        null_mask = g["is_causal_state"] == 0
        alt_mask = g["is_causal_state"] == 1
        t1e_hcsmr = g.loc[null_mask, "hcsmr_positive"].mean() if null_mask.sum() else np.nan
        t1e_ivw = g.loc[null_mask, "ivw_positive"].mean() if null_mask.sum() else np.nan
        power_hcsmr = g.loc[alt_mask, "hcsmr_positive"].mean() if alt_mask.sum() else np.nan
        power_ivw = g.loc[alt_mask, "ivw_positive"].mean() if alt_mask.sum() else np.nan
        # BF₁₀ false-positive rate at BF₁₀>10 under null
        bf_null_pos = (g.loc[null_mask, "bf10"] > 10).mean() if null_mask.sum() else np.nan
        out.append({
            **dict(zip(grp_cols, keys)),
            "n_null_states": int(null_mask.sum()),
            "n_alt_states": int(alt_mask.sum()),
            "type1_error_hcsmr": t1e_hcsmr,
            "type1_error_ivw": t1e_ivw,
            "power_hcsmr": power_hcsmr,
            "power_ivw": power_ivw,
            "bf10_fpr": bf_null_pos,
        })
    return pd.DataFrame(out)


def make_plots(summary: pd.DataFrame, out_dir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Fig 3a: T1E heatmap ---
    null_df = summary[summary["is_causal"] == False]
    if len(null_df):
        pivot_t1e = null_df.pivot_table(
            index="cell_state_freq", columns="pleiotropy_fraction",
            values="type1_error_hcsmr", aggfunc="mean",
        )
        plt.figure(figsize=(5, 4))
        sns.heatmap(pivot_t1e, annot=True, fmt=".3f", cmap="RdYlGn_r", vmin=0, vmax=0.15,
                    cbar_kws={"label": "Type-I error"})
        plt.title("Fig 3a — hcsMR type-I error\n(nominal α = 0.05 reference)")
        plt.xlabel("Pleiotropy fraction")
        plt.ylabel("Cell-state frequency")
        plt.tight_layout()
        plt.savefig(out_dir / "fig3a_t1e_heatmap.pdf")
        plt.savefig(out_dir / "fig3a_t1e_heatmap.png", dpi=200)
        plt.close()

    # --- Fig 3b: Power vs cell-state frequency, hcsMR vs IVW ---
    alt_df = summary[summary["is_causal"] == True]
    if len(alt_df):
        plt.figure(figsize=(6, 4))
        for method in ["power_hcsmr", "power_ivw"]:
            curve = alt_df.groupby("cell_state_freq")[method].mean()
            plt.plot(curve.index, curve.values, marker="o", label=method.replace("power_", "").upper())
        plt.xscale("log")
        plt.xlabel("Cell-state frequency (log scale)")
        plt.ylabel("Power")
        plt.title("Fig 3b — Power vs cell-state frequency\n(hcsMR vs bulk-IVW)")
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_dir / "fig3b_power_curve.pdf")
        plt.savefig(out_dir / "fig3b_power_curve.png", dpi=200)
        plt.close()

    # --- Fig 3f: BF₁₀ FPR ---
    if len(null_df):
        pivot_bf = null_df.pivot_table(
            index="cell_state_freq", columns="pleiotropy_fraction",
            values="bf10_fpr", aggfunc="mean",
        )
        plt.figure(figsize=(5, 4))
        sns.heatmap(pivot_bf, annot=True, fmt=".3f", cmap="RdYlGn_r", vmin=0, vmax=0.15,
                    cbar_kws={"label": "Fraction with BF₁₀ > 10 under null"})
        plt.title("Fig 3f — Bayes-factor calibration\n(should be ≤ 0.05 under standard interpretation)")
        plt.xlabel("Pleiotropy fraction")
        plt.ylabel("Cell-state frequency")
        plt.tight_layout()
        plt.savefig(out_dir / "fig3f_bf_calibration.pdf")
        plt.savefig(out_dir / "fig3f_bf_calibration.png", dpi=200)
        plt.close()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in", dest="in_path", required=True, help="results.jsonl path")
    p.add_argument("--out", required=True, help="output dir")
    args = p.parse_args()

    in_path = Path(args.in_path)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    records = load_jsonl(in_path)
    print(f"loaded {len(records)} scenario records from {in_path}")
    if not records:
        return

    df = compute_decisions(records)
    df.to_csv(out_dir / "per_state_decisions.csv", index=False)
    print(f"wrote per_state_decisions.csv ({len(df)} rows)")

    summary = summarise_by_scenario(df)
    summary.to_csv(out_dir / "scenario_summary.csv", index=False)
    print(f"wrote scenario_summary.csv ({len(summary)} scenarios)")

    print("\n=== Scenario summary ===")
    print(summary.to_string())

    make_plots(summary, out_dir)
    print(f"\nplots written to {out_dir}")


if __name__ == "__main__":
    main()
