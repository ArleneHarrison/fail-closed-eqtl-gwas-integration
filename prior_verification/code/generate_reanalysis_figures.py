"""Create data-backed supplementary figures for the conservative revision."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


METHOD_COLORS = {
    "ivw": "#0072B2",
    "weighted_median": "#009E73",
    "egger": "#D55E00",
}
METHOD_LABELS = {
    "ivw": "IVW",
    "weighted_median": "Weighted median",
    "egger": "MR-Egger",
}


def save_figure(figure: plt.Figure, path: Path) -> None:
    """Save both vector and high-resolution raster versions of a figure."""
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def make_calibration_figure(performance_path: Path, out_dir: Path) -> None:
    data = pd.read_csv(performance_path)
    null_data = data.loc[
        (data["true_theta"] == 0)
        & (data["n_instruments"].isin([3, 4, 5]))
        & data["type1_rate"].notna()
    ].copy()
    if null_data.empty:
        raise ValueError("No eligible sparse-instrument null scenarios were found.")

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.0), constrained_layout=True)
    for method in ("ivw", "weighted_median", "egger"):
        subset = null_data.loc[null_data["method"] == method]
        if subset.empty:
            continue
        for n_instruments, group in subset.groupby("n_instruments", sort=True):
            offsets = np.linspace(-0.20, 0.20, len(group))
            axes[0].scatter(
                np.full(len(group), n_instruments) + offsets,
                group["type1_rate"],
                color=METHOD_COLORS[method],
                marker={"ivw": "o", "weighted_median": "s", "egger": "^"}[method],
                alpha=0.78,
                s=26,
                linewidths=0.25,
                edgecolors="white",
                label=METHOD_LABELS[method]
                if n_instruments == subset["n_instruments"].min()
                else None,
            )
    axes[0].axhline(0.05, color="#222222", linestyle="--", linewidth=1.0, label="Nominal 0.05")
    axes[0].set(
        title="A. Null calibration in sparse-instrument scenarios",
        xlabel="Number of instruments",
        ylabel="Scenario-specific type-I error",
        xticks=[3, 4, 5],
        ylim=(0, 0.62),
    )
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    axes[0].grid(axis="y", alpha=0.22)

    summary = (
        null_data.groupby("method", as_index=False)
        .agg(mean_type1=("type1_rate", "mean"), max_type1=("type1_rate", "max"), scenarios=("type1_rate", "size"))
        .sort_values("mean_type1", ascending=True)
    )
    y_positions = np.arange(len(summary))
    axes[1].barh(
        y_positions,
        summary["mean_type1"],
        color=[METHOD_COLORS[m] for m in summary["method"]],
        alpha=0.88,
    )
    axes[1].axvline(0.05, color="#222222", linestyle="--", linewidth=1.0)
    for y, row in zip(y_positions, summary.to_dict("records")):
        axes[1].text(
            row["mean_type1"] + 0.006,
            y,
            f"mean {row['mean_type1']:.3f}; max {row['max_type1']:.3f}",
            va="center",
            fontsize=8,
        )
    axes[1].set(
        title="B. Aggregate null error across 3-5 instruments",
        xlabel="Mean scenario-specific type-I error",
        xlim=(0, 0.38),
        yticks=y_positions,
        yticklabels=[METHOD_LABELS[m] for m in summary["method"]],
    )
    axes[1].grid(axis="x", alpha=0.22)
    fig.suptitle("Baseline MR calibration audit (1,000 replicates per scenario-condition)", fontsize=12, y=1.03)
    save_figure(fig, out_dir / "figure_s1_baseline_calibration")


def read_coloc_pairs(root: Path, name: str, label: str) -> pd.DataFrame:
    data = pd.read_csv(root / name / "summary.csv")
    return data.assign(record=label)


def make_coloc_figure(coloc_root: Path, validation_path: Path, out_dir: Path) -> None:
    pairs = pd.concat(
        [
            read_coloc_pairs(coloc_root, "IFITM2_GTEx_Fibroblast_CAD", "IFITM2 / GTEx fibroblast / CAD"),
            read_coloc_pairs(coloc_root, "IFITM2_OneK1K_NK_CAD", "IFITM2 / OneK1K NK / CAD"),
            read_coloc_pairs(coloc_root, "TMEM80_Fairfax_monocyte_IFN24_CAD", "TMEM80 / Fairfax IFN24 / CAD"),
        ],
        ignore_index=True,
    )
    validation = pd.read_csv(validation_path)
    labels = [
        "CYP4V1 / Fairfax LPS24 / HF",
        "CYP4V1 / GTEx fibroblast / HF",
        "CYP4V1 / GTEx LV / HF",
        "CYP4V1 / GTEx whole blood / HF",
        "IFITM2 / GTEx fibroblast / CAD",
        "IFITM2 / OneK1K NK / CAD",
        "TMEM80 / Fairfax IFN24 / CAD",
        "TMEM80 / Fairfax naive / CAD",
    ]
    status_by_record = {
        labels[0]: "weak-IV ineligible",
        labels[1]: "no comparable credible sets",
        labels[2]: "no comparable credible sets",
        labels[3]: "no comparable credible sets",
        labels[4]: "distinct signals in all comparisons",
        labels[5]: "distinct signals in all comparisons",
        labels[6]: "distinct signals in all comparisons",
        labels[7]: "SuSiE compatibility failure",
    }
    if len(validation) != 8:
        raise ValueError("Validation-status file no longer has the expected eight audited records.")

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0), gridspec_kw={"width_ratios": [1.35, 1]}, constrained_layout=True)
    order = [
        "IFITM2 / GTEx fibroblast / CAD",
        "IFITM2 / OneK1K NK / CAD",
        "TMEM80 / Fairfax IFN24 / CAD",
    ]
    for index, record in enumerate(order):
        subset = pairs.loc[pairs["record"] == record]
        offsets = np.linspace(-0.18, 0.18, len(subset))
        axes[0].scatter(subset["PP.H3.abf"], np.full(len(subset), index) + offsets, color="#D55E00", s=24, alpha=0.75, label="H3: distinct signals" if index == 0 else None)
        axes[0].scatter(subset["PP.H4.abf"], np.full(len(subset), index) + offsets, color="#0072B2", s=24, alpha=0.75, label="H4: shared signal" if index == 0 else None)
    axes[0].set(
        title="A. Credible-set pair posteriors",
        xlabel="Posterior probability",
        xlim=(-0.02, 1.02),
        yticks=range(len(order)),
        yticklabels=order,
    )
    axes[0].grid(axis="x", alpha=0.22)
    axes[0].legend(frameon=False, fontsize=8, loc="center", bbox_to_anchor=(0.52, 0.50))

    axes[1].axis("off")
    axes[1].set_title("B. Final status of all audited records", loc="left")
    table_data = [[record, status_by_record[record]] for record in labels]
    table = axes[1].table(
        cellText=table_data,
        colLabels=["Record", "Final regional status"],
        colWidths=[0.51, 0.49],
        cellLoc="left",
        colLoc="left",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7.1)
    table.scale(1.0, 1.55)
    for (row, column), cell in table.get_celld().items():
        cell.set_edgecolor("#C8C8C8")
        if row == 0:
            cell.set_facecolor("#E8EEF5")
            cell.set_text_props(weight="bold")
        elif "distinct signals" in table_data[row - 1][1]:
            cell.set_facecolor("#FDE9E1")
        elif "failure" in table_data[row - 1][1] or "ineligible" in table_data[row - 1][1]:
            cell.set_facecolor("#F2F2F2")
    fig.suptitle("Regional colocalization audit: no record supports a shared causal signal", fontsize=12, y=1.03)
    save_figure(fig, out_dir / "figure_s2_regional_colocalization_audit")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--performance", required=True, type=Path)
    parser.add_argument("--coloc-root", required=True, type=Path)
    parser.add_argument("--validation", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    make_calibration_figure(args.performance, args.out_dir)
    make_coloc_figure(args.coloc_root, args.validation, args.out_dir)


if __name__ == "__main__":
    main()
