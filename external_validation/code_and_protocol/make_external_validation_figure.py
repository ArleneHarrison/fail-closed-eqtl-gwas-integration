#!/usr/bin/env python3
"""Build the submission figure and three-line table from real validation outputs.

The program is deliberately fail-closed: it refuses missing/extra factorial
units, missing required columns, inconsistent rerun flags, unknown states, and
implicit overwrite.  Synthetic fixtures are accepted only behind an explicit
test flag and are watermarked into a directory whose name contains
``synthetic_test``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


TISSUE_LABELS = {
    "QTD000131": "Aorta\n(QTD000131)",
    "QTD000136": "Coronary artery\n(QTD000136)",
    "QTD000251": "Heart—atrial appendage\n(QTD000251)",
    "QTD000256": "Heart—left ventricle\n(QTD000256)",
}
TISSUE_ORDER = list(TISSUE_LABELS)
OUTCOME_ORDER = ["CAD", "HF"]

SUMMARY_REQUIRED = {
    "unit_id", "region_id", "anchor_outcome", "tissue_id", "outcome_id",
    "pre_gate_status", "pre_gate_code", "gate_status", "error_codes",
    "variant_overlap_n", "ld_dimension", "ld_rank", "posterior_generated",
}
RERUN_REQUIRED = {"unit_id", "status", "byte_identical", "semantic_identical", "reason"}
REGIONS_REQUIRED = {
    "region_id", "anchor_outcome", "lead_variant", "chromosome",
    "lead_position", "region_start", "region_end",
}
GENES_REQUIRED = {
    "region_id", "gene_id", "gene_selection_status", "gene_selection_code",
    "n_shared_complete_candidates",
}
SYNTHETIC_COLUMN = "_synthetic_fixture"
OUTPUT_STEM = "external_validation_main_figure"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, label: str, required: set[str]) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{label} input does not exist: {path}")
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if frame.empty:
        raise ValueError(f"{label} input is empty: {path}")
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{label} input is missing required columns: {', '.join(missing)}")
    return frame


def parse_bool(series: pd.Series, label: str) -> pd.Series:
    mapping = {"true": True, "false": False, "1": True, "0": False}
    normalized = series.astype(str).str.strip().str.lower()
    invalid = sorted(set(normalized).difference(mapping))
    if invalid:
        raise ValueError(f"{label} contains non-boolean values: {invalid}")
    return normalized.map(mapping).astype(bool)


def parse_optional_numeric(series: pd.Series, label: str) -> pd.Series:
    stripped = series.astype(str).str.strip()
    values = pd.to_numeric(stripped.replace("", np.nan), errors="coerce")
    bad = stripped.ne("") & values.isna()
    if bad.any():
        examples = stripped.loc[bad].drop_duplicates().head(5).tolist()
        raise ValueError(f"{label} contains non-numeric values: {examples}")
    return values


def detect_synthetic(frames: Iterable[pd.DataFrame]) -> bool:
    flags: list[bool] = []
    for frame in frames:
        if SYNTHETIC_COLUMN in frame.columns:
            flags.extend(parse_bool(frame[SYNTHETIC_COLUMN], SYNTHETIC_COLUMN).tolist())
    return any(flags)


def validate_inputs(
    summary: pd.DataFrame,
    rerun: pd.DataFrame,
    regions: pd.DataFrame,
    genes: pd.DataFrame,
    expected_regions: int,
    expected_attempts: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    for label, frame, key in (
        ("summary", summary, "unit_id"),
        ("rerun", rerun, "unit_id"),
        ("regions", regions, "region_id"),
        ("genes", genes, "region_id"),
    ):
        if frame[key].eq("").any():
            raise ValueError(f"{label}.{key} contains empty identifiers")
        duplicated = frame.loc[frame[key].duplicated(keep=False), key].unique().tolist()
        if duplicated:
            raise ValueError(f"{label}.{key} contains duplicates: {duplicated[:5]}")

    if len(regions) != expected_regions:
        raise ValueError(f"expected {expected_regions} regions, observed {len(regions)}")
    if len(summary) != expected_attempts:
        raise ValueError(f"expected {expected_attempts} attempts, observed {len(summary)}")
    if len(rerun) != expected_attempts:
        raise ValueError(f"expected {expected_attempts} rerun records, observed {len(rerun)}")
    if set(summary["unit_id"]) != set(rerun["unit_id"]):
        raise ValueError("summary and rerun unit_id sets differ")
    if set(regions["region_id"]) != set(genes["region_id"]):
        raise ValueError("regions and genes region_id sets differ")
    if set(summary["region_id"]) != set(regions["region_id"]):
        raise ValueError("summary and regions region_id sets differ")
    if set(summary["tissue_id"]) != set(TISSUE_ORDER):
        raise ValueError(f"summary must contain exactly the four frozen tissue IDs: {TISSUE_ORDER}")
    if set(summary["outcome_id"]) != set(OUTCOME_ORDER):
        raise ValueError(f"summary must contain exactly the frozen outcomes: {OUTCOME_ORDER}")

    if not set(regions["anchor_outcome"]).issubset(set(OUTCOME_ORDER)):
        raise ValueError("regions.anchor_outcome contains an unknown outcome")
    if not set(summary["anchor_outcome"]).issubset(set(OUTCOME_ORDER)):
        raise ValueError("summary.anchor_outcome contains an unknown outcome")
    region_anchor = regions.set_index("region_id")["anchor_outcome"].to_dict()
    inconsistent_anchor = summary.apply(
        lambda row: row["anchor_outcome"] != region_anchor[row["region_id"]], axis=1
    )
    if inconsistent_anchor.any():
        raise ValueError("summary.anchor_outcome disagrees with regions.tsv")

    expected_units = {
        (region, tissue, outcome)
        for region in regions["region_id"]
        for tissue in TISSUE_ORDER
        for outcome in OUTCOME_ORDER
    }
    observed_units = set(summary[["region_id", "tissue_id", "outcome_id"]].itertuples(index=False, name=None))
    if observed_units != expected_units:
        raise ValueError("summary does not contain exactly one region × tissue × outcome factorial record")

    allowed_gene = {"SELECTED", "INELIGIBLE"}
    if not set(genes["gene_selection_status"]).issubset(allowed_gene):
        raise ValueError("genes.tsv contains an unknown gene_selection_status")
    allowed_gate = {"READY", "INELIGIBLE", "MISSING"}
    if not set(summary["gate_status"]).issubset(allowed_gate):
        raise ValueError("summary contains an unknown gate_status")
    allowed_pre = {"PENDING", "INELIGIBLE"}
    if not set(summary["pre_gate_status"]).issubset(allowed_pre):
        raise ValueError("summary contains an unknown pre_gate_status")
    if not set(rerun["status"]).issubset({"PASS", "FAIL"}):
        raise ValueError("rerun contains an unknown status")

    rerun = rerun.copy()
    rerun["byte_identical_bool"] = parse_bool(rerun["byte_identical"], "rerun.byte_identical")
    rerun["semantic_identical_bool"] = parse_bool(rerun["semantic_identical"], "rerun.semantic_identical")
    consistent_rerun = rerun["status"].eq("PASS") == (
        rerun["byte_identical_bool"] & rerun["semantic_identical_bool"]
    )
    if not consistent_rerun.all():
        raise ValueError("rerun.status disagrees with byte_identical/semantic_identical")

    summary = summary.copy()
    summary["posterior_generated_bool"] = parse_bool(
        summary["posterior_generated"], "summary.posterior_generated"
    )
    summary["variant_overlap_n_num"] = parse_optional_numeric(
        summary["variant_overlap_n"], "summary.variant_overlap_n"
    )
    summary["ld_dimension_num"] = parse_optional_numeric(
        summary["ld_dimension"], "summary.ld_dimension"
    )
    summary["ld_rank_num"] = parse_optional_numeric(summary["ld_rank"], "summary.ld_rank")
    negative = (
        summary[["variant_overlap_n_num", "ld_dimension_num", "ld_rank_num"]]
        .lt(0).any(axis=1)
    )
    if negative.any():
        raise ValueError("overlap/LD counts cannot be negative")
    impossible_rank = (
        summary["ld_rank_num"].notna()
        & summary["ld_dimension_num"].notna()
        & summary["ld_rank_num"].gt(summary["ld_dimension_num"])
    )
    if impossible_rank.any():
        raise ValueError("ld_rank exceeds ld_dimension")
    nonready_posterior = summary["gate_status"].ne("READY") & summary["posterior_generated_bool"]
    if nonready_posterior.any():
        raise ValueError("non-READY attempts contain posterior artifacts")
    return summary, rerun


def ensure_output_policy(out_dir: Path, synthetic: bool, allow_synthetic: bool, force: bool) -> list[Path]:
    if synthetic:
        if not allow_synthetic:
            raise ValueError("synthetic input refused; use --allow-synthetic-test-output for tests only")
        if "synthetic_test" not in out_dir.name.lower():
            raise ValueError("synthetic output directory name must contain 'synthetic_test'")
    elif allow_synthetic:
        raise ValueError("--allow-synthetic-test-output was supplied to non-synthetic inputs")
    targets = [
        out_dir / f"{OUTPUT_STEM}.png",
        out_dir / f"{OUTPUT_STEM}.svg",
        out_dir / f"{OUTPUT_STEM}.pdf",
        out_dir / "external_validation_summary_table.tsv",
        out_dir / "external_validation_summary_table.tex",
        out_dir / "external_validation_figure_derived_data.tsv",
        out_dir / "external_validation_alt_text.txt",
        out_dir / "external_validation_output_manifest.json",
    ]
    existing = [path for path in targets if path.exists()]
    if existing and not force:
        raise FileExistsError("refusing implicit overwrite: " + ", ".join(str(path) for path in existing))
    out_dir.mkdir(parents=True, exist_ok=True)
    return targets


def latex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in text)


def build_table(summary: pd.DataFrame, rerun: pd.DataFrame, regions: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    counts = Counter(summary["gate_status"])
    rerun_counts = Counter(rerun["status"])
    selected = int(genes["gene_selection_status"].eq("SELECTED").sum())
    nonready_posteriors = int(
        (summary["gate_status"].ne("READY") & summary["posterior_generated_bool"]).sum()
    )
    rows: list[dict[str, object]] = [
        {"Category": "Resource", "Item": "eQTL tissues", "Specification": ", ".join(TISSUE_ORDER), "Count": 4},
        {"Category": "Resource", "Item": "GWAS outcomes", "Specification": "CAD and HF", "Count": 2},
    ]
    for outcome in OUTCOME_ORDER:
        rows.append({
            "Category": "Design", "Item": f"{outcome}-anchored regions",
            "Specification": "Frozen GWAS-only anchor selection",
            "Count": int(regions["anchor_outcome"].eq(outcome).sum()),
        })
    rows.extend([
        {"Category": "Design", "Item": "Factorial attempts", "Specification": "regions × 4 tissues × 2 outcomes", "Count": len(summary)},
        {"Category": "Result", "Item": "Gene-selected regions", "Specification": "Shared four-tissue completeness rule", "Count": selected},
        {"Category": "Result", "Item": "Gene-ineligible regions", "Specification": "No eligible shared gene", "Count": len(genes) - selected},
        {"Category": "Result", "Item": "Gate READY attempts", "Specification": "All hard contracts passed", "Count": counts.get("READY", 0)},
        {"Category": "Result", "Item": "Gate INELIGIBLE attempts", "Specification": "Fail-closed stop", "Count": counts.get("INELIGIBLE", 0)},
        {"Category": "Result", "Item": "Missing gate records", "Specification": "No gate record available", "Count": counts.get("MISSING", 0)},
        {"Category": "Reproducibility", "Item": "Exact rerun PASS", "Specification": "Byte and semantic identity", "Count": rerun_counts.get("PASS", 0)},
        {"Category": "Reproducibility", "Item": "Exact rerun FAIL", "Specification": "Byte or semantic mismatch", "Count": rerun_counts.get("FAIL", 0)},
        {"Category": "Integrity", "Item": "Non-READY posterior artifacts", "Specification": "Required value: 0", "Count": nonready_posteriors},
    ])
    return pd.DataFrame(rows, columns=["Category", "Item", "Specification", "Count"])


def write_three_line_table(frame: pd.DataFrame, tsv_path: Path, tex_path: Path) -> None:
    frame.to_csv(tsv_path, sep="\t", index=False, lineterminator="\n")
    lines = [
        r"% Requires \usepackage{booktabs,tabularx}.",
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Frozen external technical-validation design and observed technical counts.}",
        r"\label{tab:external_validation_summary}",
        r"\begin{tabularx}{\textwidth}{@{}llXr@{}}",
        r"\toprule",
        r"Category & Item & Specification & Count \\",
        r"\midrule",
    ]
    for row in frame.itertuples(index=False):
        values = [latex_escape(value) for value in row]
        lines.append(" & ".join(values) + r" \\")
    lines.extend([r"\bottomrule", r"\end{tabularx}", r"\end{table*}", ""])
    tex_path.write_text("\n".join(lines), encoding="utf-8")


def add_panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(-0.13, 1.08, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")


def build_figure(
    summary: pd.DataFrame,
    rerun: pd.DataFrame,
    genes: pd.DataFrame,
    synthetic: bool,
) -> tuple[plt.Figure, pd.DataFrame, str]:
    gate_counts = Counter(summary["gate_status"])
    rerun_counts = Counter(rerun["status"])
    gene_eligible_regions = set(genes.loc[genes["gene_selection_status"].eq("SELECTED"), "region_id"])
    gene_eligible_attempts = int(summary["region_id"].isin(gene_eligible_regions).sum())

    derived: list[dict[str, object]] = []
    fig = plt.figure(figsize=(174 / 25.4, 150 / 25.4), layout="constrained", facecolor="white")
    grid = fig.add_gridspec(2, 2, width_ratios=[0.88, 1.12], height_ratios=[0.92, 1.08])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    # A: zero-baseline technical checkpoint counts. These are explicitly not
    # drawn as a survival funnel because rerun status has a different meaning.
    checkpoint_labels = ["All\nattempts", "Gene-eligible\nattempts", "Gate\nREADY", "Rerun\nPASS"]
    checkpoint_names = ["Planned attempts", "Gene-eligible pre-gate", "Gate READY", "Exact rerun PASS"]
    checkpoint_values = [len(summary), gene_eligible_attempts, gate_counts.get("READY", 0), rerun_counts.get("PASS", 0)]
    checkpoint_colors = ["#777777", "#56B4E9", "#0072B2", "#009E73"]
    checkpoint_hatches = ["//", "..", "", "xx"]
    bars = ax_a.bar(range(4), checkpoint_values, color=checkpoint_colors, edgecolor="black", linewidth=0.7)
    for bar, hatch, value, label in zip(bars, checkpoint_hatches, checkpoint_values, checkpoint_names):
        bar.set_hatch(hatch)
        ax_a.text(bar.get_x() + bar.get_width() / 2, value + max(checkpoint_values) * 0.025, str(value), ha="center", va="bottom", fontsize=7)
        derived.append({"panel": "A", "record_type": "checkpoint", "group_1": label, "group_2": "", "value_1": value, "value_2": "", "unit_id": ""})
    ax_a.set_xticks(range(4), checkpoint_labels)
    ax_a.set_ylabel("Attempts, n")
    ax_a.set_ylim(0, max(checkpoint_values) * 1.18 if max(checkpoint_values) else 1)
    ax_a.set_title("Technical checkpoint counts", loc="left", fontsize=8)
    ax_a.spines[["top", "right"]].set_visible(False)
    ax_a.grid(axis="y", color="#D9D9D9", linewidth=0.5)
    add_panel_label(ax_a, "A")

    # B: READY / total is printed in every cell; color only provides redundant
    # magnitude encoding.
    ready_matrix = np.zeros((len(TISSUE_ORDER), len(OUTCOME_ORDER)), dtype=float)
    total_matrix = np.zeros_like(ready_matrix, dtype=int)
    ready_n_matrix = np.zeros_like(total_matrix)
    for row_idx, tissue in enumerate(TISSUE_ORDER):
        for col_idx, outcome in enumerate(OUTCOME_ORDER):
            subset = summary.loc[(summary["tissue_id"] == tissue) & (summary["outcome_id"] == outcome)]
            total = len(subset)
            ready = int(subset["gate_status"].eq("READY").sum())
            total_matrix[row_idx, col_idx] = total
            ready_n_matrix[row_idx, col_idx] = ready
            ready_matrix[row_idx, col_idx] = ready / total if total else np.nan
            derived.append({"panel": "B", "record_type": "ready_cell", "group_1": tissue, "group_2": outcome, "value_1": ready, "value_2": total, "unit_id": ""})
    image = ax_b.pcolormesh(
        np.arange(len(OUTCOME_ORDER) + 1) - 0.5,
        np.arange(len(TISSUE_ORDER) + 1) - 0.5,
        ready_matrix,
        vmin=0, vmax=1, cmap=mpl.colormaps["Blues"], shading="flat",
    )
    ax_b.set_xlim(-0.5, len(OUTCOME_ORDER) - 0.5)
    ax_b.set_ylim(len(TISSUE_ORDER) - 0.5, -0.5)
    ax_b.set_xticks(range(len(OUTCOME_ORDER)), OUTCOME_ORDER)
    ax_b.set_yticks(range(len(TISSUE_ORDER)), [TISSUE_LABELS[tissue] for tissue in TISSUE_ORDER])
    ax_b.xaxis.tick_top()
    ax_b.tick_params(axis="both", length=0)
    for row_idx in range(len(TISSUE_ORDER)):
        for col_idx in range(len(OUTCOME_ORDER)):
            value = ready_matrix[row_idx, col_idx]
            text_color = "white" if np.isfinite(value) and value >= 0.52 else "black"
            ax_b.text(col_idx, row_idx, f"{ready_n_matrix[row_idx, col_idx]}/{total_matrix[row_idx, col_idx]}", ha="center", va="center", color=text_color, fontsize=7, fontweight="bold")
    ax_b.text(
        0.5, -0.08, "Darker fill = higher READY proportion; cells print exact READY/total",
        transform=ax_b.transAxes, ha="center", va="top", fontsize=6.2,
    )
    ax_b.set_title("Gate READY / total", loc="left", fontsize=8, pad=22)
    add_panel_label(ax_b, "B")

    # C: only complete overlap/rank diagnostics are plotted. Incomplete records
    # remain counted and are stated in the panel instead of being silently lost.
    diagnostic = summary.loc[
        summary["variant_overlap_n_num"].notna()
        & summary["ld_dimension_num"].gt(0)
        & summary["ld_rank_num"].notna()
    ].copy()
    diagnostic["rank_fraction"] = diagnostic["ld_rank_num"] / diagnostic["ld_dimension_num"]
    marker_by_outcome = {"CAD": "o", "HF": "^"}
    color_by_outcome = {"CAD": "#0072B2", "HF": "#D55E00"}
    for outcome in OUTCOME_ORDER:
        for status in ["READY", "INELIGIBLE", "MISSING"]:
            subset = diagnostic.loc[(diagnostic["outcome_id"] == outcome) & (diagnostic["gate_status"] == status)]
            if subset.empty:
                continue
            face = color_by_outcome[outcome] if status == "READY" else "none"
            ax_c.scatter(
                subset["variant_overlap_n_num"], subset["rank_fraction"], s=20,
                marker=marker_by_outcome[outcome], facecolors=face,
                edgecolors=color_by_outcome[outcome], linewidths=0.75, alpha=0.72,
            )
    for row in diagnostic.itertuples(index=False):
        derived.append({"panel": "C", "record_type": "attempt_diagnostic", "group_1": row.outcome_id, "group_2": row.gate_status, "value_1": row.variant_overlap_n_num, "value_2": row.rank_fraction, "unit_id": row.unit_id})
    positive_x = diagnostic.loc[diagnostic["variant_overlap_n_num"] > 0, "variant_overlap_n_num"]
    use_log = len(positive_x) == len(diagnostic) and not positive_x.empty and positive_x.max() / positive_x.min() >= 100
    if use_log:
        ax_c.set_xscale("log")
        x_label = "Ordered variant overlap, n (log scale)"
    else:
        x_label = "Ordered variant overlap, n"
    ax_c.set_xlabel(x_label)
    ax_c.set_ylabel("LD numerical rank / dimension")
    ax_c.set_ylim(-0.03, 1.05)
    ax_c.spines[["top", "right"]].set_visible(False)
    ax_c.grid(color="#D9D9D9", linewidth=0.5)
    omitted_n = len(summary) - len(diagnostic)
    ax_c.text(
        0.98, 0.97, f"Complete: {len(diagnostic)}/{len(summary)}; incomplete: {omitted_n}",
        transform=ax_c.transAxes, fontsize=6, ha="right", va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.0},
    )
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color_by_outcome["CAD"], markeredgecolor=color_by_outcome["CAD"], label="CAD; READY", markersize=5),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=color_by_outcome["HF"], markeredgecolor=color_by_outcome["HF"], label="HF; READY", markersize=5),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor="#555555", label="Non-READY; open", markersize=5),
    ]
    ax_c.legend(handles=handles, loc="lower right", frameon=False, fontsize=6, handletextpad=0.4)
    ax_c.set_title("Ordered overlap versus LD rank", loc="left", fontsize=8)
    add_panel_label(ax_c, "C")

    # D: exact rerun counts, zero baseline, direct labels and hatch redundancy.
    rerun_labels = ["PASS", "FAIL"]
    rerun_values = [rerun_counts.get(label, 0) for label in rerun_labels]
    rerun_colors = ["#009E73", "#CC79A7"]
    rerun_hatches = ["", "///"]
    bars = ax_d.barh(range(2), rerun_values, color=rerun_colors, edgecolor="black", linewidth=0.7)
    for bar, hatch, value, label in zip(bars, rerun_hatches, rerun_values, rerun_labels):
        bar.set_hatch(hatch)
        ax_d.text(value + max(rerun_values + [1]) * 0.02, bar.get_y() + bar.get_height() / 2, str(value), va="center", fontsize=7)
        derived.append({"panel": "D", "record_type": "rerun_status", "group_1": label, "group_2": "", "value_1": value, "value_2": len(rerun), "unit_id": ""})
    ax_d.set_yticks(range(2), rerun_labels)
    ax_d.invert_yaxis()
    ax_d.set_xlabel("Attempts, n")
    ax_d.set_xlim(0, max(rerun_values + [1]) * 1.12)
    ax_d.spines[["top", "right"]].set_visible(False)
    ax_d.grid(axis="x", color="#D9D9D9", linewidth=0.5)
    ax_d.text(0.98, 0.08, "PASS requires byte and\nsemantic identity", transform=ax_d.transAxes, ha="right", va="bottom", fontsize=6.5)
    ax_d.set_title("Deterministic gate rerun", loc="left", fontsize=8)
    add_panel_label(ax_d, "D")

    if synthetic:
        fig.text(0.5, 0.5, "SYNTHETIC TEST ONLY — NOT FOR SUBMISSION", ha="center", va="center", rotation=30, fontsize=18, color="#A00000", alpha=0.22, fontweight="bold")

    fig.suptitle("Internally frozen multi-tissue cardiovascular technical evaluation", fontsize=9, fontweight="bold")
    alt_text = (
        f"Four-panel technical-validation figure for {len(summary)} attempts across four GTEx tissues "
        f"and CAD/HF. {gene_eligible_attempts} attempts were gene-eligible before gating; "
        f"{gate_counts.get('READY', 0)} were READY, {gate_counts.get('INELIGIBLE', 0)} were "
        f"INELIGIBLE, and {gate_counts.get('MISSING', 0)} lacked a gate record. The heat map "
        f"prints READY/total for every tissue-outcome cell. The diagnostic scatter displays "
        f"ordered overlap against LD numerical-rank fraction for {len(diagnostic)} complete records "
        f"and explicitly notes {omitted_n} incomplete records. Exact rerun passed for "
        f"{rerun_counts.get('PASS', 0)} and failed for {rerun_counts.get('FAIL', 0)} attempts."
    )
    return fig, pd.DataFrame(derived), alt_text


def write_manifest(out_dir: Path, inputs: dict[str, Path], synthetic: bool) -> None:
    output_paths = sorted(
        path for path in out_dir.iterdir()
        if path.is_file() and path.name != "external_validation_output_manifest.json"
    )
    png_path = out_dir / f"{OUTPUT_STEM}.png"
    from PIL import Image
    with Image.open(png_path) as image:
        pixels = {"width": image.width, "height": image.height, "dpi": list(image.info.get("dpi", (None, None)))}
    manifest = {
        "schema_version": "external_validation_figure_1.0",
        "data_origin": "SYNTHETIC_TEST_ONLY" if synthetic else "REAL_RUNNER_OUTPUT",
        "input_sha256": {label: sha256_file(path) for label, path in inputs.items()},
        "output_sha256": {path.name: sha256_file(path) for path in output_paths},
        "png_metadata": pixels,
        "figure_width_mm": 174,
        "figure_height_mm": 150,
        "png_export_dpi": 600,
        "transformations": [
            "counts grouped directly from input records",
            "READY proportions are READY count divided by observed factorial cell total",
            "LD rank fraction is numerical rank divided by LD dimension",
            "records without complete overlap/rank diagnostics remain explicitly counted as incomplete",
        ],
        "software": {
            "python": sys.version.split()[0], "platform": platform.platform(),
            "matplotlib": mpl.__version__, "numpy": np.__version__, "pandas": pd.__version__,
        },
    }
    (out_dir / "external_validation_output_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    root.add_argument("--summary", type=Path, required=True, help="runner summarize TSV (attempts plus gate diagnostics)")
    root.add_argument("--rerun", type=Path, required=True, help="compare-reruns TSV")
    root.add_argument("--regions", type=Path, required=True, help="frozen regions.tsv")
    root.add_argument("--genes", type=Path, required=True, help="frozen genes.tsv")
    root.add_argument("--out-dir", type=Path, required=True)
    root.add_argument("--expected-regions", type=int, default=24)
    root.add_argument("--expected-attempts", type=int, default=192)
    root.add_argument("--force", action="store_true", help="explicitly replace named outputs")
    root.add_argument("--allow-synthetic-test-output", action="store_true", help="tests only; forces watermark and synthetic_test directory name")
    return root


def main() -> None:
    args = parser().parse_args()
    if args.expected_regions <= 0 or args.expected_attempts <= 0:
        raise ValueError("expected counts must be positive")
    summary = read_tsv(args.summary, "summary", SUMMARY_REQUIRED)
    rerun = read_tsv(args.rerun, "rerun", RERUN_REQUIRED)
    regions = read_tsv(args.regions, "regions", REGIONS_REQUIRED)
    genes = read_tsv(args.genes, "genes", GENES_REQUIRED)
    frames = [summary, rerun, regions, genes]
    synthetic = detect_synthetic(frames)
    ensure_output_policy(args.out_dir, synthetic, args.allow_synthetic_test_output, args.force)
    summary, rerun = validate_inputs(
        summary, rerun, regions, genes, args.expected_regions, args.expected_attempts
    )

    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 8,
        "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5, "axes.linewidth": 0.6, "lines.linewidth": 0.8,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "savefig.transparent": False,
    })
    fig, derived, alt_text = build_figure(summary, rerun, genes, synthetic)
    png_path = args.out_dir / f"{OUTPUT_STEM}.png"
    svg_path = args.out_dir / f"{OUTPUT_STEM}.svg"
    pdf_path = args.out_dir / f"{OUTPUT_STEM}.pdf"
    fig.savefig(png_path, dpi=600, bbox_inches=None, facecolor="white")
    # Matplotlib's PNG backend writes RGBA even against an opaque figure patch.
    # Flatten only that redundant alpha channel onto the declared white
    # background; scientific marks and pixel dimensions remain unchanged.
    from PIL import Image
    with Image.open(png_path) as png:
        dpi = png.info.get("dpi", (600, 600))
        rgb = Image.new("RGB", png.size, "white")
        if png.mode == "RGBA":
            rgb.paste(png, mask=png.getchannel("A"))
        else:
            rgb.paste(png.convert("RGB"))
        rgb.save(png_path, format="PNG", dpi=dpi, optimize=False)
    fig.savefig(svg_path, bbox_inches=None, facecolor="white")
    fig.savefig(pdf_path, bbox_inches=None, facecolor="white")
    plt.close(fig)

    derived.to_csv(
        args.out_dir / "external_validation_figure_derived_data.tsv",
        sep="\t", index=False, lineterminator="\n",
    )
    (args.out_dir / "external_validation_alt_text.txt").write_text(alt_text + "\n", encoding="utf-8")
    table = build_table(summary, rerun, regions, genes)
    write_three_line_table(
        table,
        args.out_dir / "external_validation_summary_table.tsv",
        args.out_dir / "external_validation_summary_table.tex",
    )
    write_manifest(args.out_dir, {
        "summary": args.summary, "rerun": args.rerun,
        "regions": args.regions, "genes": args.genes,
    }, synthetic)
    print(json.dumps({
        "status": "PASS", "data_origin": "SYNTHETIC_TEST_ONLY" if synthetic else "REAL_RUNNER_OUTPUT",
        "output_directory": str(args.out_dir), "figure_stem": OUTPUT_STEM,
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
