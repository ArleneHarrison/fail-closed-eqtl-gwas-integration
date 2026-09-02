"""Construct the cell-state instrument-variable (IV) library for hcsMR.

Inputs:
  - The integrated atlas with `cell_state` and `cell_type` columns in obs.
  - A sc-eQTL summary statistics table (OneK1K format expected; columns:
    snp_id, chrom, pos, gene, cell_type, beta, se, p).
  - GTEx bulk eQTL (optional fallback for non-immune cardiac cell states).

Outputs:
  - `iv_library.parquet`: a tidy table with one row per (cell_state, gene, SNP)
    IV, including beta, se, F-statistic, sourceflag (sc-eqtl vs bulk-fallback).
  - `cell_state_iv_map.json`: per-cell-state list of (gene, SNP) tuples.

The cross-source cell-state mapping (e.g., OneK1K monocyte → cardiac
macrophage) is handled by `STATE_TO_SCEQTL_CELLTYPE` below — a small, curated
mapping that the researcher should review.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# Curated mapping: cardiac cell-state name (or substring) → OneK1K cell type.
# Cardiac cell states whose name STARTS with one of these prefixes are mapped to
# the listed OneK1K cell types for sc-eQTL IV ingestion.
STATE_TO_SCEQTL_CELLTYPE = {
    "macrophage": ["Mono_C", "Mono_NC", "DC"],
    "foamcell_macrophage": ["Mono_C", "Mono_NC"],
    "monocyte": ["Mono_C", "Mono_NC"],
    "T_cell": ["CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET", "CD8_S100B"],
    "B_cell": ["B_Mem", "B_IN", "Plasma"],
    "NK_cell": ["NK", "NK_R"],
    "mast_cell": [],  # not in OneK1K — will fall back to GTEx
    # All non-immune cardiac states fall back to GTEx bulk
    "cardiomyocyte": [],
    "stressed_cardiomyocyte": [],
    "atrial_cardiomyocyte": [],
    "ventricular_cardiomyocyte": [],
    "pacemaker_cardiomyocyte": [],
    "fibroblast": [],
    "activated_fibroblast": [],
    "endothelial": [],
    "lymphatic_endothelial": [],
    "smooth_muscle": [],
    "pericyte": [],
    "adipocyte": [],
    "neuronal": [],
    "unknown": [],
}

# Cell states mapped to GTEx tissue eQTL as a downweighted fallback
STATE_TO_GTEX_TISSUE = {
    "cardiomyocyte": "Heart_Left_Ventricle",
    "stressed_cardiomyocyte": "Heart_Left_Ventricle",
    "atrial_cardiomyocyte": "Heart_Atrial_Appendage",
    "ventricular_cardiomyocyte": "Heart_Left_Ventricle",
    "pacemaker_cardiomyocyte": "Heart_Atrial_Appendage",
    "fibroblast": "Heart_Left_Ventricle",
    "activated_fibroblast": "Heart_Left_Ventricle",
    "endothelial": "Artery_Coronary",
    "lymphatic_endothelial": "Artery_Coronary",
    "smooth_muscle": "Artery_Coronary",
    "pericyte": "Artery_Coronary",
    "mast_cell": "Heart_Left_Ventricle",
    "adipocyte": "Heart_Left_Ventricle",
    "neuronal": "Heart_Left_Ventricle",
}


def load_onek1k(eqtl_path: str | Path, *, fdr_threshold: float = 0.05) -> pd.DataFrame:
    """Load OneK1K combined eQTL table.  Format: 14-cell-type tidy table."""
    eqtl_path = Path(eqtl_path)
    logger.info("loading OneK1K from %s", eqtl_path)
    df = pd.read_csv(eqtl_path, sep="\t", compression="infer", low_memory=False)
    logger.info("loaded %d rows; columns: %s", len(df), df.columns.tolist())
    # Heuristic normalisation of column names — adapt to the actual schema
    rename = {}
    for c in df.columns:
        cl = c.lower()
        if cl in ("snp", "snp_id", "rsid", "variant_id"):
            rename[c] = "snp"
        elif cl in ("chr", "chrom", "chromosome"):
            rename[c] = "chrom"
        elif cl in ("pos", "position", "bp"):
            rename[c] = "pos"
        elif cl in ("gene", "gene_id", "ensembl_id", "symbol", "gene_name", "ensembl"):
            rename[c] = "gene"
        elif cl in ("celltype", "cell_type", "cell-type", "cluster"):
            rename[c] = "cell_type"
        elif cl in ("beta", "effect", "slope"):
            rename[c] = "beta"
        elif cl in ("se", "stderr", "std_err", "slope_se"):
            rename[c] = "se"
        elif cl in ("p", "pval", "p_value", "pvalue"):
            rename[c] = "p"
        elif cl in ("fdr", "qval", "q_value"):
            rename[c] = "fdr"
    df = df.rename(columns=rename)
    # Filter to significant hits (column may not exist)
    if "fdr" in df.columns:
        df = df[df["fdr"] < fdr_threshold].copy()
        logger.info("after FDR<%.2f filter: %d rows", fdr_threshold, len(df))
    elif "p" in df.columns:
        df = df[df["p"] < 5e-8].copy()
        logger.info("after p<5e-8 filter: %d rows", len(df))
    return df


def load_gtex_signif(path: str | Path, tissue: str) -> pd.DataFrame:
    """Load a GTEx significant_variant_gene_pairs file for a single tissue."""
    path = Path(path)
    logger.info("loading GTEx %s from %s", tissue, path)
    df = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    df["tissue"] = tissue
    # GTEx column names: gene_id, variant_id (chr_pos_ref_alt_b38), pval_nominal, slope, slope_se, ...
    rename = {
        "gene_id": "gene",
        "variant_id": "snp",
        "pval_nominal": "p",
        "slope": "beta",
        "slope_se": "se",
    }
    return df.rename(columns={k: v for k, v in rename.items() if k in df.columns})


def f_statistic(beta: np.ndarray, se: np.ndarray) -> np.ndarray:
    return np.where(se > 0, (beta / se) ** 2, 0.0)


def build_iv_library(
    cell_states: Iterable[str],
    state_to_type: dict[str, str],
    onek1k_df: pd.DataFrame | None,
    gtex_dict: dict[str, pd.DataFrame] | None,
    *,
    f_threshold: float = 10.0,
    min_ivs_per_state_gene: int = 3,
    bulk_downweight: float = 2.0,
) -> pd.DataFrame:
    """Assemble a tidy (cell_state, gene, snp, beta, se, F, source) table.

    For each cell_state c:
      - If state has matching OneK1K cell types, use those sc-eQTLs.
      - Else if state has matching GTEx tissue, use those bulk eQTLs with se *
        bulk_downweight (Bayesian "weakly informative").
    """
    rows = []
    for state in cell_states:
        parental = state_to_type.get(state, state).split("__")[0]
        sc_targets = STATE_TO_SCEQTL_CELLTYPE.get(parental, [])
        gtex_tissue = STATE_TO_GTEX_TISSUE.get(parental)

        used_sc = False
        if onek1k_df is not None and sc_targets and "cell_type" in onek1k_df.columns:
            sub = onek1k_df[onek1k_df["cell_type"].isin(sc_targets)]
            if len(sub):
                used_sc = True
                logger.info("state %s: %d sc-eQTL IVs from OneK1K %s", state, len(sub), sc_targets)
                for _, r in sub.iterrows():
                    beta = float(r.get("beta", np.nan))
                    se = float(r.get("se", np.nan))
                    if np.isnan(beta) or np.isnan(se) or se <= 0:
                        continue
                    F = (beta / se) ** 2
                    if F < f_threshold:
                        continue
                    rows.append({
                        "cell_state": state,
                        "cell_type": parental,
                        "gene": r.get("gene"),
                        "snp": r.get("snp"),
                        "beta": beta,
                        "se": se,
                        "F": float(F),
                        "source": "sc-eQTL_OneK1K",
                        "source_cell_type": r.get("cell_type"),
                    })

        if not used_sc and gtex_tissue and gtex_dict is not None and gtex_tissue in gtex_dict:
            gx = gtex_dict[gtex_tissue]
            if len(gx):
                logger.info("state %s: falling back to GTEx %s (%d eQTLs)", state, gtex_tissue, len(gx))
                for _, r in gx.iterrows():
                    beta = float(r.get("beta", np.nan))
                    se = float(r.get("se", np.nan)) * bulk_downweight  # downweight!
                    if np.isnan(beta) or np.isnan(se) or se <= 0:
                        continue
                    F = (beta / se) ** 2
                    if F < f_threshold:
                        continue
                    rows.append({
                        "cell_state": state,
                        "cell_type": parental,
                        "gene": r.get("gene"),
                        "snp": r.get("snp"),
                        "beta": beta,
                        "se": se,
                        "F": float(F),
                        "source": "GTEx_bulk_downweighted",
                        "source_cell_type": gtex_tissue,
                    })

    df = pd.DataFrame(rows)
    if df.empty:
        logger.warning("IV library is empty!")
        return df
    # Filter (cell_state, gene) groups with too few IVs
    grp = df.groupby(["cell_state", "gene"]).size().reset_index(name="n_iv")
    keep = grp[grp["n_iv"] >= min_ivs_per_state_gene][["cell_state", "gene"]]
    df = df.merge(keep, on=["cell_state", "gene"], how="inner")
    logger.info("final IV library: %d IVs across %d (cell_state, gene) pairs",
                len(df), df.groupby(["cell_state", "gene"]).ngroups)
    return df


def write_iv_library(df: pd.DataFrame, out_path: str | Path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)
    logger.info("wrote IV library to %s (%d rows)", out_path, len(df))
