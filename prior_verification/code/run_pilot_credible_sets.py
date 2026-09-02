"""Real-data hcsMR pilot using EBI eQTL Catalogue credible_sets format.

Inputs:
  - HERMES HF GWAS sumstats (zipped txt)
  - Directory of credible_sets.tsv.gz files (OneK1K + GTEx + supporting)

For each cell state (= one credible_sets file), build the IV set, merge with
GWAS by rsid, and fit hcsMR. Then aggregate across all cell states.

Cell-state structure (curated, hierarchical):
  parental cell type → cell state list

  macrophage:
    - Alasoo_naive, Alasoo_IFNg, Alasoo_Salmonella, Alasoo_IFNg_Salmonella
  monocyte:
    - Fairfax_naive, Fairfax_IFN24, Fairfax_LPS2, Fairfax_LPS24
    - BLUEPRINT_naive
    - OneK1K_CD14_Mono, OneK1K_CD16_Mono
  NK:
    - OneK1K_NK, OneK1K_NK_CD56bright, OneK1K_NK_Proliferating
  T_cell:
    - OneK1K_CD4_Naive, OneK1K_CD4_CTL, OneK1K_CD4_TCM, OneK1K_CD4_TEM
    - OneK1K_CD8_Naive, OneK1K_CD8_TCM, OneK1K_CD8_TEM
    - OneK1K_Treg, OneK1K_MAIT, OneK1K_gdT, OneK1K_dnT
  B_cell:
    - OneK1K_B_naive, OneK1K_B_memory, OneK1K_B_intermediate, OneK1K_Plasmablast
  DC:
    - OneK1K_cDC2, OneK1K_pDC
  HSPC:
    - OneK1K_HSPC
  cardiac_bulk_LV:
    - GTEx_Heart_Left_Ventricle
  cardiac_bulk_AA:
    - GTEx_Heart_Atrial_Appendage
  vascular_coronary:
    - GTEx_Artery_Coronary
  vascular_aorta:
    - GTEx_Artery_Aorta
  fibroblast:
    - GTEx_Cells_Cultured_fibroblasts
  blood:
    - GTEx_Whole_Blood

Total: 31 cell states across 11 parental types.
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("JAX_PLATFORMS", "cpu")
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from hcsmr.gwas_io import load_hermes_hf
from hcsmr.model import fit_hcsmr, bayes_factor_zero

# State -> file mapping (relative to data root)
CELL_STATES = [
    # (state_name, parental_type, relative_path)
    ("Alasoo_macrophage_naive", "macrophage", "sceqtl/Supporting_credible_sets/Alasoo_macrophage_naive_QTD000001.credible_sets.tsv.gz"),
    ("Alasoo_macrophage_IFNg", "macrophage", "sceqtl/Supporting_credible_sets/Alasoo_macrophage_IFNg_QTD000006.credible_sets.tsv.gz"),
    ("Alasoo_macrophage_Salmonella", "macrophage", "sceqtl/Supporting_credible_sets/Alasoo_macrophage_Salmonella_QTD000011.credible_sets.tsv.gz"),
    ("Alasoo_macrophage_IFNg_Salmonella", "macrophage", "sceqtl/Supporting_credible_sets/Alasoo_macrophage_IFNg_Salmonella_QTD000016.credible_sets.tsv.gz"),
    ("BLUEPRINT_monocyte_naive", "monocyte", "sceqtl/Supporting_credible_sets/BLUEPRINT_monocyte_naive_QTD000021.credible_sets.tsv.gz"),
    ("Fairfax_monocyte_naive", "monocyte", "sceqtl/Supporting_credible_sets/Fairfax_monocyte_naive_QTD000081.credible_sets.tsv.gz"),
    ("Fairfax_monocyte_IFN24", "monocyte", "sceqtl/Supporting_credible_sets/Fairfax_monocyte_IFN24_QTD000082.credible_sets.tsv.gz"),
    ("Fairfax_monocyte_LPS2", "monocyte", "sceqtl/Supporting_credible_sets/Fairfax_monocyte_LPS2_QTD000083.credible_sets.tsv.gz"),
    ("Fairfax_monocyte_LPS24", "monocyte", "sceqtl/Supporting_credible_sets/Fairfax_monocyte_LPS24_QTD000084.credible_sets.tsv.gz"),
    ("OneK1K_CD14_Mono", "monocyte", "sceqtl/OneK1K_credible_sets/OneK1K_CD14_Mono_QTD000609.credible_sets.tsv.gz"),
    ("OneK1K_CD16_Mono", "monocyte", "sceqtl/OneK1K_credible_sets/OneK1K_CD16_Mono_QTD000610.credible_sets.tsv.gz"),
    ("OneK1K_NK", "NK", "sceqtl/OneK1K_credible_sets/OneK1K_NK_QTD000620.credible_sets.tsv.gz"),
    ("OneK1K_NK_CD56bright", "NK", "sceqtl/OneK1K_credible_sets/OneK1K_NK_CD56bright_QTD000621.credible_sets.tsv.gz"),
    ("OneK1K_NK_Proliferating", "NK", "sceqtl/OneK1K_credible_sets/OneK1K_NK_Proliferating_QTD000622.credible_sets.tsv.gz"),
    ("OneK1K_CD4_Naive", "T_cell", "sceqtl/OneK1K_credible_sets/OneK1K_CD4_Naive_QTD000612.credible_sets.tsv.gz"),
    ("OneK1K_CD4_CTL", "T_cell", "sceqtl/OneK1K_credible_sets/OneK1K_CD4_CTL_QTD000611.credible_sets.tsv.gz"),
    ("OneK1K_CD4_TCM", "T_cell", "sceqtl/OneK1K_credible_sets/OneK1K_CD4_TCM_QTD000613.credible_sets.tsv.gz"),
    ("OneK1K_CD4_TEM", "T_cell", "sceqtl/OneK1K_credible_sets/OneK1K_CD4_TEM_QTD000614.credible_sets.tsv.gz"),
    ("OneK1K_CD8_Naive", "T_cell", "sceqtl/OneK1K_credible_sets/OneK1K_CD8_Naive_QTD000615.credible_sets.tsv.gz"),
    ("OneK1K_CD8_TCM", "T_cell", "sceqtl/OneK1K_credible_sets/OneK1K_CD8_TCM_QTD000616.credible_sets.tsv.gz"),
    ("OneK1K_CD8_TEM", "T_cell", "sceqtl/OneK1K_credible_sets/OneK1K_CD8_TEM_QTD000617.credible_sets.tsv.gz"),
    ("OneK1K_Treg", "T_cell", "sceqtl/OneK1K_credible_sets/OneK1K_Treg_QTD000625.credible_sets.tsv.gz"),
    ("OneK1K_MAIT", "T_cell", "sceqtl/OneK1K_credible_sets/OneK1K_MAIT_QTD000619.credible_sets.tsv.gz"),
    ("OneK1K_gdT", "T_cell", "sceqtl/OneK1K_credible_sets/OneK1K_gdT_QTD000628.credible_sets.tsv.gz"),
    ("OneK1K_dnT", "T_cell", "sceqtl/OneK1K_credible_sets/OneK1K_dnT_QTD000627.credible_sets.tsv.gz"),
    ("OneK1K_B_naive", "B_cell", "sceqtl/OneK1K_credible_sets/OneK1K_B_naive_QTD000608.credible_sets.tsv.gz"),
    ("OneK1K_B_memory", "B_cell", "sceqtl/OneK1K_credible_sets/OneK1K_B_memory_QTD000607.credible_sets.tsv.gz"),
    ("OneK1K_B_intermediate", "B_cell", "sceqtl/OneK1K_credible_sets/OneK1K_B_intermediate_QTD000606.credible_sets.tsv.gz"),
    ("OneK1K_Plasmablast", "B_cell", "sceqtl/OneK1K_credible_sets/OneK1K_Plasmablast_QTD000623.credible_sets.tsv.gz"),
    ("OneK1K_cDC2", "DC", "sceqtl/OneK1K_credible_sets/OneK1K_cDC2_QTD000626.credible_sets.tsv.gz"),
    ("OneK1K_pDC", "DC", "sceqtl/OneK1K_credible_sets/OneK1K_pDC_QTD000629.credible_sets.tsv.gz"),
    ("GTEx_Heart_LV", "cardiac_bulk", "bulk_eqtl/GTEx_v8_credible_sets/GTEx_Heart_Left_Ventricle_QTD000256.credible_sets.tsv.gz"),
    ("GTEx_Heart_AA", "cardiac_bulk", "bulk_eqtl/GTEx_v8_credible_sets/GTEx_Heart_Atrial_Appendage_QTD000251.credible_sets.tsv.gz"),
    ("GTEx_Artery_Coronary", "vascular", "bulk_eqtl/GTEx_v8_credible_sets/GTEx_Artery_Coronary_QTD000136.credible_sets.tsv.gz"),
    ("GTEx_Artery_Aorta", "vascular", "bulk_eqtl/GTEx_v8_credible_sets/GTEx_Artery_Aorta_QTD000131.credible_sets.tsv.gz"),
    ("GTEx_Fibroblast", "fibroblast", "bulk_eqtl/GTEx_v8_credible_sets/GTEx_Cells_Cultured_fibroblasts_QTD000216.credible_sets.tsv.gz"),
    ("GTEx_Whole_Blood", "blood", "bulk_eqtl/GTEx_v8_credible_sets/GTEx_Whole_Blood_QTD000356.credible_sets.tsv.gz"),
]


def load_credible_sets(path):
    """Return a DataFrame with columns: gene_id, rsid, variant, pip, beta, se, z, pvalue, cs_id."""
    df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--hermes", required=True)
    p.add_argument("--data-root", required=True, help="root containing sceqtl/ and bulk_eqtl/")
    p.add_argument("--out", required=True)
    p.add_argument("--top-n-genes", type=int, default=300, help="genes ranked by total IV count")
    p.add_argument("--min-ivs-per-state", type=int, default=3)
    p.add_argument("--n-warmup", type=int, default=300)
    p.add_argument("--n-samples", type=int, default=300)
    p.add_argument("--n-chains", type=int, default=1)
    p.add_argument("--pip-threshold", type=float, default=0.01)
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    log = logging.getLogger("pilot")
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    data_root = Path(args.data_root)

    # ----- Stage A: load HERMES -----
    log.info("Loading HERMES HF...")
    gwas = load_hermes_hf(args.hermes)
    log.info("HERMES: %d rows; cols=%s", len(gwas), gwas.columns.tolist())
    # Normalise rsid column name (HERMES uses 'snp' from gwas_io._normalise_gwas)
    if "snp" not in gwas.columns:
        log.error("HERMES has no rsid-like column; abort")
        return
    # HERMES uses 'b' for effect; gwas_io._normalise_gwas doesn't map 'b'.
    beta_col = "beta" if "beta" in gwas.columns else ("b" if "b" in gwas.columns else None)
    if beta_col is None:
        log.error("GWAS missing both 'beta' and 'b' columns; abort")
        return
    gwas_lite = gwas[["snp", beta_col, "se", "p"]].rename(columns={
        "snp": "rsid", beta_col: "gwas_beta", "se": "gwas_se", "p": "gwas_p"
    })
    gwas_lite["gwas_beta"] = pd.to_numeric(gwas_lite["gwas_beta"], errors="coerce")
    gwas_lite["gwas_se"] = pd.to_numeric(gwas_lite["gwas_se"], errors="coerce")
    gwas_lite = gwas_lite.dropna(subset=["gwas_beta", "gwas_se"])
    log.info("HERMES (clean): %d rows", len(gwas_lite))

    # ----- Stage B: load each cell state's credible_sets and merge with GWAS -----
    merged_per_state = {}
    for state_name, parental, rel in CELL_STATES:
        path = data_root / rel
        if not path.exists():
            log.warning("missing %s", path)
            continue
        try:
            cs = load_credible_sets(path)
        except Exception as e:
            log.warning("load failed for %s: %s", state_name, e)
            continue
        # Filter by pip
        cs = cs[cs["pip"] >= args.pip_threshold].copy()
        if len(cs) == 0:
            log.info("state %s: 0 IVs after pip filter, skipping", state_name)
            continue
        # Keep only the TOP-PIP SNP per credible set per gene — this gives one
        # LD-independent IV per credible-set, which is what we want for MR
        cs = cs.sort_values("pip", ascending=False).groupby(["gene_id", "cs_id"]).head(1)
        merged = cs.merge(gwas_lite, on="rsid", how="inner")
        if len(merged) < args.min_ivs_per_state:
            log.info("state %s: only %d IVs after GWAS merge, skipping", state_name, len(merged))
            continue
        merged["state"] = state_name
        merged["parental"] = parental
        merged_per_state[state_name] = merged
        log.info("state %-32s (%s): %d IVs after pip>=%.2f and GWAS merge",
                 state_name, parental, len(merged), args.pip_threshold)

    if not merged_per_state:
        log.error("no cell state has any IVs after merge — abort")
        return

    # Save per-state merged tables
    merged_all = pd.concat(merged_per_state.values(), ignore_index=True)
    merged_all.to_csv(out / "merged_iv_gwas.csv.gz", index=False, compression="gzip")
    log.info("merged_all: %d rows across %d states", len(merged_all), len(merged_per_state))

    # ----- Stage C: identify genes with IVs across many states -----
    gene_state_counts = merged_all.groupby(["gene_id"])["state"].nunique().sort_values(ascending=False)
    log.info("top 10 genes by state-coverage:\n%s", gene_state_counts.head(10).to_string())
    top_genes = gene_state_counts.head(args.top_n_genes).index.tolist()

    # ----- Stage D: build the cell-state hierarchy -----
    state_list = sorted(merged_per_state.keys())
    parental_list = sorted(set(s[1] for s in CELL_STATES if s[0] in state_list))
    state_to_type = np.array([parental_list.index(merged_per_state[s].iloc[0]["parental"])
                               for s in state_list], dtype=np.int32)
    n_types = len(parental_list)
    log.info("cell states: %d; parental types: %d", len(state_list), n_types)

    # ----- Stage E: fit hcsMR per gene -----
    results = []
    fit_count = 0; skip_count = 0
    MAX_IVS_PER_GENE = 40   # hard cap to keep model size tractable
    for gi, g in enumerate(top_genes):
        gdf = merged_all[merged_all["gene_id"] == g]
        if len(gdf) < 3:
            skip_count += 1
            continue
        # Cap at top MAX_IVS_PER_GENE by pip (avoid huge J)
        if len(gdf) > MAX_IVS_PER_GENE:
            gdf = gdf.sort_values("pip", ascending=False).head(MAX_IVS_PER_GENE)
        # Build (J, C) arrays
        snps = gdf["rsid"].unique().tolist()
        J = len(snps); C = len(state_list)
        beta_hat = np.full((J, C), np.nan, dtype=np.float32)
        se_X = np.full((J, C), np.nan, dtype=np.float32)
        gamma_hat = np.zeros(J, dtype=np.float32)
        se_Y = np.zeros(J, dtype=np.float32)
        iv_in_state = np.zeros((J, C), dtype=bool)
        rsid_to_j = {snp: j for j, snp in enumerate(snps)}
        for _, r in gdf.iterrows():
            j = rsid_to_j[r["rsid"]]
            gamma_hat[j] = float(r["gwas_beta"])
            se_Y[j] = float(r["gwas_se"])
            s_idx = state_list.index(r["state"])
            try:
                beta_hat[j, s_idx] = float(r["beta"])
                se_X[j, s_idx] = float(r["se"])
                iv_in_state[j, s_idx] = True
            except (ValueError, TypeError):
                continue
        active = iv_in_state.any(axis=0).sum()
        if active < 1:
            skip_count += 1
            continue
        try:
            fit = fit_hcsmr(
                beta_hat=beta_hat, se_X=se_X,
                gamma_hat=gamma_hat, se_Y=se_Y,
                state_to_type=state_to_type, iv_in_state=iv_in_state,
                n_warmup=args.n_warmup, n_samples=args.n_samples, n_chains=args.n_chains,
                chain_method="sequential", progress_bar=False, seed=gi,
            )
        except Exception as e:
            log.warning("fit failed for %s: %s", g, e)
            skip_count += 1
            continue

        theta = fit["samples"]["theta"]
        for ci, st in enumerate(state_list):
            if iv_in_state[:, ci].sum() == 0:
                continue
            tt = theta[:, ci]
            bf = bayes_factor_zero(tt)
            results.append({
                "gene_id": g,
                "cell_state": st,
                "parental_type": parental_list[state_to_type[ci]],
                "n_iv": int(iv_in_state[:, ci].sum()),
                "post_mean": float(tt.mean()),
                "post_sd": float(tt.std()),
                "post_q025": float(np.quantile(tt, 0.025)),
                "post_q500": float(np.quantile(tt, 0.5)),
                "post_q975": float(np.quantile(tt, 0.975)),
                "BF10": float(bf),
            })
        fit_count += 1
        if fit_count % 10 == 0:
            log.info("[%d/%d] fitted (skipped %d so far); writing partial CSV", fit_count, len(top_genes), skip_count)
            pd.DataFrame(results).to_csv(out / "pilot_results.csv", index=False)

    df = pd.DataFrame(results)
    df.to_csv(out / "pilot_results.csv", index=False)
    log.info("DONE: %d fits, %d skipped; %d (gene, state) results", fit_count, skip_count, len(df))

    if len(df):
        sig = df[(df["BF10"] > 10) & ~((df["post_q025"] < 0) & (df["post_q975"] > 0))]
        log.info("\nHits with BF10>10 and CI excluding 0: %d / %d (%.1f%%)",
                 len(sig), len(df), 100 * len(sig) / len(df))
        log.info("\nTop 20 by BF10:\n%s",
                 sig.nlargest(20, "BF10").to_string())
        # Per cell-state summary
        per_state = df.groupby("cell_state").agg(
            n_fits=("gene_id", "count"),
            n_hits=("BF10", lambda x: int((x > 10).sum())),
            mean_BF=("BF10", "mean"),
        ).reset_index().sort_values("n_hits", ascending=False)
        per_state.to_csv(out / "per_state_summary.csv", index=False)
        log.info("\nPer cell-state summary:\n%s", per_state.head(15).to_string())


if __name__ == "__main__":
    main()
