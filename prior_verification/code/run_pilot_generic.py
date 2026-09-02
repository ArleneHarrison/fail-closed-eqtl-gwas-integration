"""Generic hcsMR pilot for any CVD outcome.

Outcomes supported:
  hf  : HERMES (zip with 'b', 'se', 'p', 'SNP' columns)
  af  : Nielsen 2018 AF (zip with 'Effect', 'StdErr', 'P-value', 'MarkerName')
  is  : GIGASTROKE (gzipped tsv with 'beta', 'standard_error', 'p_value', no rsid — uses chr:pos)

Memory-frugal: gc.collect() after every fit; smaller defaults.
"""
from __future__ import annotations

import argparse
import gc
import gzip
import logging
import os
import sys
import time
import zipfile
from pathlib import Path

os.environ.setdefault("JAX_PLATFORMS", "cpu")
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from hcsmr.model import fit_hcsmr, bayes_factor_zero

# Same cell-state table as run_pilot_credible_sets.py
CELL_STATES = [
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


def load_hf_hermes(path):
    """Returns DataFrame with rsid, beta, se, p columns."""
    with zipfile.ZipFile(path) as z:
        inner = [n for n in z.namelist() if n.endswith(".txt")][0]
        with z.open(inner) as f:
            df = pd.read_csv(f, sep="\t", low_memory=False)
    return df.rename(columns={"SNP": "rsid", "b": "beta", "se": "se", "p": "p"})[["rsid", "beta", "se", "p"]]


def load_af_nielsen(path):
    """Nielsen 2018 AF: MarkerName / chr / pos / Effect / StdErr / P-value."""
    with zipfile.ZipFile(path) as z:
        inner = [n for n in z.namelist() if n.endswith(".txt") and "README" not in n][0]
        with z.open(inner) as f:
            df = pd.read_csv(f, sep="\t", low_memory=False)
    return df.rename(columns={"MarkerName": "rsid", "Effect": "beta", "StdErr": "se", "P-value": "p"})[["rsid", "beta", "se", "p"]]


def load_is_gigastroke(path):
    """GIGASTROKE: chromosome / base_pair_location / beta / standard_error / p_value."""
    df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    df["rsid"] = df["chromosome"].astype(str) + ":" + df["base_pair_location"].astype(str)
    return df.rename(columns={"beta": "beta", "standard_error": "se", "p_value": "p"})[["rsid", "beta", "se", "p"]]


def load_cad_aragam(path):
    """CAD Aragam 2022 GCST90132314: p_value / chromosome / base_pair_location / beta / standard_error / markername (chr:pos_REF_ALT)."""
    df = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    df["rsid"] = df["chromosome"].astype(str) + ":" + df["base_pair_location"].astype(str)
    return df.rename(columns={"beta": "beta", "standard_error": "se", "p_value": "p"})[["rsid", "beta", "se", "p"]]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gwas-format", required=True, choices=["hf", "af", "is", "cad"])
    p.add_argument("--gwas-path", required=True)
    p.add_argument("--data-root", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--top-n-genes", type=int, default=60)
    p.add_argument("--min-ivs-per-state", type=int, default=3)
    p.add_argument("--max-ivs-per-gene", type=int, default=30)
    p.add_argument("--n-warmup", type=int, default=200)
    p.add_argument("--n-samples", type=int, default=200)
    p.add_argument("--n-chains", type=int, default=1)
    p.add_argument("--pip-threshold", type=float, default=0.05)
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    log = logging.getLogger("pilot")
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    data_root = Path(args.data_root)

    log.info("Loading %s GWAS from %s", args.gwas_format, args.gwas_path)
    if args.gwas_format == "hf":
        gwas = load_hf_hermes(args.gwas_path)
    elif args.gwas_format == "af":
        gwas = load_af_nielsen(args.gwas_path)
    elif args.gwas_format == "is":
        gwas = load_is_gigastroke(args.gwas_path)
    else:
        gwas = load_cad_aragam(args.gwas_path)
    gwas = gwas.rename(columns={"beta": "gwas_beta", "se": "gwas_se", "p": "gwas_p"})
    for col in ("gwas_beta", "gwas_se", "gwas_p"):
        gwas[col] = pd.to_numeric(gwas[col], errors="coerce")
    gwas = gwas.dropna(subset=["gwas_beta", "gwas_se"])
    log.info("GWAS loaded: %d rows", len(gwas))

    # Load credible sets and merge
    merged_per_state = {}
    for state_name, parental, rel in CELL_STATES:
        path = data_root / rel
        if not path.exists():
            continue
        try:
            cs = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
        except Exception as e:
            log.warning("load failed for %s: %s", state_name, e)
            continue
        cs = cs[cs["pip"] >= args.pip_threshold].copy()
        if args.gwas_format in ("is", "cad"):
            cs["chrpos"] = cs["variant"].str.replace("chr", "").str.split("_").str[:2].apply(lambda x: ":".join(x))
            cs["rsid_for_merge"] = cs["chrpos"]
        else:
            cs["rsid_for_merge"] = cs["rsid"]
        cs = cs.sort_values("pip", ascending=False).groupby(["gene_id", "cs_id"]).head(1)
        merged = cs.merge(gwas, left_on="rsid_for_merge", right_on="rsid", how="inner", suffixes=("_eqtl", "_gwas"))
        if len(merged) < args.min_ivs_per_state:
            continue
        merged["state"] = state_name
        merged["parental"] = parental
        merged_per_state[state_name] = merged
        log.info("state %-32s (%-12s): %d IVs", state_name, parental, len(merged))

    if not merged_per_state:
        log.error("no states have IVs — abort")
        return

    merged_all = pd.concat(merged_per_state.values(), ignore_index=True)
    merged_all.to_csv(out / "merged_iv_gwas.csv.gz", index=False, compression="gzip")
    log.info("merged_all: %d rows across %d states", len(merged_all), len(merged_per_state))

    # Find genes with most state coverage
    gene_state_counts = merged_all.groupby(["gene_id"])["state"].nunique().sort_values(ascending=False)
    top_genes = gene_state_counts.head(args.top_n_genes).index.tolist()

    # Build hierarchy
    state_list = sorted(merged_per_state.keys())
    parental_list = sorted(set(s[1] for s in CELL_STATES if s[0] in state_list))
    state_to_type = np.array([parental_list.index(merged_per_state[s].iloc[0]["parental"])
                               for s in state_list], dtype=np.int32)
    n_types = len(parental_list)

    # Fit per gene
    results = []
    fit_count = 0
    for gi, g in enumerate(top_genes):
        gdf = merged_all[merged_all["gene_id"] == g]
        if len(gdf) < 3:
            continue
        if len(gdf) > args.max_ivs_per_gene:
            gdf = gdf.sort_values("pip", ascending=False).head(args.max_ivs_per_gene)
        snps = gdf["rsid_for_merge"].unique().tolist()
        J = len(snps); C = len(state_list)
        # The "beta" column from credible_sets clashed with our renaming; use "beta_eqtl"
        beta_col = "beta_eqtl" if "beta_eqtl" in gdf.columns else "beta"
        se_col = "se_eqtl" if "se_eqtl" in gdf.columns else "se"
        beta_hat = np.full((J, C), np.nan, dtype=np.float32)
        se_X = np.full((J, C), np.nan, dtype=np.float32)
        gamma_hat = np.zeros(J, dtype=np.float32)
        se_Y = np.zeros(J, dtype=np.float32)
        iv_in_state = np.zeros((J, C), dtype=bool)
        rsid_to_j = {snp: j for j, snp in enumerate(snps)}
        for _, r in gdf.iterrows():
            j = rsid_to_j[r["rsid_for_merge"]]
            gamma_hat[j] = float(r["gwas_beta"])
            se_Y[j] = float(r["gwas_se"])
            s_idx = state_list.index(r["state"])
            try:
                beta_hat[j, s_idx] = float(r[beta_col])
                se_X[j, s_idx] = float(r[se_col])
                iv_in_state[j, s_idx] = True
            except (ValueError, TypeError, KeyError):
                continue
        active = iv_in_state.any(axis=0).sum()
        if active < 1:
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
            log.warning("fit %d/%d failed for %s: %s", gi + 1, len(top_genes), g, e)
            gc.collect()
            continue

        theta = fit["samples"]["theta"]
        for ci, st in enumerate(state_list):
            if iv_in_state[:, ci].sum() == 0:
                continue
            tt = theta[:, ci]
            try:
                bf = bayes_factor_zero(tt)
            except Exception:
                bf = float("nan")
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
        # Memory hygiene
        del fit, theta
        gc.collect()
        if fit_count % 5 == 0:
            log.info("[%d/%d] %d fits done", fit_count, len(top_genes), fit_count)
            pd.DataFrame(results).to_csv(out / "pilot_results.csv", index=False)

    df = pd.DataFrame(results)
    df.to_csv(out / "pilot_results.csv", index=False)
    log.info("DONE: %d fits, %d (gene,state) results", fit_count, len(df))


if __name__ == "__main__":
    main()
