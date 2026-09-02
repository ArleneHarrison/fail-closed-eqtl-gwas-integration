"""Real-data hcsMR pilot: HF (HERMES) using OneK1K immune-cell IVs + GTEx cardiac bulk fallback.

The pilot is intentionally minimal in scope:
  - Load HERMES sumstats
  - Load OneK1K combined eQTL
  - Define a few cardiac immune cell states (macrophage, monocyte, T-cell, NK)
  - Build IV library at gene × cell-state granularity
  - Merge with GWAS
  - Fit hcsMR for each of a small set of cardiac-relevant genes (e.g., 30
    "gold standard" CVD genes + 70 highest-IV-count genes)
  - Write per-(gene, cell_state) posterior table

The goal of this pilot is to validate the end-to-end pipeline on REAL data
before scaling to the full 4-outcome × ~30 cell-state run.  Expected runtime:
~1-3 hours on 112 CPU cores with chain_method="parallel".
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.setdefault("JAX_PLATFORMS", "cpu")
os.environ.setdefault("XLA_FLAGS", "--xla_force_host_platform_device_count=8")
os.environ.setdefault("OMP_NUM_THREADS", "4")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from hcsmr.gwas_io import load_hermes_hf
from hcsmr.iv_library import (
    load_onek1k,
    STATE_TO_SCEQTL_CELLTYPE,
)


GOLD_HF_GENES = [
    "TTN", "LMNA", "MYH7", "MYH6", "NPPA", "NPPB", "ACTC1", "BAG3",
    "FLNC", "RBM20", "PLN", "JPH2", "SOS1", "TBX5", "PCSK9", "LDLR",
    "APOB", "APOE", "ANGPTL3", "ANGPTL4", "SORT1", "LPA", "IL6R",
    "ABCG5", "ABCG8", "MIA3", "GP6", "F2", "F11", "F13B", "EDN1",
    "HDAC9", "FOXF2", "PITX2", "ZFHX3", "KCNN3", "GJA1", "SCN5A",
    "SCN10A", "PRDM16", "TBX18", "ISL1", "BMP4", "HCN4", "SHOX2",
    "GRK5", "ADRB1", "ADRB2", "AGTR1", "NOS3",
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--hermes", required=True, help="HERMES HF zip or txt")
    p.add_argument("--onek1k", required=True, help="OneK1K combined eqtl tsv.gz")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--n-genes", type=int, default=100,
                   help="number of genes to fit (top by IV count, plus gold standard)")
    p.add_argument("--n-warmup", type=int, default=500)
    p.add_argument("--n-samples", type=int, default=500)
    p.add_argument("--n-chains", type=int, default=2)
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    log = logging.getLogger("pilot")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # --- Stage A: load sumstats ---
    gwas = load_hermes_hf(args.hermes)
    log.info("HERMES loaded: %d rows; columns: %s", len(gwas), gwas.columns.tolist())

    # --- Stage B: load OneK1K sc-eQTL ---
    sceqtl = load_onek1k(args.onek1k)
    log.info("OneK1K loaded: %d rows; columns: %s", len(sceqtl), sceqtl.columns.tolist())

    # --- Stage C: define pilot cell-state hierarchy (5 immune states) ---
    # For the pilot we use OneK1K cell type labels directly as 'cell states',
    # mapped to the cardiac immune parental types.
    pilot_states = {
        "Mono_C":   "monocyte",
        "Mono_NC":  "monocyte",
        "DC":       "macrophage",
        "NK":       "NK_cell",
        "CD4_NC":   "T_cell",
        "CD8_NC":   "T_cell",
        "B_Mem":    "B_cell",
    }
    log.info("pilot cell states: %s", list(pilot_states.keys()))

    # --- Stage D: gene selection ---
    if "gene" not in sceqtl.columns:
        log.error("OneK1K has no 'gene' column — abort")
        return
    gene_counts = sceqtl.groupby("gene").size().sort_values(ascending=False)
    gold_present = [g for g in GOLD_HF_GENES if g in gene_counts.index]
    log.info("gold-standard HF genes with sc-eQTL: %d/%d", len(gold_present), len(GOLD_HF_GENES))
    top_by_count = gene_counts.head(args.n_genes).index.tolist()
    pilot_genes = sorted(set(gold_present + top_by_count))
    log.info("pilot will fit %d genes", len(pilot_genes))

    # --- Stage E: build IV → GWAS-merged table per (gene, cell_state) ---
    # OneK1K-style columns assumed: snp, gene, cell_type, beta, se
    needed_cols = {"snp", "gene", "cell_type", "beta", "se"}
    missing = needed_cols - set(sceqtl.columns)
    if missing:
        log.error("OneK1K missing columns %s; cannot continue", missing)
        return

    # Restrict to pilot cell types and pilot genes
    sceqtl_p = sceqtl[sceqtl["cell_type"].isin(pilot_states) & sceqtl["gene"].isin(pilot_genes)].copy()
    log.info("after pilot restriction: %d rows", len(sceqtl_p))

    # Merge with GWAS on snp
    merged = sceqtl_p.merge(
        gwas[["snp", "beta", "se", "p"]].rename(
            columns={"beta": "gwas_beta", "se": "gwas_se", "p": "gwas_p"}
        ),
        on="snp", how="inner",
    )
    log.info("merged with HERMES: %d rows", len(merged))
    merged.to_csv(out / "merged_iv_gwas.csv.gz", index=False, compression="gzip")

    # --- Stage F: prepare per-gene array inputs for hcsMR ---
    # For each gene g, we build:
    #   beta_hat[J, C]: sc-eQTL effect per SNP per cell-state
    #   se_X[J, C]:     SE
    #   gamma_hat[J]:   GWAS effect
    #   se_Y[J]:        GWAS SE
    #   iv_in_state[J, C]: SNP j is an IV for state c iff sc-eQTL exists for (j, c)
    # In OneK1K cell_type = cell_state mapping, each SNP is naturally tied to
    # one OneK1K cell type ⇒ iv_in_state has a single True per row.

    from hcsmr.model import fit_hcsmr, bayes_factor_zero

    state_list = list(pilot_states.keys())
    state_to_type = np.array([sorted(set(pilot_states.values())).index(pilot_states[s]) for s in state_list],
                              dtype=np.int32)
    n_types = int(state_to_type.max()) + 1
    log.info("state_to_type: %s; n_types=%d", state_to_type, n_types)

    results = []
    for gi, g in enumerate(pilot_genes):
        gdf = merged[merged["gene"] == g]
        if len(gdf) < 5:
            continue
        snps = gdf["snp"].unique().tolist()
        J = len(snps)
        C = len(state_list)
        beta_hat = np.full((J, C), np.nan, dtype=np.float32)
        se_X = np.full((J, C), np.nan, dtype=np.float32)
        gamma_hat = np.zeros(J, dtype=np.float32)
        se_Y = np.zeros(J, dtype=np.float32)
        iv_in_state = np.zeros((J, C), dtype=bool)
        for ji, snp in enumerate(snps):
            rows = gdf[gdf["snp"] == snp]
            # GWAS values are constant across rows for the same SNP
            r0 = rows.iloc[0]
            gamma_hat[ji] = float(r0["gwas_beta"])
            se_Y[ji] = float(r0["gwas_se"])
            for _, r in rows.iterrows():
                ct = r["cell_type"]
                if ct in state_list:
                    ci = state_list.index(ct)
                    beta_hat[ji, ci] = float(r["beta"])
                    se_X[ji, ci] = float(r["se"])
                    iv_in_state[ji, ci] = True

        try:
            fit = fit_hcsmr(
                beta_hat=beta_hat,
                se_X=se_X,
                gamma_hat=gamma_hat,
                se_Y=se_Y,
                state_to_type=state_to_type,
                iv_in_state=iv_in_state,
                n_warmup=args.n_warmup,
                n_samples=args.n_samples,
                n_chains=args.n_chains,
                chain_method="sequential",
                progress_bar=False,
                seed=gi,
            )
        except Exception as e:
            log.warning("fit failed for %s: %s", g, e)
            continue

        theta = fit["samples"]["theta"]  # (S, C)
        for ci, st in enumerate(state_list):
            tt = theta[:, ci]
            bf = bayes_factor_zero(tt)
            results.append({
                "gene": g,
                "cell_state": st,
                "parental_cell_type": pilot_states[st],
                "n_iv": int(iv_in_state[:, ci].sum()),
                "post_mean": float(tt.mean()),
                "post_sd": float(tt.std()),
                "post_q025": float(np.quantile(tt, 0.025)),
                "post_q500": float(np.quantile(tt, 0.5)),
                "post_q975": float(np.quantile(tt, 0.975)),
                "BF10": float(bf),
                "is_gold_standard": int(g in GOLD_HF_GENES),
            })
        if gi % 10 == 0:
            log.info("[%d/%d] fitted gene=%s", gi + 1, len(pilot_genes), g)
            pd.DataFrame(results).to_csv(out / "hcsmr_pilot_results.csv", index=False)

    df = pd.DataFrame(results)
    df.to_csv(out / "hcsmr_pilot_results.csv", index=False)
    log.info("DONE — %d (gene, state) results to %s", len(df),
             out / "hcsmr_pilot_results.csv")


if __name__ == "__main__":
    main()
