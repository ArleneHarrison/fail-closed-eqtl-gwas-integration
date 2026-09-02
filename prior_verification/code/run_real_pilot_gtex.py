"""Real-data hcsMR pilot using ONLY GTEx bulk-tissue eQTL (no OneK1K).

This is the fallback pilot for the case where OneK1K is unreachable (AWS Sydney
SSL blocked from CN).  We use GTEx Heart_Left_Ventricle / Heart_Atrial_Appendage
/ Artery_Coronary as three tissue-level "cell types" and demonstrate the hcsMR
framework on them.  The cell-state resolution is therefore at tissue level
rather than within-cell-type sub-state level — methodologically less compelling
but still validates the end-to-end pipeline.

For HF outcome (HERMES) we:
  1. Load HERMES sumstats
  2. Load GTEx LV / AA / CA signif variant-gene pairs
  3. For each gene with IVs in ≥ 2 tissues, fit hcsMR jointly across the 3
     tissues (here tissue = pseudo-cell-state) with partial pooling over a
     single 'cardiovascular tissue' parental type.
  4. Output: gene × tissue posterior table

Output: SERVER_ACCOUNT_ROOT/hcsmr-cvd/results/pilot_gtex_hf/hcsmr_pilot_results.csv
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
from hcsmr.model import fit_hcsmr, bayes_factor_zero


GOLD_HF_GENES = [
    "TTN", "LMNA", "MYH7", "MYH6", "NPPA", "NPPB", "ACTC1", "BAG3",
    "FLNC", "RBM20", "PLN", "JPH2", "SOS1", "TBX5", "PCSK9", "LDLR",
    "APOB", "APOE", "ANGPTL3", "ANGPTL4", "SORT1", "LPA", "IL6R",
    "ABCG5", "ABCG8", "MIA3", "GP6", "F2", "F11", "F13B", "EDN1",
    "HDAC9", "FOXF2", "PITX2", "ZFHX3", "KCNN3", "GJA1", "SCN5A",
    "SCN10A", "PRDM16", "TBX18", "ISL1", "BMP4", "HCN4", "SHOX2",
    "GRK5", "ADRB1", "ADRB2", "AGTR1", "NOS3",
]


def load_gtex_signif_variant_gene_pairs(path: str | Path, tissue_name: str) -> pd.DataFrame:
    """Load GTEx v8 signif_variant_gene_pairs file.

    Columns: variant_id, gene_id, tss_distance, ma_samples, ma_count, maf,
             pval_nominal, slope, slope_se, pval_nominal_threshold,
             min_pval_nominal, pval_beta
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    df["tissue"] = tissue_name
    # Parse variant_id like chr1_153818236_C_T_b38 → rsid lookup not available; we
    # use the GTEx variant_id directly as the SNP key.  We'll need to match to
    # GWAS by rs (later) OR by chr:pos.  For now, just retain.
    return df


def gtex_variant_to_chrpos(variant_id: str) -> tuple[str | None, int | None]:
    """Parse a GTEx v8 variant_id (chr1_12345_A_G_b38) into (chrom, pos)."""
    try:
        parts = variant_id.split("_")
        chrom = parts[0].replace("chr", "")
        pos = int(parts[1])
        return chrom, pos
    except Exception:
        return None, None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--hermes", required=True)
    p.add_argument("--gtex-lv", required=True)
    p.add_argument("--gtex-aa", default=None)
    p.add_argument("--gtex-ca", default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--n-genes", type=int, default=200)
    p.add_argument("--n-warmup", type=int, default=400)
    p.add_argument("--n-samples", type=int, default=400)
    p.add_argument("--n-chains", type=int, default=2)
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    log = logging.getLogger("pilot_gtex")
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    # --- Stage A: load HERMES ---
    gwas = load_hermes_hf(args.hermes)
    log.info("HERMES: %d rows; cols=%s", len(gwas), gwas.columns.tolist())
    # Build chr:pos key if available
    gwas["chrpos"] = gwas["chrom"].astype(str).str.replace("chr", "") + ":" + gwas["pos"].astype("Int64").astype(str)

    # --- Stage B: load GTEx for each tissue ---
    tissue_files = {"Heart_Left_Ventricle": args.gtex_lv}
    if args.gtex_aa:
        tissue_files["Heart_Atrial_Appendage"] = args.gtex_aa
    if args.gtex_ca:
        tissue_files["Artery_Coronary"] = args.gtex_ca

    gtex_all = []
    for tname, fp in tissue_files.items():
        if not Path(fp).exists():
            log.warning("missing %s file %s, skipping", tname, fp)
            continue
        gx = load_gtex_signif_variant_gene_pairs(fp, tname)
        gx["chrom"], gx["pos"] = zip(*gx["variant_id"].apply(gtex_variant_to_chrpos))
        gx["chrpos"] = gx["chrom"].astype(str) + ":" + gx["pos"].astype(str)
        # Strip Ensembl version suffix to align with HERMES gene names (HERMES may use symbols)
        gx["gene_id_stripped"] = gx["gene_id"].astype(str).str.split(".").str[0]
        gtex_all.append(gx)
        log.info("loaded %s: %d eQTLs", tname, len(gx))

    if not gtex_all:
        log.error("no GTEx data loaded — abort")
        return
    gtex_df = pd.concat(gtex_all, ignore_index=True)
    log.info("total GTEx eQTLs across tissues: %d", len(gtex_df))

    # --- Stage C: gene-pair filtering and merge ---
    # We need to align GTEx gene IDs (Ensembl, e.g. ENSG00000155657 for TTN) with
    # whatever HERMES uses.  HERMES has 'snp' column with rsids; we merge on chr:pos.
    if "chrpos" not in gwas.columns:
        log.error("HERMES has no chrpos column — gene-symbol-based merge not yet implemented")
        return

    merged = gtex_df.merge(
        gwas[["chrpos", "snp", "beta", "se", "p"]].rename(
            columns={"beta": "gwas_beta", "se": "gwas_se", "p": "gwas_p"}
        ),
        on="chrpos", how="inner",
    )
    log.info("merged GTEx with HERMES by chr:pos -> %d rows", len(merged))
    if len(merged) == 0:
        log.error("zero merged rows — chr:pos mismatch likely; check sumstats schema")
        return
    merged.to_csv(out / "merged_iv_gwas.csv.gz", index=False, compression="gzip")

    # --- Stage D: per-gene fits ---
    # Top genes by IV count (genes with most variants across all tissues)
    gene_counts = merged.groupby("gene_id_stripped").size().sort_values(ascending=False)
    top_genes = gene_counts.head(args.n_genes).index.tolist()
    log.info("top %d genes by IV count: e.g. %s", len(top_genes), top_genes[:5])

    tissue_list = list(tissue_files.keys())
    state_to_type = np.zeros(len(tissue_list), dtype=np.int32)  # one parental type

    results = []
    for gi, g in enumerate(top_genes):
        gdf = merged[merged["gene_id_stripped"] == g]
        if len(gdf) < 5:
            continue
        snps = gdf["chrpos"].unique().tolist()
        J = len(snps)
        C = len(tissue_list)
        beta_hat = np.full((J, C), np.nan, dtype=np.float32)
        se_X = np.full((J, C), np.nan, dtype=np.float32)
        gamma_hat = np.zeros(J, dtype=np.float32)
        se_Y = np.zeros(J, dtype=np.float32)
        iv_in_state = np.zeros((J, C), dtype=bool)
        for ji, snp in enumerate(snps):
            rows = gdf[gdf["chrpos"] == snp]
            r0 = rows.iloc[0]
            gamma_hat[ji] = float(r0["gwas_beta"])
            se_Y[ji] = float(r0["gwas_se"])
            for _, r in rows.iterrows():
                tissue = r["tissue"]
                if tissue in tissue_list:
                    ci = tissue_list.index(tissue)
                    beta_hat[ji, ci] = float(r["slope"])
                    se_X[ji, ci] = float(r["slope_se"])
                    iv_in_state[ji, ci] = True

        # require ≥ 1 IV in ≥ 2 tissues  OR  ≥ 3 IVs in one tissue
        active_tissues = iv_in_state.any(axis=0).sum()
        if active_tissues < 1:
            continue
        if active_tissues == 1 and iv_in_state.sum() < 3:
            continue

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

        theta = fit["samples"]["theta"]
        for ci, st in enumerate(tissue_list):
            tt = theta[:, ci]
            bf = bayes_factor_zero(tt)
            results.append({
                "gene": g,
                "tissue": st,
                "n_iv": int(iv_in_state[:, ci].sum()),
                "post_mean": float(tt.mean()),
                "post_sd": float(tt.std()),
                "post_q025": float(np.quantile(tt, 0.025)),
                "post_q500": float(np.quantile(tt, 0.5)),
                "post_q975": float(np.quantile(tt, 0.975)),
                "BF10": float(bf),
            })
        if gi % 5 == 0:
            log.info("[%d/%d] %s fitted", gi + 1, len(top_genes), g)
            pd.DataFrame(results).to_csv(out / "hcsmr_pilot_gtex_results.csv", index=False)

    df = pd.DataFrame(results)
    df.to_csv(out / "hcsmr_pilot_gtex_results.csv", index=False)
    log.info("DONE — %d (gene, tissue) results saved", len(df))

    # Quick summary
    sig = df[(df["BF10"] > 10) & ((df["post_q025"] > 0) | (df["post_q975"] < 0))]
    log.info("hits with BF10>10 and CrI excluding 0: %d", len(sig))
    log.info("top 10 hits:\n%s", sig.sort_values("BF10", ascending=False).head(10).to_string())


if __name__ == "__main__":
    main()
