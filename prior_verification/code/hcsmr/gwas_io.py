"""GWAS sumstats loaders / harmonisers for hcsMR."""
from __future__ import annotations

import logging
import zipfile
from pathlib import Path
import io

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def load_hermes_hf(path: str | Path) -> pd.DataFrame:
    """Load HERMES Jan 2019 HF GWAS (Shah/Aragam 2020, in CVDKP-distributed zip)."""
    path = Path(path)
    logger.info("loading HERMES HF from %s", path)
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as z:
            inner = [n for n in z.namelist() if n.endswith(".txt")][0]
            with z.open(inner) as f:
                df = pd.read_csv(f, sep="\t", low_memory=False)
    else:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    logger.info("HERMES columns: %s; rows: %d", df.columns.tolist(), len(df))
    return _normalise_gwas(df, outcome="HF")


def load_gigastroke(path: str | Path) -> pd.DataFrame:
    """Load GIGASTROKE (Mishra 2022) any-stroke EUR GWAS — GCST90104540 EBI FTP."""
    path = Path(path)
    logger.info("loading GIGASTROKE from %s", path)
    df = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    logger.info("GIGASTROKE columns: %s; rows: %d", df.columns.tolist(), len(df))
    return _normalise_gwas(df, outcome="IS")


def load_aragam_cad(path: str | Path) -> pd.DataFrame:
    """Load Aragam 2022 CAD GWAS (GCST90132314, 3GB tsv)."""
    path = Path(path)
    logger.info("loading CAD Aragam 2022 from %s", path)
    df = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    logger.info("CAD columns: %s; rows: %d", df.columns.tolist(), len(df))
    return _normalise_gwas(df, outcome="CAD")


def load_nielsen_af(path: str | Path) -> pd.DataFrame:
    """Load Nielsen 2018 AF GWAS (broadinstitute AF_HRC_GWAS_ALLv11 zip)."""
    path = Path(path)
    logger.info("loading Nielsen AF from %s", path)
    with zipfile.ZipFile(path) as z:
        cands = [n for n in z.namelist() if n.lower().endswith((".txt", ".tsv", ".csv"))]
        if not cands:
            cands = [n for n in z.namelist() if not n.endswith("/")]
        if not cands:
            raise FileNotFoundError(f"no usable file in {path}")
        with z.open(cands[0]) as f:
            df = pd.read_csv(f, sep="\t", low_memory=False)
    logger.info("Nielsen AF columns: %s; rows: %d", df.columns.tolist(), len(df))
    return _normalise_gwas(df, outcome="AF")


def _normalise_gwas(df: pd.DataFrame, *, outcome: str) -> pd.DataFrame:
    """Standardise common GWAS column names → snp, chrom, pos, ea, oa, eaf, beta, se, p."""
    rename = {}
    for c in df.columns:
        cl = c.lower()
        if cl in ("snp", "snp_id", "rsid", "rs_id", "variant_id", "marker_name", "rs"):
            rename[c] = "snp"
        elif cl in ("chr", "chrom", "chromosome", "#chr"):
            rename[c] = "chrom"
        elif cl in ("pos", "position", "bp", "base_pair_location"):
            rename[c] = "pos"
        elif cl in ("effect_allele", "ea", "a1", "allele1", "alt"):
            rename[c] = "ea"
        elif cl in ("other_allele", "oa", "a2", "allele2", "ref"):
            rename[c] = "oa"
        elif cl in ("eaf", "freq", "frq", "maf", "a1freq", "effect_allele_frequency"):
            rename[c] = "eaf"
        elif cl in ("beta", "effect", "log_or", "log_odds", "or_beta"):
            rename[c] = "beta"
        elif cl in ("se", "stderr", "std_err", "standard_error"):
            rename[c] = "se"
        elif cl in ("p", "pval", "p_value", "pvalue", "p_bolt_lmm", "p_score"):
            rename[c] = "p"
        elif cl in ("n", "n_total", "n_samples"):
            rename[c] = "n"
    df = df.rename(columns=rename)
    df["outcome"] = outcome
    # Coerce numeric
    for col in ("beta", "se", "p", "eaf", "pos", "n"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def merge_iv_with_gwas(iv_df: pd.DataFrame, gwas_df: pd.DataFrame) -> pd.DataFrame:
    """Merge IV table (sc-eQTL side) with GWAS sumstats (outcome side) on SNP."""
    if "snp" not in iv_df.columns or "snp" not in gwas_df.columns:
        raise KeyError("both IV and GWAS must have a 'snp' column for merge")
    g = gwas_df[["snp", "beta", "se", "p"]].rename(columns={"beta": "gwas_beta",
                                                              "se": "gwas_se",
                                                              "p": "gwas_p"})
    merged = iv_df.merge(g, on="snp", how="inner")
    logger.info("merged %d IVs with %d GWAS variants -> %d joint rows",
                len(iv_df), len(gwas_df), len(merged))
    return merged
