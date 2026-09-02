"""Atlas integration and cell-state hierarchy construction for hcsMR.

Loads the Heart Cell Atlas v2 (and optionally Tabula Sapiens Heart), QCs cells,
performs batch-aware integration with scVI (preferred) or Harmony (fallback)
when scvi-tools is available, then builds a multi-resolution cell-state
hierarchy via Leiden clustering and consensus pooling.

Output: an integrated `anndata` object with a `cell_state` column at multiple
resolutions, suitable as input to the IV library construction step.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s [%(levelname)s] %(message)s",
    )


def load_hca_global(path: str | Path):
    """Load the HCA-Heart v2 Global_raw.h5ad file (8.5 GB on disk)."""
    import anndata as ad
    logger.info("loading HCA-Heart from %s", path)
    a = ad.read_h5ad(path, backed=None)
    logger.info("loaded: %d cells x %d genes", a.n_obs, a.n_vars)
    return a


def qc_minimal(adata, *, min_genes: int = 200, max_pct_mito: float = 10.0):
    """Apply a conservative QC filter to a raw atlas.

    HCA-Heart already has its own QC done; this is a defensive double-check.
    """
    import scanpy as sc
    logger.info("running QC: min_genes=%d max_pct_mito=%.1f", min_genes, max_pct_mito)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=10)
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    n_before = adata.n_obs
    if "pct_counts_mt" in adata.obs.columns:
        adata = adata[adata.obs["pct_counts_mt"] < max_pct_mito].copy()
    logger.info("QC: kept %d/%d cells", adata.n_obs, n_before)
    return adata


def normalise_and_subset(adata, *, n_hvg: int = 4000, subset_n_cells: Optional[int] = None,
                         random_state: int = 0):
    """Normalise to 1e4 total, log1p, find HVGs, optionally subset for speed."""
    import scanpy as sc

    if subset_n_cells is not None and adata.n_obs > subset_n_cells:
        logger.info("subsetting to %d cells for development speed", subset_n_cells)
        rng = np.random.default_rng(random_state)
        idx = rng.choice(adata.n_obs, size=subset_n_cells, replace=False)
        adata = adata[np.sort(idx)].copy()

    # Save raw counts before normalisation (needed by scVI)
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, batch_key=_batch_key(adata))
    return adata


def _batch_key(adata) -> Optional[str]:
    """Try to find a sensible batch column in adata.obs."""
    for candidate in ("donor", "donor_id", "sample", "sample_id", "batch", "study", "dataset", "Sample"):
        if candidate in adata.obs.columns:
            return candidate
    return None


def integrate_harmony(adata, *, n_pcs: int = 30, batch_key: Optional[str] = None):
    """Batch-aware integration via Harmony (CPU-friendly fallback to scVI)."""
    import scanpy as sc

    if batch_key is None:
        batch_key = _batch_key(adata)
    if batch_key is None:
        logger.warning("no batch key found; skipping integration")
        sc.pp.pca(adata, n_comps=n_pcs, mask_var="highly_variable")
        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15)
        return adata

    logger.info("integrating via Harmony on batch=%s", batch_key)
    sc.pp.pca(adata, n_comps=n_pcs, mask_var="highly_variable")
    sc.external.pp.harmony_integrate(adata, key=batch_key)
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15)
    return adata


def integrate_scvi(adata, *, batch_key: Optional[str] = None, n_latent: int = 30,
                    max_epochs: int = 100, use_gpu: bool = False):
    """Batch-aware integration via scVI."""
    import scvi
    if batch_key is None:
        batch_key = _batch_key(adata)
    if batch_key is None:
        raise ValueError("scVI requires a batch_key")
    logger.info("integrating via scVI on batch=%s n_latent=%d max_epochs=%d use_gpu=%s",
                batch_key, n_latent, max_epochs, use_gpu)
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key, layer="counts")
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(max_epochs=max_epochs, accelerator="cpu" if not use_gpu else "gpu")
    adata.obsm["X_scvi"] = model.get_latent_representation()
    import scanpy as sc
    sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=15)
    return adata, model


def multi_resolution_leiden(adata, *, resolutions=(0.3, 0.5, 1.0, 1.5), random_state: int = 0):
    """Multi-resolution Leiden clustering → cell-state labels at each resolution."""
    import scanpy as sc
    for r in resolutions:
        key = f"cell_state_res{r}"
        logger.info("running Leiden at resolution %.2f", r)
        sc.tl.leiden(adata, resolution=r, key_added=key, random_state=random_state)
    return adata


def annotate_cell_types(adata, *, marker_gene_table: Optional[dict[str, list[str]]] = None,
                         res_key: str = "cell_state_res0.5"):
    """Assign a parental cell-type label to each cell state by marker enrichment.

    `marker_gene_table` is a {cell_type: [marker genes]} mapping.  If None, use
    a default cardiac panel.
    """
    import scanpy as sc

    if marker_gene_table is None:
        marker_gene_table = DEFAULT_CARDIAC_MARKERS

    # Score each cell against each cell-type panel
    for ct, markers in marker_gene_table.items():
        valid = [g for g in markers if g in adata.var_names]
        if not valid:
            continue
        sc.tl.score_genes(adata, valid, score_name=f"score_{ct}")

    # For each Leiden cluster, take the highest-scoring cell-type
    if res_key not in adata.obs.columns:
        raise KeyError(f"resolution key {res_key} not in obs")
    cluster_to_type: dict[str, str] = {}
    for cl in adata.obs[res_key].unique():
        mask = adata.obs[res_key] == cl
        means = {ct: float(adata.obs.loc[mask, f"score_{ct}"].mean())
                 for ct in marker_gene_table if f"score_{ct}" in adata.obs.columns}
        if not means:
            cluster_to_type[cl] = "unknown"
            continue
        cluster_to_type[cl] = max(means, key=means.get)
    adata.obs["cell_type"] = adata.obs[res_key].map(cluster_to_type).astype("category")
    return adata, cluster_to_type


DEFAULT_CARDIAC_MARKERS = {
    "cardiomyocyte": ["TNNT2", "TNNI3", "MYH7", "MYH6", "MYL2", "MYL7", "ACTC1", "NPPA", "NPPB", "PLN", "RYR2", "TTN"],
    "stressed_cardiomyocyte": ["NPPA", "NPPB", "ANKRD1", "MYH7", "ACTA2", "FOS", "JUN", "ATF3"],
    "atrial_cardiomyocyte": ["NPPA", "MYL7", "KCNA5", "HEY1"],
    "ventricular_cardiomyocyte": ["MYL2", "MYH7", "IRX4"],
    "pacemaker_cardiomyocyte": ["HCN4", "SHOX2", "ISL1", "TBX18", "BMP4"],
    "fibroblast": ["DCN", "COL1A1", "COL1A2", "COL3A1", "PDGFRA", "POSTN", "FAP"],
    "activated_fibroblast": ["POSTN", "FAP", "ACTA2", "TAGLN", "CTHRC1", "COMP"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "EGFL7", "TIE1", "KDR"],
    "smooth_muscle": ["ACTA2", "MYH11", "TAGLN", "CNN1", "MYL9", "MYLK"],
    "pericyte": ["RGS5", "PDGFRB", "ACTA2", "MCAM", "ABCC9", "KCNJ8"],
    "macrophage": ["CD68", "CD163", "C1QA", "C1QB", "MRC1", "MARCO", "LYVE1"],
    "foamcell_macrophage": ["TREM2", "OLR1", "PLIN2", "ABCA1", "ABCG1", "FABP4"],
    "monocyte": ["CD14", "FCN1", "VCAN", "S100A8", "S100A9", "LYZ"],
    "T_cell": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC", "CD4", "CD8A"],
    "B_cell": ["CD19", "MS4A1", "CD79A", "CD79B"],
    "NK_cell": ["NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"],
    "mast_cell": ["KIT", "TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "adipocyte": ["FABP4", "ADIPOQ", "LEP", "PLIN1", "PLIN4"],
    "neuronal": ["PLP1", "MPZ", "S100B", "SOX10", "NRXN1"],
    "lymphatic_endothelial": ["PROX1", "LYVE1", "PDPN", "FLT4"],
}


def build_cell_state_hierarchy(adata, *, res_main: float = 0.5, sub_res: float = 1.5):
    """Build a two-level cell-state hierarchy: cell type (broad) × cell state (fine).

    The 'cell_type' column is already assigned by annotate_cell_types().  We
    then re-cluster within each cell type at higher resolution to get cell
    states nested under their parental cell type.
    """
    import scanpy as sc

    if "cell_type" not in adata.obs.columns:
        raise KeyError("annotate_cell_types() must be called before build_cell_state_hierarchy()")

    cell_state = adata.obs["cell_type"].astype(str).copy()
    for ct in adata.obs["cell_type"].cat.categories:
        idx = adata.obs.index[adata.obs["cell_type"] == ct]
        if len(idx) < 200:
            cell_state.loc[idx] = f"{ct}__0"
            continue
        sub = adata[idx].copy()
        sc.pp.neighbors(sub, use_rep=("X_scvi" if "X_scvi" in adata.obsm else
                                       "X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca"),
                        n_neighbors=10)
        sc.tl.leiden(sub, resolution=sub_res, key_added="sub")
        for sublabel in sub.obs["sub"].unique():
            sub_idx = sub.obs.index[sub.obs["sub"] == sublabel]
            cell_state.loc[sub_idx] = f"{ct}__{sublabel}"
    adata.obs["cell_state"] = cell_state.astype("category")

    # Build state→type mapping
    state_to_type: dict[str, str] = {}
    for s in adata.obs["cell_state"].cat.categories:
        parental = s.split("__")[0]
        state_to_type[s] = parental
    return adata, state_to_type
