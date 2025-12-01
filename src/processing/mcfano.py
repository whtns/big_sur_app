"""Processing helpers derived from the Shiny app.

This module provides functions to:
- load a 10x h5 AnnData file
- process the AnnData (filtering, normalization)
- run mcFano feature selection (via BigSur.feature_selection.mcfano_feature_selection)
- compute default and user-selected highly-variable gene masks
- run PCA/neighbors/leiden/UMAP for the chosen masks
- produce a summary dict and a UMAP matplotlib Figure (or SVG)

These are pure functions with no Shiny dependency so they can be used
from other parts of this repository.
"""

from __future__ import annotations

import os
from io import BytesIO
from typing import Dict, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import ipdb
from utils import check_csr_data_for_int_floats

# imports for BigSur feature selection
from BigSur.feature_selection import mcfano_feature_selection
import sys
from BigSur.correlations import calculate_correlations

__all__ = [
    "load_10x_h5",
    "process_adata",
    "apply_feature_selection_and_reductions",
    "get_summary",
    "fig_to_svg",
]


def load_10x_h5(path: str):
    """Load a 10x H5 file into an AnnData object.

    This mirrors the original app's use of `sc.read_10x_h5`.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    return adata


def process_adata(adata, min_genes: int = 400, min_cells: int = 3):
    """Basic preprocessing: filter cells/genes, copy counts layer, normalize and log1p."""
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    # ensure layers exists and save raw counts
    try:
        adata.layers["counts"] = adata.X.copy()
    except Exception:
        # some backed modes or sparse may need cast
        adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def _ensure_file_id(adata) -> str:
    file_id = adata.uns.get("_shiny_file_id", None)
    if file_id is None:
        file_id = str(hash(adata.obs_names[0])) if adata.n_obs > 0 else "default"
        adata.uns["_shiny_file_id"] = file_id
    return file_id


def apply_feature_selection_and_reductions(
    adata,
    mcfano_cutoff: float = 6.0,
    pvalue_cutoff: float = 0.01,
    default_mcFano_threshold: float = 4.0,
    default_pvalue_threshold: float = 0.01,
    cacache: Optional[dict] = None,
) -> Tuple[object, dict]:
    """Run mcFano feature selection and compute HV masks and reductions.

    Modifies `adata` in-place. Returns the adata and a small dict of keys used.
    `cacache` is an optional dict used to store default mask per-file-id (mirrors Shiny session cache).
    """
    if cacache is None:
        cacache = {}
    
    if getattr(adata, "raw", None) is not None:
        try:
            adata = adata.raw.to_adata()
            sc.pp.filter_cells(adata, min_genes=200)
            sc.pp.filter_genes(adata, min_cells=3)
            adata.layers['counts'] = adata.X
            adata.raw = adata.X.copy()
        except Exception as e:
            print(f"[WARN] failed to convert adata.raw to AnnData, continuing with original adata: {e}")
    else:
        print("[DEBUG] adata.raw not present; using adata as-is for preprocessing")
        sc.pp.filter_genes(adata, min_cells=3)
        if check_csr_data_for_int_floats(adata.X):
            adata.layers['counts'] = adata.X.astype(np.int32)
    ipdb.set_trace()
    # run mcFano-based feature selection (expected to populate adata.var['mc_Fano'] and 'FDR_adj_pvalue')
    mcfano_feature_selection(adata, layer="counts", min_mcfano_cutoff = mcfano_cutoff, p_val_cutoff = pvalue_cutoff)

    # default mask (one-time per file)
    file_id = _ensure_file_id(adata)
    if file_id not in cacache:
        adata.var["highly_variable_default"] = False
        default_cutoffs = np.logical_and(
            adata.var.get("mc_Fano", np.zeros(len(adata.var_names))) > default_mcFano_threshold,
            adata.var.get("FDR_adj_pvalue", np.ones(len(adata.var_names))) < default_pvalue_threshold,
        )
        adata.var.loc[default_cutoffs, "highly_variable_default"] = True
        cacache[file_id] = adata.var["highly_variable_default"].copy()
    else:
        adata.var["highly_variable_default"] = cacache[file_id].copy()

    # reductions for default if needed
    if adata.var["highly_variable_default"].sum() > 0 and "X_umap_default_cutoffs" not in adata.obsm:
        sc.pp.pca(
            adata,
            mask_var="highly_variable_default",
            key_added="X_pca_default_cutoffs",
            svd_solver="arpack",
            n_comps=50,
        )
        sc.pp.neighbors(adata, use_rep="X_pca_default_cutoffs", key_added="neighbors_default_cutoffs")
        sc.tl.leiden(adata, neighbors_key="neighbors_default_cutoffs", key_added="leiden_default_cutoffs")
        sc.tl.umap(adata, neighbors_key="neighbors_default_cutoffs", key_added="X_umap_default_cutoffs")

    # user mask and reductions
    adata.var["highly_variable_user"] = False
    user_cutoffs = np.logical_and(
        adata.var.get("mc_Fano", np.zeros(len(adata.var_names))) > mcfano_cutoff,
        adata.var.get("FDR_adj_pvalue", np.ones(len(adata.var_names))) < pvalue_cutoff,
    )
    adata.var.loc[user_cutoffs, "highly_variable_user"] = True

    if adata.var["highly_variable_user"].sum() > 0:
        sc.pp.pca(
            adata,
            mask_var="highly_variable_user",
            key_added="X_pca_user_cutoffs",
            svd_solver="arpack",
            n_comps=50,
        )
        sc.pp.neighbors(adata, use_rep="X_pca_user_cutoffs", key_added="neighbors_user_cutoffs")
        sc.tl.leiden(adata, neighbors_key="neighbors_user_cutoffs", key_added="leiden_user_cutoffs")
        sc.tl.umap(adata, neighbors_key="neighbors_user_cutoffs", key_added="X_umap_user_cutoffs")

    keys = {
        "umap_default_key": "X_umap_default_cutoffs",
        "leiden_default_key": "leiden_default_cutoffs",
        "umap_user_key": "X_umap_user_cutoffs",
        "leiden_user_key": "leiden_user_cutoffs",
    }
    return adata, keys


def get_summary(adata) -> Dict[str, int]:
    return {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_hv_default": int(adata.var.get("highly_variable_default", pd.Series(dtype=bool)).sum()),
        "n_hv_user": int(adata.var.get("highly_variable_user", pd.Series(dtype=bool)).sum()),
    }

def fig_to_svg(fig: plt.Figure) -> str:
    buf = BytesIO()
    fig.savefig(buf, format="svg")
    buf.seek(0)
    svg = buf.read().decode("utf-8")
    plt.close(fig)
    return svg

