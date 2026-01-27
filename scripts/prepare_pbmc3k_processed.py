#!/usr/bin/env python3
"""Prepare pbmc3k_processed.h5ad from pbmc3k_raw.h5ad using app's pipeline.

Usage: run inside the project root or in-container. Writes output to /app/data/selected_datasets/
"""
import os
import sys
from pathlib import Path

# ensure src is on path when running from project root
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

import scanpy as sc
from processing.processing_functions import preprocess_data, do_PCA, do_neighborhood_graph, do_UMAP, do_clustering


def find_raw():
    candidates = [
        Path(os.environ.get('SELECTED_DATASETS_PATH', '/app/data/selected_datasets')) / 'pbmc3k_raw.h5ad',
        Path('/app/data') / 'pbmc3k_raw.h5ad',
        Path('/Library/WebServer/Documents/BigSuR/selected_datasets') / 'pbmc3k_raw.h5ad',
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def main():
    raw_path = find_raw()
    if raw_path is None:
        print('ERROR: pbmc3k_raw.h5ad not found. Checked SELECTED_DATASETS_PATH, /app/data, and legacy path.')
        sys.exit(2)

    print('Loading', raw_path)
    adata = sc.read_h5ad(str(raw_path))

    session_ID = 'pbmc3k_processed_prepare'

    # Run preprocessing same as app
    adata = preprocess_data(session_ID, adata, min_cells=2, min_genes=200, max_genes=10000, target_sum=1e6, flavor='cell_ranger', n_top_genes=2000)

    # PCA, neighbors, UMAP, clustering
    adata = do_PCA(session_ID, adata, n_comps=50, random_state=0)
    adata = do_neighborhood_graph(session_ID, adata, method='standard', n_neighbors=20, random_state=0)
    adata = do_UMAP(session_ID, adata, n_dim_proj=2, random_state=0)
    adata = do_clustering(session_ID, adata, resolution=0.5, random_state=0)

    out_dir = Path(os.environ.get('SELECTED_DATASETS_PATH', '/app/data/selected_datasets'))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / 'pbmc3k_processed.h5ad'
    print('Writing processed dataset to', out_path)
    adata.write(str(out_path))

    print('Done.')


if __name__ == '__main__':
    main()
