import os
import json
import sys

from processing.processing_functions import (
    downsample_adata,
    preprocess_data,
    do_PCA,
    do_neighborhood_graph,
    do_UMAP,
    do_clustering,
)
from helper_functions import cache_adata
import scanpy as sc


def main():
    session_ID = "recalc_demo"
    path = "/app/data/pbmc3k_processed.h5ad"
    out = {"session": session_ID, "steps": []}

    if not os.path.exists(path):
        print("missing input file", path)
        sys.exit(2)

    try:
        ad = sc.read_h5ad(path)
        out['steps'].append({"loaded": True, "shape": getattr(ad, 'shape', None)})
        cache_adata(session_ID, ad)
        out['steps'].append({"cached": True})

        ad = downsample_adata(session_ID, ad)
        out['steps'].append({"downsample": True})

        ad = preprocess_data(session_ID, ad, n_top_genes=300)
        out['steps'].append({"preprocess": True, "shape": getattr(ad, 'shape', None)})

        ad = do_PCA(session_ID, ad, n_comps=8)
        out['steps'].append({"pca": True, "obsm": list(ad.obsm.keys())})

        ad = do_neighborhood_graph(session_ID, ad, n_neighbors=8)
        out['steps'].append({"neighbors": True, "obsp": list(ad.obsp.keys())})

        ad = do_UMAP(session_ID, ad)
        out['steps'].append({"umap": True, "obsm": list(ad.obsm.keys())})

        ad = do_clustering(session_ID, ad)
        out['steps'].append({"clustering": True, "obs_cols": list(ad.obs.columns)[:8]})

    except Exception as e:
        out['error'] = str(e)
        print('error', e)

    # ensure tmp exists
    try:
        os.makedirs('/app/tmp', exist_ok=True)
        with open('/app/tmp/recalc_result.json', 'w') as f:
            json.dump(out, f)
    except Exception as e:
        print('failed to write result', e)

    print('done')


if __name__ == '__main__':
    main()
