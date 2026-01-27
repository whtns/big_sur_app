import anndata as ad
import scanpy as sc
import os
import sys

# Add src to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from processing.mcfano import apply_feature_selection_and_reductions

def main():
    """
    This script loads the pbmc3k_processed dataset, recalculates highly variable
    genes (HVGs) using the same logic as the "Highly Variable Genes" tab in the app,
    and saves the result to a new AnnData file.
    """
    input_path = "/app/data/selected_datasets/pbmc3k_processed.h5ad"
    output_path = "/app/data/selected_datasets/pbmc3k_hvg.h5ad"
    
    # Default cutoff values from hvg_callbacks.py
    mcfano_cutoff = 1.0
    pvalue_cutoff = 0.05
    default_mcfano_cutoff = 0.9
    default_pvalue_cutoff = 0.05

    print(f"Loading dataset from {input_path}...")
    if not os.path.exists(input_path):
        print(f"ERROR: Input file not found at {input_path}")
        sys.exit(1)
        
    adata = sc.read_h5ad(input_path)
    print("Dataset loaded successfully.")
    
    print(f"Recalculating highly variable genes with mcFano > {mcfano_cutoff} and p-value < {pvalue_cutoff}...")

    # Ensure raw counts layer exists for mcFano
    if 'counts' not in adata.layers:
        print("Copying .X to .layers['counts'] for mcFano processing.")
        adata.layers['counts'] = adata.X.copy()

    # Apply the same feature selection logic as the app
    adata, _ = apply_feature_selection_and_reductions(
        adata,
        mcfano_cutoff=mcfano_cutoff,
        pvalue_cutoff=pvalue_cutoff,
        default_mcFano_threshold=default_mcfano_cutoff,
        default_pvalue_threshold=default_pvalue_cutoff,
        cacache=None,
    )
    
    n_hvgs = adata.var['highly_variable_user'].sum()
    print(f"Found {n_hvgs} highly variable genes.")

    print(f"Saving processed dataset to {output_path}...")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path)
    print("Done.")

if __name__ == "__main__":
    main()
