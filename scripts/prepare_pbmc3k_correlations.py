import anndata as ad
import scanpy as sc
import os
import sys

# Add src to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

# The app dynamically imports calculate_correlations, we'll do the same.
try:
    from correlations.correlation_functions import calculate_correlations
except ImportError:
    try:
        from BigSur.correlations import calculate_correlations
    except ImportError:
        print("ERROR: calculate_correlations is not available.")
        print("Install the BigSur package or ensure the local `src/correlations` module is present.")
        sys.exit(1)

def main():
    """
    This script loads the pbmc3k_hvg dataset, computes gene-gene correlations
    for the highly variable genes, and saves the results to disk, similar to
    the 'Correlations' tab in the app.
    """
    input_path = "/app/data/selected_datasets/pbmc3k_hvg.h5ad"
    output_dir = "/app/data/correlation_results"

    print(f"Loading dataset from {input_path}...")
    if not os.path.exists(input_path):
        print(f"ERROR: Input file not found at {input_path}")
        print("Please run the `prepare_pbmc3k_hvg.py` script first.")
        sys.exit(1)
        
    adata = sc.read_h5ad(input_path)
    print("Dataset loaded successfully.")

    # Check for HVG genes
    if 'highly_variable_user' not in adata.var.columns:
        print("ERROR: No highly variable genes found in the dataset ('highly_variable_user' column missing).")
        sys.exit(1)
    
    hvg_mask = adata.var['highly_variable_user']
    n_hvgs = hvg_mask.sum()

    if n_hvgs == 0:
        print("ERROR: No highly variable genes are marked as True in the 'highly_variable_user' column.")
        sys.exit(1)

    print(f"Subsetting AnnData to {n_hvgs} highly variable genes.")
    adata_hvg = adata[:, hvg_mask].copy()

    os.makedirs(output_dir, exist_ok=True)
    # The calculate_correlations function expects the path to end with a slash
    write_out_path = os.path.join(output_dir, '')

    print(f"Computing correlations and saving results to {output_dir}...")
    
    # Ensure a 'counts' layer exists, as required by the correlation function
    if 'counts' not in adata_hvg.layers:
        print("WARNING: 'counts' layer not found. Using .X as raw counts.")
        adata_hvg.layers['counts'] = adata_hvg.X.copy()

    calculate_correlations(
        adata_hvg,
        layer='counts',
        write_out=write_out_path,
        verbose=1
    )
    
    mcPCCs_path = os.path.join(output_dir, 'mcPCCs.npz')
    pvalues_path = os.path.join(output_dir, 'BH_corrected_pvalues.npz')

    if os.path.exists(mcPCCs_path) and os.path.exists(pvalues_path):
        print("Correlation calculation complete.")
        print(f"  - mcPCCs matrix saved to: {mcPCCs_path}")
        print(f"  - P-values saved to: {pvalues_path}")
    else:
        print("ERROR: Correlation calculation finished, but output files were not found.")
        sys.exit(1)

    print("Done.")

if __name__ == "__main__":
    main()
