import anndata as ad
import scanpy as sc
import os
import sys
import argparse

# Add src to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

# Register lightweight stubs for UI/cache modules that otherwise import the
# main `app` and initialize DB/cache. This mirrors the approach used for
# prepare_pbmc3k_hvg_local.py so this script can run outside the container.
import types
import os as _os
if 'status.status_functions' not in sys.modules:
    status_mod = types.ModuleType('status')
    status_mod.__path__ = [os.path.join(os.path.dirname(__file__), '../src/status')]
    status_funcs = types.ModuleType('status.status_functions')

    def cache_state(session_ID, state=None, key=None, val=None):
        if state is None:
            return {}
        return state

    def cache_progress(session_ID, progress=None):
        return progress

    def cache_history(session_ID, history=None):
        return [] if history is None else history

    def build_state_table(session_ID):
        return []

    def build_history_table(session_ID):
        return []

    status_funcs.cache_state = cache_state
    status_funcs.cache_progress = cache_progress
    status_funcs.cache_history = cache_history
    status_funcs.build_state_table = build_state_table
    status_funcs.build_history_table = build_history_table

    status_layouts = types.ModuleType('status.status_layouts')
    def status_layout():
        return None
    status_layouts.status_layout = status_layout

    status_components = types.ModuleType('status.status_components')
    def status_progress():
        return None
    def status_history():
        return None
    def status_state():
        return None
    status_components.status_progress = status_progress
    status_components.status_history = status_history
    status_components.status_state = status_state

    sys.modules['status'] = status_mod
    sys.modules['status.status_functions'] = status_funcs
    sys.modules['status.status_layouts'] = status_layouts
    sys.modules['status.status_components'] = status_components

# Minimal plotting stubs to avoid optional UI/visual deps during local runs
if 'plotting.multi_color_scale' not in sys.modules:
    pm = types.ModuleType('plotting')
    pm.__path__ = [os.path.join(os.path.dirname(__file__), '../src/plotting')]
    pm_sub = types.ModuleType('plotting.multi_color_scale')

    class MultiColorScale:
        def __init__(self, n_levels=None):
            self.scale_dict = {str(i): "#%02x%02x%02x" % (int(255 * i / 10), 0, 0) for i in range(11)}
            self.value_dict = {v: k for k, v in self.scale_dict.items()}
            self.scale_list = list(self.scale_dict.items())
            self.value_list = list(self.value_dict.items())

        def get_scale_list(self):
            return self.scale_list

        def get_scale_dict(self):
            return self.scale_dict

        def get_value_list(self):
            return self.value_list

        def get_value_dict(self):
            return self.value_dict

    pm_sub.MultiColorScale = MultiColorScale
    pm_params = types.ModuleType('plotting.plotting_parameters')
    pm_params.scale = 250
    pm_params.scaleratio = 1.0
    pm_params.pt_expression_scaleratio = 0.5
    pm_params.violin_expression_scaleratio = 1.5
    pm_params.margin = {"r": 50, "l": 50, "t": 50, "b": 50}
    pm_params.plot_height = 450
    pm_params.point_line_width_2d = 0.5
    pm_params.point_line_width_3d = 0.5
    pm_params.point_size_2d = 7
    pm_params.point_size_3d = 2.5
    pm_params.point_size_pt_trend = 2
    pm_params.min_opacity = 0.15
    pm_params.max_opacity = 1

    sys.modules['plotting'] = pm
    sys.modules['plotting.multi_color_scale'] = pm_sub
    sys.modules['plotting.plotting_parameters'] = pm_params


# The app dynamically imports calculate_correlations; try local first and
# fallback to BigSur package (mirrors original script behavior).
try:
    from correlations.correlation_functions import calculate_correlations
except ImportError:
    try:
        from BigSur.correlations import calculate_correlations
    except ImportError:
        calculate_correlations = None


def main():
    parser = argparse.ArgumentParser(description="Compute gene-gene correlations for pbmc3k (local).")
    parser.add_argument("--input", "-i", default="data/selected_datasets/pbmc3k_hvg.h5ad",
                        help="Input .h5ad with HVGs (default: data/selected_datasets/pbmc3k_hvg.h5ad)")
    parser.add_argument("--outdir", "-o", default="data/correlation_results",
                        help="Output directory for correlation results (default: data/correlation_results)")
    parser.add_argument("--layer", default="counts", help="Layer to use for correlations (default: counts)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    args = parser.parse_args()

    if calculate_correlations is None:
        print("ERROR: calculate_correlations is not available.")
        print("Install the BigSur package or ensure the local `src/correlations` module is present.")
        sys.exit(1)

    input_path = args.input
    output_dir = args.outdir
    layer = args.layer
    verbose = args.verbose

    print(f"Loading dataset from {input_path}...")
    if not os.path.exists(input_path):
        print(f"ERROR: Input file not found at {input_path}")
        print("Please run the prepare_pbmc3k_hvg_local.py script first.")
        sys.exit(1)

    adata = sc.read_h5ad(input_path)
    print("Dataset loaded successfully.")

    if 'highly_variable_user' not in adata.var.columns:
        print("ERROR: No highly variable genes found in the dataset ('highly_variable_user' column missing).")
        sys.exit(1)

    hvg_mask = adata.var['highly_variable_user']
    n_hvgs = int(hvg_mask.sum())

    if n_hvgs == 0:
        print("ERROR: No highly variable genes are marked as True in the 'highly_variable_user' column.")
        sys.exit(1)

    print(f"Subsetting AnnData to {n_hvgs} highly variable genes.")
    adata_hvg = adata[:, hvg_mask].copy()

    os.makedirs(output_dir, exist_ok=True)
    write_out_path = os.path.join(output_dir, '')

    print(f"Computing correlations and saving results to {output_dir}...")

    if layer not in adata_hvg.layers:
        print(f"WARNING: '{layer}' layer not found. Using .X as raw counts.")
        adata_hvg.layers[layer] = adata_hvg.X.copy()

    calculate_correlations(
        adata_hvg,
        layer=layer,
        write_out=write_out_path,
        verbose=1 if verbose else 0,
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
