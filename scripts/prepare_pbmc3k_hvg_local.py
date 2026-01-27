import anndata as ad
import scanpy as sc
import os
import sys
import argparse

# Add src to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

# Provide a lightweight stub for status.status_functions to avoid importing
# `app` (which initializes the DB/cache) when running this script locally.
import types
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

    # Minimal status layouts/components stubs
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

# Provide a minimal plotting.multi_color_scale stub if the real dependency
# (package 'colour') is not available when running locally.
if 'plotting.multi_color_scale' not in sys.modules:
    pm = types.ModuleType('plotting')
    pm.__path__ = [os.path.join(os.path.dirname(__file__), '../src/plotting')]
    pm_sub = types.ModuleType('plotting.multi_color_scale')

    # Minimal plotting.plotting_parameters stub
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
    sys.modules['plotting'] = pm
    sys.modules['plotting.multi_color_scale'] = pm_sub
    sys.modules['plotting.plotting_parameters'] = pm_params

from processing.mcfano import apply_feature_selection_and_reductions


def main():
    """
    Local-friendly version of prepare_pbmc3k_hvg.py.
    Accepts CLI args and uses workspace-relative paths by default.
    """
    parser = argparse.ArgumentParser(description="Recalculate HVGs for pbmc3k (local).")
    parser.add_argument("--input", "-i", default="data/selected_datasets/pbmc3k_processed.h5ad",
                        help="Input .h5ad (default: data/selected_datasets/pbmc3k_processed.h5ad)")
    parser.add_argument("--output", "-o", default="data/selected_datasets/pbmc3k_hvg.h5ad",
                        help="Output .h5ad (default: data/selected_datasets/pbmc3k_hvg.h5ad)")
    parser.add_argument("--mcfano", type=float, default=1.0, help="mcFano cutoff (default: 1.0)")
    parser.add_argument("--pvalue", type=float, default=0.05, help="p-value cutoff (default: 0.05)")
    args = parser.parse_args()

    input_path = args.input
    output_path = args.output

    mcfano_cutoff = args.mcfano
    pvalue_cutoff = args.pvalue
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

    n_hvgs = int(adata.var['highly_variable_user'].sum()) if 'highly_variable_user' in adata.var else 0
    print(f"Found {n_hvgs} highly variable genes.")

    print(f"Saving processed dataset to {output_path}...")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path)
    print("Done.")


if __name__ == "__main__":
    main()
