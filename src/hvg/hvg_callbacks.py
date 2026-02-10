import os

import dash
from dash import html
from dash.dependencies import Input, Output, State

from helper_functions import cache_adata, cache_history, load_selected_dataset
from app import app

from processing.mcfano import apply_feature_selection_and_reductions
import traceback
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
from plotting.plotting_parameters import margin, point_size_2d, plot_height
from utils import check_csr_data_for_int_floats


def _generate_hvg_plots(adata, umap_select, hvg_max_points, keys):
    """Helper function to generate UMAP plots from an AnnData object."""
    print("[DEBUG] Inside _generate_hvg_plots function.")
    try:
        def _sample_idx(n, max_points):
            max_points = int(max_points) if max_points is not None else 10000
            if max_points <= 0 or max_points >= n:
                return np.arange(n)
            rng = np.random.default_rng(0)
            return rng.choice(n, size=max_points, replace=False)

        def get_coords_and_clusters(u_key, l_key):
            print(f"[DEBUG] Getting coords for UMAP key: {u_key}, Leiden key: {l_key}")
            if not u_key or u_key not in adata.obsm:
                print(f"[ERROR] UMAP key '{u_key}' not found in adata.obsm. Available keys: {list(adata.obsm.keys())}")
                return None, None
            coords = np.asarray(adata.obsm[u_key])
            if coords.shape[1] < 2:
                coords = np.hstack([coords, np.zeros((coords.shape[0], 1))])
            
            clusters = None
            if l_key and l_key in adata.obs.columns:
                clusters = adata.obs[l_key]
                print(f"[DEBUG] Found clusters in adata.obs['{l_key}'].")
            else:
                print(f"[WARN] Leiden key '{l_key}' not found in adata.obs. Available columns: {list(adata.obs.columns)}")

            return coords[:, :2], clusters

        default_leiden = keys.get("leiden_default_key")
        user_leiden = keys.get("leiden_user_key")
        coords_def, clusters_def = get_coords_and_clusters(keys.get("umap_default_key"), default_leiden)
        coords_user, clusters_user = get_coords_and_clusters(keys.get("umap_user_key"), user_leiden)

        fig = make_subplots(rows=1, cols=2, subplot_titles=["Default-cutoff UMAP", "User-cutoff UMAP"])

        def add_panel(col, coords, clusters, title):
            print(f"[DEBUG] Adding panel {col} ('{title}'). Coords shape: {coords.shape if coords is not None else 'None'}")
            if coords is None:
                fig.add_annotation(text="No UMAP available", xref="paper", yref="paper",
                                   x=0.25 if col == 1 else 0.75, y=0.5, showarrow=False)
                return
            idx = _sample_idx(coords.shape[0], hvg_max_points)
            df = pd.DataFrame({'x': coords[idx, 0], 'y': coords[idx, 1]})
            if clusters is not None:
                cluster_vals = clusters.values[idx]
                print(f"[DEBUG] Panel {col}: Plotting {len(df)} points with {len(pd.unique(cluster_vals))} clusters.")
                for u in pd.unique(cluster_vals):
                    mask = cluster_vals == u
                    fig.add_trace(go.Scattergl(x=df.x[mask], y=df.y[mask], mode='markers',
                                               marker=dict(size=point_size_2d), name=str(u), showlegend=(col == 1)), row=1, col=col)
            else:
                print(f"[DEBUG] Panel {col}: Plotting {len(df)} points without clusters.")
                fig.add_trace(go.Scattergl(x=df.x, y=df.y, mode='markers', marker=dict(size=point_size_2d),
                                           name=title, showlegend=False), row=1, col=col)

        add_panel(1, coords_def, clusters_def, 'Default')
        add_panel(2, coords_user, clusters_user, 'User')
        fig.update_layout(height=plot_height, showlegend=True, autosize=False)
        print("[DEBUG] Successfully generated plots.")
        return fig
    except Exception as e:
        print(f"[ERROR] Failed to generate HVG plots: {e}")
        import traceback; traceback.print_exc()
        return go.Figure()


@app.callback(
    Output('hvg-initial-load-trigger', 'data'),
    [Input('main_tabs', 'active_tab')],
    [State('session-id', 'children')]
)
def initial_hvg_load(active_tab, session_ID):
    """
    This callback automatically loads the pre-calculated pbmc3k_hvg.h5ad
    dataset into the current session when the 'Highly Variable Genes' tab
    is selected for the first time, ONLY if no other dataset is loaded.
    """
    print(f"[DEBUG] initial_hvg_load called: active_tab={active_tab}, session_ID={session_ID}")

    # Only run when the HVG tab is selected
    if active_tab != 'hvg_tab':
        return dash.no_update

    # Check if a dataset is already in the cache for this session
    adata = cache_adata(session_ID)
    if adata is not None:
        print(f"[DEBUG] A dataset is already in cache for session {session_ID}. Skipping default HVG load.")
        # Trigger recalc_hvgs with existing data, which will force a plot refresh
        return {'status': 'existing_data'}

    # If no data is cached, load the default pbmc3k_hvg dataset via the template cache
    # (fast Zarr copy) rather than parsing the raw .h5ad file on every new session.
    print(f"[{session_ID}] No dataset loaded. Loading default pre-calculated HVG dataset.")
    try:
        hvg_adata = load_selected_dataset(session_ID, "00002")  # 00002 == pbmc3k_hvg
        if hvg_adata is None:
            print(f"[{session_ID}] load_selected_dataset returned None for pbmc3k_hvg")
            return dash.no_update
        print(f"[DEBUG] Loaded default AnnData: obs={hvg_adata.n_obs}, var={hvg_adata.n_vars}")
        cache_history(session_ID, history="Loaded default pre-calculated HVG dataset (pbmc3k_hvg)")
    except Exception as e:
        print(f"[{session_ID}] Failed to load default pre-calculated HVG dataset: {e}")
        import traceback; traceback.print_exc()
        return dash.no_update

    print(f"[DEBUG] Default HVG dataset loaded and cached for session {session_ID}")
    return {'status': 'loaded_default'}


@app.callback(
    [Output('hvg_list', 'children'), Output('hvg_status', 'children'), Output('hvg_umap_graph', 'figure')],
    [Input('hvg-initial-load-trigger', 'data'), Input('recalc_hvgs_button', 'n_clicks'), Input('hvg_UMAP_dropdown', 'value'), Input('hvg_n_dims_radio', 'value'), Input('hvg_umap_select', 'value')],
    [State('session-id', 'children'), State('mcfano_cutoff', 'value'), State('pvalue_cutoff', 'value'), State('hvg_max_points', 'value'), State('default_mcfano_cutoff', 'value'), State('default_pvalue_cutoff', 'value')]
)
def recalc_hvgs(initial_load_trigger, n_clicks, processing_plot_type, n_dim_proj_plot, umap_select, session_ID, mcfano_cutoff, pvalue_cutoff, hvg_max_points, default_mcfano_cutoff, default_pvalue_cutoff):
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update, dash.no_update

    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    adata = cache_adata(session_ID)
    if adata is None:
        return "", "No dataset cached for this session", go.Figure()

    keys = {
        "umap_default_key": "X_umap_default_cutoffs",
        "leiden_default_key": "leiden_default_cutoffs",
        "umap_user_key": "X_umap_user_cutoffs",
        "leiden_user_key": "leiden_user_cutoffs",
    }

    # --- PATH 1: Initial load or UI change with pre-calculated data ---
    if trigger_id != 'recalc_hvgs_button':
        if 'highly_variable_user' in adata.var.columns and 'X_umap_user_cutoffs' in adata.obsm:
            print("[INFO] Pre-calculated HVG data found. Generating plots directly.")
            fig = _generate_hvg_plots(adata, umap_select, hvg_max_points, keys)
            hvgs = list(adata.var[adata.var['highly_variable_user']].index)
            status = f"Displaying {len(hvgs)} pre-calculated HVGs."
            hvg_list = html.Ul([html.Li(g) for g in hvgs])
            return hvg_list, status, fig

    # --- PATH 2: Explicit recalculation required ---
    print("[INFO] Recalculating HVGs based on UI controls.")
    try:
        if 'counts' not in adata.layers:
            adata.layers['counts'] = adata.X.copy()
        
        adata, returned_keys = apply_feature_selection_and_reductions(
            adata,
            mcfano_cutoff=float(mcfano_cutoff),
            pvalue_cutoff=float(pvalue_cutoff),
            default_mcFano_threshold=float(default_mcfano_cutoff),
            default_pvalue_threshold=float(default_pvalue_cutoff),
        )
        cache_adata(session_ID, adata)
        
        fig = _generate_hvg_plots(adata, umap_select, hvg_max_points, returned_keys)
        hvgs = list(adata.var[adata.var['highly_variable_user']].index)
        status = f"Found {len(hvgs)} HVGs with current cutoffs."
        hvg_list = html.Ul([html.Li(g) for g in hvgs])
        
        return hvg_list, status, fig

    except Exception as e:
        print(f"[ERROR] HVG recalculation failed: {e}")
        return "Error during recalculation.", str(e), go.Figure()
