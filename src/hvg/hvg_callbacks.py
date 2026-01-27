import anndata as ad
import os

import dash
from dash import html
from dash.dependencies import Input, Output, State

from helper_functions import cache_adata, cache_history
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

    # If no data is cached, load the default pbmc3k_hvg dataset
    hvg_adata_path = "/app/data/selected_datasets/pbmc3k_hvg.h5ad"
    print(f"[DEBUG] No dataset in cache. Checking for default HVG file at {hvg_adata_path}")
    if not os.path.exists(hvg_adata_path):
        print(f"[DEBUG] Default HVG file not found at {hvg_adata_path}")
        return dash.no_update

    print(f"[{session_ID}] No dataset loaded. Loading default pre-calculated HVG dataset from {hvg_adata_path}")
    try:
        hvg_adata = ad.read_h5ad(hvg_adata_path)
        print(f"[DEBUG] Loaded default AnnData: obs={hvg_adata.n_obs}, var={hvg_adata.n_vars}")
        cache_adata(session_ID, hvg_adata)
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
    print(f"[DEBUG] recalc_hvgs called: triggered={ctx.triggered}, session_ID={session_ID}")
    if not ctx.triggered:
        print("[DEBUG] No trigger, returning no update.")
        return dash.no_update, dash.no_update, dash.no_update

    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    adata = cache_adata(session_ID)
    print(f"[DEBUG] cache_adata in recalc_hvgs: {type(adata)}")
    if adata is None:
        print("[DEBUG] No AnnData in cache for session.")
        return "", "No dataset cached for this session", {}

    # Define a default keys dictionary to use if we skip recalculation
    keys = {
        "umap_default_key": "X_umap_default_cutoffs",
        "leiden_default_key": "leiden_default_cutoffs",
        "umap_user_key": "X_umap_user_cutoffs",
        "leiden_user_key": "leiden_user_cutoffs",
    }

    try:
        # Check if the trigger is the initial load and if data is pre-calculated
        recalculate = True
        if trigger_id == 'hvg-initial-load-trigger':
            required_var = 'highly_variable_user' in adata.var.columns
            required_obsm = 'X_umap_user_cutoffs' in adata.obsm
            if required_var and required_obsm:
                print("[DEBUG] Pre-calculated HVG data found in AnnData. Skipping recalculation.")
                recalculate = False

        if recalculate:
            print("[DEBUG] Recalculating HVGs based on UI controls.")
            mcfano_cutoff = float(mcfano_cutoff) if mcfano_cutoff is not None else 1.0
            pvalue_cutoff = float(pvalue_cutoff) if pvalue_cutoff is not None else 0.05
            print(f"[DEBUG] Using mcfano_cutoff={mcfano_cutoff}, pvalue_cutoff={pvalue_cutoff}")

            # Ensure raw counts layer exists for mcFano
            if 'counts' not in getattr(adata, 'layers', {}):
                print("[DEBUG] 'counts' layer missing, copying .X to .layers['counts']")
                adata.layers['counts'] = adata.X.copy()

            # Run mcFano selection and reductions
            default_mcFano = float(default_mcfano_cutoff) if default_mcfano_cutoff is not None else 0.9
            default_pval = float(default_pvalue_cutoff) if default_pvalue_cutoff is not None else 0.05
            adata, keys = apply_feature_selection_and_reductions(
                adata,
                mcfano_cutoff=mcfano_cutoff,
                pvalue_cutoff=pvalue_cutoff,
                default_mcFano_threshold=default_mcFano,
                default_pvalue_threshold=default_pval,
                cacache=None,
            )
            print(f"[DEBUG] apply_feature_selection_and_reductions returned keys: {keys}")
            cache_history(session_ID, history=(f"Recalculated HVGs with mcFano={mcfano_cutoff}, p={pvalue_cutoff}"))
            cache_adata(session_ID, adata)
        
        # --- Plotting logic (runs for both pre-calculated and recalculated data) ---

        hvgs = []
        if 'highly_variable_user' in adata.var.columns:
            hvgs = list(adata.var[adata.var['highly_variable_user']].index)
        print(f"[DEBUG] Found {len(hvgs)} HVGs to display.")

        hvg_list = html.Ul([html.Li(g) for g in hvgs]) if hvgs else html.Div("No HVGs found.")

        default_key = keys.get('umap_default_key')
        user_key = keys.get('umap_user_key')
        chosen_key = user_key if umap_select == 'user' else default_key
        print(f"[DEBUG] UMAP keys for plotting: default={default_key}, user={user_key}, chosen={chosen_key}")

        try:
            def _sample_idx(n, max_points):
                try:
                    max_points = int(max_points) if max_points is not None else 10000
                except Exception:
                    max_points = 10000
                if max_points <= 0 or max_points >= n:
                    return np.arange(n)
                rng = np.random.default_rng(0)
                return rng.choice(n, size=max_points, replace=False)
            def get_coords_and_clusters(u_key, l_key):
                if (u_key is None) or (u_key not in adata.obsm):
                    print(f"[DEBUG] UMAP key {u_key} not in adata.obsm")
                    return None, None
                coords = np.asarray(adata.obsm[u_key])
                if coords.ndim == 1:
                    coords = coords.reshape(-1, 1)
                if coords.ndim != 2:
                    print(f"[DEBUG] UMAP coords for {u_key} not 2D")
                    return None, None
                if coords.shape[1] < 2:
                    coords = np.hstack([coords, np.zeros((coords.shape[0], 1))])
                coords = coords[:, :2].astype(float)
                clusters = None
                if (l_key is not None) and (l_key in adata.obs):
                    clusters = adata.obs[l_key].astype(str)
                return coords, clusters

            default_leiden = keys.get('leiden_default_key') if keys else None
            user_leiden = keys.get('leiden_user_key') if keys else None

            coords_def, clusters_def = get_coords_and_clusters(default_key, default_leiden)
            coords_user, clusters_user = get_coords_and_clusters(user_key, user_leiden)

            fig = make_subplots(rows=1, cols=2, subplot_titles=["Default-cutoff UMAP", "User-cutoff UMAP"])

            def add_panel(col, coords, clusters, title):
                if coords is None:
                    print(f"[DEBUG] No UMAP available for panel {col} ({title})")
                    fig.add_annotation(text="No UMAP available", xref="paper", yref="paper",
                                       x=0.25 if col == 1 else 0.75, y=0.5, showarrow=False)
                    return
                idx = _sample_idx(coords.shape[0], hvg_max_points)
                sampled_coords = coords[idx]
                df = pd.DataFrame({'x': sampled_coords[:, 0], 'y': sampled_coords[:, 1]})
                if clusters is not None:
                    cluster_vals = clusters.values[idx] if hasattr(clusters, 'values') else np.asarray(clusters)[idx]
                    unique = list(pd.Series(cluster_vals).unique())
                    for u in unique:
                        mask = (cluster_vals == u)
                        fig.add_trace(go.Scattergl(x=df['x'][mask], y=df['y'][mask], mode='markers',
                                                   marker=dict(size=point_size_2d), name=str(u), showlegend=(col == 1)), row=1, col=col)
                else:
                    fig.add_trace(go.Scattergl(x=df['x'], y=df['y'], mode='markers', marker=dict(size=point_size_2d),
                                               name=title, showlegend=False), row=1, col=col)

            add_panel(1, coords_def, clusters_def, 'Default')
            add_panel(2, coords_user, clusters_user, 'User')

            fig.update_layout(height=plot_height, showlegend=True, autosize=False)
        except Exception as e:
            print(f"[DEBUG] Exception while building UMAP plots: {e}")
            fig = {'data': [], 'layout': {'height': plot_height, 'autosize': False}}

        try:
            if processing_plot_type in ['leiden_n']:
                coords, clusters = get_coords_and_clusters(chosen_key, user_leiden if umap_select == 'user' else default_leiden)
                if coords is None:
                    print(f"[DEBUG] No coords for chosen UMAP {chosen_key}")
                    return hvg_list, f"Selected {len(hvgs)} HVGs (user mask)", fig
                traces = []
                obs = adata.obs
                idx = _sample_idx(coords.shape[0], hvg_max_points)
                sampled_coords = coords[idx]
                sampled_obs = obs.iloc[idx]
                for i, val in enumerate(sorted(sampled_obs['leiden_n'].unique())):
                    mask = sampled_obs['leiden_n'] == val
                    b = sampled_coords[mask.values]
                    traces.append(go.Scattergl(x=b[:, 0], y=b[:, 1], text=("Cell ID: " + sampled_obs.loc[mask, 'cell_ID']), mode='markers', marker={'size': point_size_2d}, name=("Cluster " + str(val))))
                single_fig = {'data': traces, 'layout': dict(xaxis={'title': 'UMAP 1'}, yaxis={'title': 'UMAP 2'}, margin=margin, hovermode='closest', height=plot_height, autosize=False)}
                return hvg_list, f"Selected {len(hvgs)} HVGs (user mask)", single_fig
            elif processing_plot_type in ['total_counts', 'log1p_total_counts', 'n_genes']:
                coords, clusters = get_coords_and_clusters(chosen_key, None)
                if coords is None:
                    print(f"[DEBUG] No coords for chosen UMAP {chosen_key}")
                    return hvg_list, f"Selected {len(hvgs)} HVGs (user mask)", fig
                selected_gene = processing_plot_type
                values = None
                if selected_gene in adata.obs.columns:
                    values = adata.obs[selected_gene].values
                elif hasattr(adata, 'obs_vector'):
                    try:
                        values = adata.obs_vector(selected_gene)
                    except Exception:
                        values = None
                idx = _sample_idx(coords.shape[0], hvg_max_points)
                df = pd.DataFrame({'x': coords[idx, 0], 'y': coords[idx, 1], 'val': values[idx] if values is not None else None})
                trace = go.Scattergl(x=df['x'], y=df['y'], mode='markers', marker={'size': point_size_2d, 'color': df['val'], 'colorscale': 'viridis', 'colorbar': dict(title=str(selected_gene))})
                single_fig = {'data': [trace], 'layout': dict(xaxis={'title': 'UMAP 1'}, yaxis={'title': 'UMAP 2'}, margin=margin, hovermode='closest', height=plot_height, autosize=False)}
                return hvg_list, f"Selected {len(hvgs)} HVGs (user mask)", single_fig
        except Exception as e:
            print(f"[DEBUG] Exception while building single-panel plot: {e}")
            pass

        status = f"Selected {len(hvgs)} HVGs (user mask)"
        print(f"[DEBUG] Returning status: {status}")
        return hvg_list, status, fig

    except Exception as e:
        print(f"[DEBUG] Exception in recalc_hvgs: {e}")
        import traceback; traceback.print_exc()
        return "", f"HVG calculation failed: {e}", {}
