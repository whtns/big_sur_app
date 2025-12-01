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
    [Output('hvg_list', 'children'), Output('hvg_status', 'children'), Output('hvg_umap_graph', 'figure')],
    [Input('recalc_hvgs_button', 'n_clicks'), Input('hvg_UMAP_dropdown', 'value'), Input('hvg_n_dims_radio', 'value'), Input('hvg_umap_select', 'value')],
    [State('session-id', 'children'), State('mcfano_cutoff', 'value'), State('pvalue_cutoff', 'value'), State('hvg_max_points', 'value'), State('default_mcfano_cutoff', 'value'), State('default_pvalue_cutoff', 'value')]
)
def recalc_hvgs(n_clicks, processing_plot_type, n_dim_proj_plot, umap_select, session_ID, mcfano_cutoff, pvalue_cutoff, hvg_max_points, default_mcfano_cutoff, default_pvalue_cutoff):
    # This callback is triggered both by the Recalculate HVGs button and by
    # UMAP selection controls. If nothing has been clicked and controls are
    # uninitialized, bail out.
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update, dash.no_update

    adata = cache_adata(session_ID)
    if adata is None:
        return "", "No dataset cached for this session", {}

    try:
        mcfano_cutoff = float(mcfano_cutoff) if mcfano_cutoff is not None else 1.0
        pvalue_cutoff = float(pvalue_cutoff) if pvalue_cutoff is not None else 0.05

        # Ensure raw counts layer exists for mcFano (some adata may not have been preprocessed)
        try:
            if 'counts' not in getattr(adata, 'layers', {}):
                try:
                    adata.layers['counts'] = adata.X.copy()
                except Exception:
                    # last-ditch: convert to dense then copy
                    try:
                        adata.layers['counts'] = adata.X.A.copy() if hasattr(adata.X, 'A') else np.asarray(adata.X).copy()
                    except Exception:
                        pass
        except Exception:
            # ignore layering failures; the called function will raise a clearer error
            pass

        # Run mcFano selection and user/default masks + reductions
        try:
            # coerce default cutoff
            try:
                default_mcFano = float(default_mcfano_cutoff) if default_mcfano_cutoff is not None else 0.9
            except Exception:
                default_mcFano = 0.9
            try:
                default_pval = float(default_pvalue_cutoff) if default_pvalue_cutoff is not None else 0.05
            except Exception:
                default_pval = 0.05
            adata, keys = apply_feature_selection_and_reductions(
                adata,
                mcfano_cutoff=mcfano_cutoff,
                pvalue_cutoff=pvalue_cutoff,
                default_mcFano_threshold=default_mcFano,
                default_pvalue_threshold=default_pval,
                cacache=None,
            )
        except Exception as e:
            tb = traceback.format_exc()
            # include available layer/obs/var keys to help debugging
            layers = list(getattr(adata, 'layers', {}).keys()) if adata is not None else []
            obs_cols = list(getattr(adata, 'obs', []).columns) if (adata is not None and hasattr(adata, 'obs')) else []
            var_cols = list(getattr(adata, 'var', []).columns) if (adata is not None and hasattr(adata, 'var')) else []
            err_msg = f"HVG calculation failed: {e} | layers={layers} | obs_cols={obs_cols[:10]} | var_cols={var_cols[:10]}"
            print(tb)
            return "", err_msg, {}

        cache_history(session_ID, history=(f"Selected HVGs via mcFano cutoff={mcfano_cutoff}, p={pvalue_cutoff}"))
        cache_adata(session_ID, adata)

        # get user-selected HVGs
        hvgs = []
        if 'highly_variable_user' in adata.var.columns:
            hvgs = list(adata.var[adata.var['highly_variable_user']].index)

        if len(hvgs) == 0:
            hvg_list = html.Div("No HVGs found with these cutoffs")
        else:
            items = [html.Li(g) for g in hvgs]
            hvg_list = html.Ul(items)

        # Decide which UMAP keys to use
        default_key = keys.get('umap_default_key') if keys else None
        user_key = keys.get('umap_user_key') if keys else None
        # pick requested UMAP (preference: chosen, then fallback)
        chosen_key = user_key if umap_select == 'user' else default_key

        # Build a two-panel Plotly subplot (Default | User) UMAP for comparison
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
                    return None, None
                coords = np.asarray(adata.obsm[u_key])
                if coords.ndim == 1:
                    coords = coords.reshape(-1, 1)
                if coords.ndim != 2:
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
                    # annotate empty panel
                    fig.add_annotation(text="No UMAP available", xref="paper", yref="paper",
                                       x=0.25 if col == 1 else 0.75, y=0.5, showarrow=False)
                    return
                # apply downsampling for performance
                idx = _sample_idx(coords.shape[0], hvg_max_points)
                sampled_coords = coords[idx]
                df = pd.DataFrame({'x': sampled_coords[:, 0], 'y': sampled_coords[:, 1]})
                if clusters is not None:
                    # clusters may be a pandas Series; get sampled cluster values
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
        except Exception:
            fig = {'data': [], 'layout': {'height': plot_height, 'autosize': False}}

        # If the user requested a specific display type (clusters or expression),
        # allow the user selection in `processing_plot_type` to override the two-panel view.
        try:
            if processing_plot_type in ['leiden_n']:
                # return the cluster plot for the chosen UMAP selection
                coords, clusters = get_coords_and_clusters(chosen_key, user_leiden if umap_select == 'user' else default_leiden)
                if coords is None:
                    return hvg_list, f"Selected {len(hvgs)} HVGs (user mask)", fig
                # build a single-panel figure similar to plot_UMAP with downsampling
                traces = []
                obs = adata.obs
                # sample indices first
                idx = _sample_idx(coords.shape[0], hvg_max_points)
                sampled_coords = coords[idx]
                sampled_obs = obs.iloc[idx]
                # iterate cluster labels present in sampled_obs
                for i, val in enumerate(sorted(sampled_obs['leiden_n'].unique())):
                    mask = sampled_obs['leiden_n'] == val
                    b = sampled_coords[mask.values]
                    traces.append(go.Scattergl(x=b[:, 0], y=b[:, 1], text=("Cell ID: " + sampled_obs.loc[mask, 'cell_ID']), mode='markers', marker={'size': point_size_2d}, name=("Cluster " + str(val))))
                single_fig = {'data': traces, 'layout': dict(xaxis={'title': 'UMAP 1'}, yaxis={'title': 'UMAP 2'}, margin=margin, hovermode='closest', height=plot_height, autosize=False)}
                return hvg_list, f"Selected {len(hvgs)} HVGs (user mask)", single_fig
            elif processing_plot_type in ['total_counts', 'log1p_total_counts', 'n_genes']:
                # expression plot on chosen UMAP
                coords, clusters = get_coords_and_clusters(chosen_key, None)
                if coords is None:
                    return hvg_list, f"Selected {len(hvgs)} HVGs (user mask)", fig
                selected_gene = processing_plot_type
                # get expression vector; prefer obs column, otherwise use obs_vector
                values = None
                if selected_gene in adata.obs.columns:
                    values = adata.obs[selected_gene].values
                elif hasattr(adata, 'obs_vector'):
                    try:
                        values = adata.obs_vector(selected_gene)
                    except Exception:
                        values = None
                # sample indices
                idx = _sample_idx(coords.shape[0], hvg_max_points)
                df = pd.DataFrame({'x': coords[idx, 0], 'y': coords[idx, 1], 'val': values[idx] if values is not None else None})
                trace = go.Scattergl(x=df['x'], y=df['y'], mode='markers', marker={'size': point_size_2d, 'color': df['val'], 'colorscale': 'viridis', 'colorbar': dict(title=str(selected_gene))})
                single_fig = {'data': [trace], 'layout': dict(xaxis={'title': 'UMAP 1'}, yaxis={'title': 'UMAP 2'}, margin=margin, hovermode='closest', height=plot_height, autosize=False)}
                return hvg_list, f"Selected {len(hvgs)} HVGs (user mask)", single_fig
        except Exception:
            pass

        status = f"Selected {len(hvgs)} HVGs (user mask)"
        return hvg_list, status, fig

    except Exception as e:
        return "", f"HVG calculation failed: {e}", {}
