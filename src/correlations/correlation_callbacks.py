"""
Dash callbacks for correlation network analysis interface.
"""

import dash
from dash import Input, Output, State, callback_context
from dash import dcc, html
import dash_bootstrap_components as dbc
import pickle
import base64
import io
import os
from scipy.sparse import load_npz, save_npz, csr_matrix
import json
import numpy as np
import pandas as pd
from celery.result import AsyncResult

from helper_functions import *
from correlations.correlation_functions import *
from correlations.correlation_plotting import *
from correlations.correlation_layouts import create_correlation_statistics_table
from tasks.tasks import compute_correlations_task
from tasks.celery import task_queue
import plotly.graph_objs as go

# Prefer local correlations implementation; fall back to external BigSuR if needed.
BIGSUR_AVAILABLE = False
calculate_correlations = None
try:
    from correlations.correlation_functions import calculate_correlations
    BIGSUR_AVAILABLE = True
except Exception:
    try:
        from BigSur.correlations import calculate_correlations
        BIGSUR_AVAILABLE = True
    except Exception:
        BIGSUR_AVAILABLE = False
        import logging
        logging.getLogger(__name__).warning(
            "BigSuR.correlations not available; correlation computation from session will be disabled"
        )


def register_correlation_callbacks(app):
    """Register all callbacks for correlation network analysis."""
    
    @app.callback(
        [Output('correlation_data_store', 'data'),
         Output('correlation_graph_store', 'data'),
         Output('correlation_var_rank_store', 'data'),
         Output('correlation_communities_store', 'data'),
         Output('community_0_dropdown', 'options'),
         Output('community_1_dropdown', 'options'),
         Output('correlation_load_status', 'children', allow_duplicate=True)],
        [Input('main_tabs', 'active_tab')],
        [State('session-id', 'children')],
        prevent_initial_call=True
    )
    def initial_correlation_load(active_tab, session_ID):
        """
        Automatically load pre-computed correlation data when the Correlations tab
        is selected.
        """
        if active_tab != 'correlation_tab':
            raise dash.exceptions.PreventUpdate

        correlation_dir = "/app/data/correlation_results"
        mcPCCs_path = os.path.join(correlation_dir, 'mcPCCs.npz')
        pvalues_path = os.path.join(correlation_dir, 'BH_corrected_pvalues.npz')

        if not (os.path.exists(mcPCCs_path) and os.path.exists(pvalues_path)):
            return [dash.no_update] * 6 + [html.P("Pre-computed correlation data not found.", className="text-muted")]

        adata = cache_adata(session_ID)
        if adata is None:
            return [dash.no_update] * 6 + [html.P("Load a dataset before viewing correlations.", className="text-warning")]
        
        try:
            # Fast path: use the pre-warmed pickle built at startup.
            correlation_cache_path = os.path.join(correlation_dir, "correlation_cache.pkl")
            if os.path.exists(correlation_cache_path):
                print(f"[{session_ID}] Using pre-warmed correlation cache.")
                with open(correlation_cache_path, "rb") as f:
                    G, var_with_rank, communities_dict = pickle.load(f)
            else:
                print(f"[{session_ID}] Loading pre-computed correlation data (no cache)...")
                G, var_with_rank, communities_dict = load_correlation_data(
                    session_ID,
                    mcPCCs_path=mcPCCs_path,
                    pvalues_path=pvalues_path,
                )
            
            session_dir = os.path.join(save_analysis_path, str(session_ID))
            graph_pickle_path = os.path.join(session_dir, "correlation_graph.pkl")
            with open(graph_pickle_path, 'wb') as f:
                pickle.dump(G, f)

            var_rank_dict = var_with_rank.to_dict('index')
            community_options = [
                {'label': f'Community {i} ({len(communities_dict[i])} genes)', 'value': i}
                for i in sorted(communities_dict.keys())
            ]
            communities_dict_serializable = {
                str(k): v.tolist() for k, v in communities_dict.items()
            }
            
            status_msg = html.Div([
                html.I(className="fas fa-check-circle text-success me-2"),
                f"Loaded pre-computed correlations: {len(G.nodes())} genes, {len(communities_dict)} communities"
            ], className="text-success")

            return (
                {'loaded': True, 'source': 'pre-computed'},
                {'graph_path': graph_pickle_path},
                var_rank_dict,
                communities_dict_serializable,
                community_options,
                community_options,
                status_msg
            )
        except Exception as e:
            print(f"[{session_ID}] Error loading pre-computed correlations: {e}")
            return [dash.no_update] * 6 + [html.Div(f"Error: {e}", className="text-danger")]


    # Combined callback for both triggering and polling correlation computation
    @app.callback(
        [Output('correlation_task_id_store', 'data'),
         Output('correlation_load_status', 'children'),
         Output('correlation_progress_interval', 'disabled'),
         Output('correlation_data_store', 'data', allow_duplicate=True),
         Output('correlation_graph_store', 'data', allow_duplicate=True),
         Output('correlation_var_rank_store', 'data', allow_duplicate=True),
         Output('correlation_communities_store', 'data', allow_duplicate=True),
         Output('community_0_dropdown', 'options', allow_duplicate=True),
         Output('community_1_dropdown', 'options', allow_duplicate=True)],
        [Input('compute_correlations_button', 'n_clicks'),
         Input('correlation_progress_interval', 'n_intervals'),
         Input('correlation_task_id_store', 'data')],
        [State('session-id', 'children')],
        prevent_initial_call=True
    )
    def combined_correlation_callback(n_clicks, n_intervals, task_data, session_ID):
            ctx = dash.callback_context
            if not ctx.triggered:
                return (dash.no_update, "", True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
            trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]

            # If triggered by compute_correlations_button, start computation
            if trigger_id == 'compute_correlations_button':
                if n_clicks is None:
                    return (None, "", True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
                if not BIGSUR_AVAILABLE:
                    error_msg = html.Div([
                        html.I(className="fas fa-exclamation-circle text-danger me-2"),
                        "BigSuR correlation module not available"
                    ], className="text-danger")
                    return (None, error_msg, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
                try:
                    adata = cache_adata(session_ID)
                    if adata is None:
                        error_msg = html.Div("No dataset in current session", className="text-danger")
                        return (None, error_msg, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
                    if 'highly_variable_user' not in adata.var.columns:
                        error_msg = html.Div([
                            html.I(className="fas fa-exclamation-circle text-warning me-2"),
                            "No HVG analysis found. Please run HVG analysis first."
                        ], className="text-warning")
                        return (None, error_msg, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
                    hvg_count = adata.var['highly_variable_user'].sum()
                    if hvg_count == 0:
                        error_msg = html.Div("No highly variable genes selected", className="text-warning")
                        return (None, error_msg, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
                    task = compute_correlations_task.delay(session_ID)
                    status_msg = html.Div([
                        html.I(className="fas fa-spinner fa-spin text-info me-2"),
                        f"Computing correlations for {hvg_count} HVGs..."
                    ], className="text-info")
                    return ({'task_id': task.id, 'session_ID': session_ID}, status_msg, False, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
                except Exception as e:
                    error_msg = html.Div([
                        html.I(className="fas fa-exclamation-circle text-danger me-2"),
                        f"Error: {str(e)}"
                    ], className="text-danger")
                    return (None, error_msg, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])

            # Otherwise, poll task progress
            if task_data is None or 'task_id' not in task_data:
                return (dash.no_update, "", True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
            task_id = task_data['task_id']
            session_ID = task_data.get('session_ID')
            try:
                task_result = task_queue.AsyncResult(task_id)
                state = task_result.state
                if state == 'PENDING':
                    return (dash.no_update, "", False, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
                elif state == 'PROGRESS':
                    try:
                        status_msg = task_result.info.get('status', 'Computing...') if task_result.info else 'Computing...'
                    except Exception:
                        status_msg = 'Computing...'
                    progress_display = html.Div([
                        html.I(className="fas fa-spinner fa-spin text-info me-2"),
                        status_msg
                    ], className="text-info")
                    return (dash.no_update, progress_display, False, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
                elif state == 'SUCCESS':
                    result = task_result.result
                    if result is None:
                        error_msg = html.Div([
                            html.I(className="fas fa-exclamation-circle text-danger me-2"),
                            "Task completed but returned no data"
                        ], className="text-danger")
                        return (dash.no_update, error_msg, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
                    mcPCCs_path = result['mcPCCs_path']
                    pvalues_path = result['pvalues_path']
                    G, var_with_rank, communities_dict = load_correlation_data(
                        session_ID,
                        mcPCCs_path=mcPCCs_path,
                        pvalues_path=pvalues_path,
                        pvalue_threshold=0.001,
                        correlation_threshold=0.1
                    )
                    session_dir = os.path.join(save_analysis_path, str(session_ID))
                    graph_pickle_path = os.path.join(session_dir, "correlation_graph.pkl")
                    with open(graph_pickle_path, 'wb') as f:
                        pickle.dump(G, f)
                    var_rank_dict = var_with_rank.to_dict('index')
                    community_options = [
                        {'label': f'Community {i} ({len(communities_dict[i])} genes)', 'value': i}
                        for i in sorted(communities_dict.keys())
                    ]
                    communities_dict_serializable = {
                        str(k): v.tolist() if isinstance(v, np.ndarray) else list(v)
                        for k, v in communities_dict.items()
                    }
                    success_msg = html.Div([
                        html.I(className="fas fa-check-circle text-success me-2"),
                        f"Computed correlations: {len(G.nodes())} genes, {len(communities_dict)} communities"
                    ], className="text-success")
                    return (
                        dash.no_update, success_msg, True,
                        {'loaded': True, 'source': 'computed'},
                        {'graph_path': graph_pickle_path},
                        var_rank_dict,
                        communities_dict_serializable,
                        community_options,
                        community_options
                    )
                elif state == 'FAILURE':
                    try:
                        error_msg_text = str(task_result.result) if task_result.result else 'Unknown error'
                    except Exception:
                        error_msg_text = 'Task failed with unknown error'
                    error_msg = html.Div([
                        html.I(className="fas fa-exclamation-circle text-danger me-2"),
                        f"Computation failed: {error_msg_text}"
                    ], className="text-danger")
                    return (dash.no_update, error_msg, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
                return (dash.no_update, "", False, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])
            except Exception as e:
                error_msg = html.Div([
                    html.I(className="fas fa-exclamation-circle text-danger me-2"),
                    f"Error polling task: {str(e)}"
                ], className="text-danger")
                return (dash.no_update, error_msg, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], [])

    # Update plots and statistics once data is available or controls change
    @app.callback(
        [Output('centroid_graph', 'figure'),
         Output('gene_network_graph', 'figure'),
         Output('correlation_statistics_table', 'children')],
        [Input('update_correlation_network_button', 'n_clicks'),
         Input('community_0_dropdown', 'value'),
         Input('community_1_dropdown', 'value'),
         Input('k_centroids_slider', 'value'),
         Input('k_network_slider', 'value'),
         Input('correlation_n_top_genes_slider', 'value'),
         Input('positive_edge_cutoff_slider', 'value'),
         Input('negative_edge_cutoff_slider', 'value'),
         Input('calculate_spring_checkbox', 'value'),
         Input('recalculate_position_checkbox', 'value'),
         Input('correlation_graph_store', 'data'),
         Input('correlation_var_rank_store', 'data'),
         Input('correlation_communities_store', 'data')],
        prevent_initial_call=True
    )
    def update_correlation_visuals(
        _n_clicks,
        community_0,
        community_1,
        k_centroids,
        k_network,
        n_top_genes,
        positive_edge_cutoff,
        negative_edge_cutoff,
        calculate_spring,
        recalc_position,
        graph_store,
        var_rank_store,
        communities_store
    ):
        # Ensure we have required data
        if not graph_store or 'graph_path' not in graph_store:
            empty_fig = go.Figure()
            return empty_fig, empty_fig, html.P("No correlation data loaded", className="text-muted")

        graph_path = graph_store['graph_path']
        if not os.path.exists(graph_path) or var_rank_store is None or communities_store is None:
            empty_fig = go.Figure()
            return empty_fig, empty_fig, html.P("No correlation data loaded", className="text-muted")

        # Load graph
        with open(graph_path, 'rb') as f:
            G = pickle.load(f)

        # Restore var_with_rank dataframe
        var_with_rank = pd.DataFrame.from_dict(var_rank_store, orient='index')
        if 'Community' in var_with_rank.columns:
            var_with_rank['Community'] = var_with_rank['Community'].astype(int)

        # Restore communities dict with integer keys
        communities_dict = {int(k): np.array(v) for k, v in communities_store.items()}

        # Rebuild communities with top-gene cap if requested
        communities_dict = make_communities_dict(var_with_rank, n_top_genes=n_top_genes)

        # Ensure selected communities exist; default to first available
        available_keys = sorted(communities_dict.keys())
        if not available_keys:
            empty_fig = go.Figure()
            return empty_fig, empty_fig, html.P("No communities detected", className="text-muted")

        community_0 = community_0 if community_0 in available_keys else available_keys[0]
        community_1 = community_1 if community_1 in available_keys else available_keys[0]

        centroid_fig, network_fig = plot_correlation_network(
            G,
            var_with_rank,
            communities_dict,
            k_centroids=k_centroids,
            k_network=k_network,
            seed=0,
            n_top_genes=n_top_genes,
            positive_edge_weight_cutoff=positive_edge_cutoff,
            negative_edge_weight_cutoff=negative_edge_cutoff,
            community_0=community_0,
            community_1=community_1,
            calculate_spring=bool(calculate_spring),
            recalculate_position=bool(recalc_position)
        )

        stats_table = create_correlation_statistics_table(var_with_rank, communities_dict)

        return centroid_fig, network_fig, stats_table

