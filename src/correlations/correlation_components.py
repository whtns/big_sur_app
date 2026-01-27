"""
Dash components for correlation network analysis interface.
"""

from dash import dcc, html
import dash_bootstrap_components as dbc


def create_correlation_controls():
    """Create control panel for correlation network parameters."""
    return dbc.Card([
        dbc.CardHeader("Network Parameters"),
        dbc.CardBody([
            # Centroid layout parameter
            html.Label("Centroid Spring Constant (k):", className="mt-2"),
            dcc.Slider(
                id='k_centroids_slider',
                min=0.001,
                max=2.0,
                step=0.001,
                value=0.1,
                marks={0.001: '0.001', 0.5: '0.5', 1.0: '1.0', 2.0: '2.0'},
                tooltip={"placement": "bottom", "always_visible": False}
            ),
            
            # Network layout parameter
            html.Label("Network Spring Constant (k):", className="mt-3"),
            dcc.Slider(
                id='k_network_slider',
                min=0.001,
                max=1.0,
                step=0.001,
                value=0.1,
                marks={0.001: '0.001', 0.25: '0.25', 0.5: '0.5', 1.0: '1.0'},
                tooltip={"placement": "bottom", "always_visible": False}
            ),
            
            # Number of top genes
            html.Label("Top Genes per Community:", className="mt-3"),
            dcc.Slider(
                id='correlation_n_top_genes_slider',
                min=10,
                max=200,
                step=5,
                value=100,
                marks={10: '10', 50: '50', 100: '100', 150: '150', 200: '200'},
                tooltip={"placement": "bottom", "always_visible": False}
            ),
            
            # Positive edge cutoff
            html.Label("Positive Edge Cutoff:", className="mt-3"),
            dcc.Slider(
                id='positive_edge_cutoff_slider',
                min=0.0,
                max=1.0,
                step=0.01,
                value=0.1,
                marks={0: '0', 0.25: '0.25', 0.5: '0.5', 0.75: '0.75', 1.0: '1.0'},
                tooltip={"placement": "bottom", "always_visible": False}
            ),
            
            # Negative edge cutoff
            html.Label("Negative Edge Cutoff:", className="mt-3"),
            dcc.Slider(
                id='negative_edge_cutoff_slider',
                min=-1.0,
                max=0.0,
                step=0.01,
                value=-0.1,
                marks={-1.0: '-1.0', -0.75: '-0.75', -0.5: '-0.5', -0.25: '-0.25', 0: '0'},
                tooltip={"placement": "bottom", "always_visible": False}
            ),
            
            # Community selectors
            html.Label("Community 1:", className="mt-3"),
            dcc.Dropdown(
                id='community_0_dropdown',
                options=[],
                value=0,
                clearable=False
            ),
            
            html.Label("Community 2:", className="mt-3"),
            dcc.Dropdown(
                id='community_1_dropdown',
                options=[],
                value=0,
                clearable=False
            ),
            
            # Layout options
            html.Div([
                dbc.Checkbox(
                    id='calculate_spring_checkbox',
                    label="Use Spring Layout for Centroids",
                    value=True,
                    className="mt-3"
                ),
                dbc.Checkbox(
                    id='recalculate_position_checkbox',
                    label="Recalculate Network Layout",
                    value=True,
                    className="mt-2"
                ),
            ]),
            
            # Update button
            dbc.Button(
                "Update Visualization",
                id='update_correlation_network_button',
                color="primary",
                className="mt-3 w-100"
            )
        ])
    ], className="mb-3")


def create_correlation_upload_controls():
    """Create upload controls for correlation data."""
    return dbc.Card([
        dbc.CardHeader("Load Correlation Data"),
        dbc.CardBody([            # Option to use current session data
            dbc.Button(
                "Compute from Current Session",
                id='compute_correlations_button',
                color="info",
                className="mb-3 w-100",
                outline=True
            ),
            html.Div(id='correlation_load_status', className="mb-3"),
            
            html.Hr(),
            html.P("Or upload pre-computed correlation files:", className="text-muted"),
                        html.Label("Upload mcPCCs (.npz):"),
            dcc.Upload(
                id='upload_mcPCCs',
                children=html.Div([
                    'Drag and Drop or ',
                    html.A('Select mcPCCs File')
                ]),
                style={
                    'width': '100%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px 0'
                },
                multiple=False
            ),
            
            html.Label("Upload BH-corrected p-values (.npz):", className="mt-2"),
            dcc.Upload(
                id='upload_pvalues',
                children=html.Div([
                    'Drag and Drop or ',
                    html.A('Select p-values File')
                ]),
                style={
                    'width': '100%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px 0'
                },
                multiple=False
            ),
            
            html.Label("P-value Threshold:", className="mt-3"),
            dcc.Input(
                id='pvalue_threshold_input',
                type='number',
                value=0.001,
                min=0,
                max=1,
                step=0.0001,
                style={'width': '100%'}
            ),
            
            html.Label("Correlation Threshold:", className="mt-3"),
            dcc.Input(
                id='correlation_threshold_input',
                type='number',
                value=0.1,
                min=0,
                max=1,
                step=0.01,
                style={'width': '100%'}
            ),
            
            dbc.Button(
                "Load & Analyze",
                id='load_correlation_data_button',
                color="success",
                className="mt-3 w-100"
            )
        ])
    ], className="mb-3")


correlation_network_components = html.Div([
    dcc.Store(id='correlation_data_store'),
    dcc.Store(id='correlation_graph_store'),
    dcc.Store(id='correlation_var_rank_store'),
    dcc.Store(id='correlation_communities_store'),
    dcc.Store(id='hvg_genes_store'),  # Store HVG genes list
    dcc.Store(id='correlation_task_id_store'),  # Store task ID for polling
    dcc.Interval(id='correlation_progress_interval', interval=500, disabled=True),  # Poll every 500ms
])
