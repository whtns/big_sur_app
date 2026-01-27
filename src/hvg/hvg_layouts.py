from dash import html, dcc
import dash_bootstrap_components as dbc


def hvg_layout():
    return dbc.Container(fluid=True, children=[
        dcc.Store(id='hvg-initial-load-trigger'),
        dbc.Row(children=[
            dbc.Col(children=[
                html.H3("Highly Variable Genes (HVG)"),
                html.P("Use this panel to inspect and (re)compute HVGs for the current dataset."),
                html.Div([
                    html.Label("mcFano cutoff (mcfano_cutoff):"),
                    dcc.Input(id='mcfano_cutoff', type='number', value=1.0, min=0.0, step=0.1)
                ], style={'marginBottom': 6}),
                html.Div([
                    html.Label("Default mcFano cutoff (default_mcfano_cutoff):"),
                    dcc.Input(id='default_mcfano_cutoff', type='number', value=0.9, min=0.0, step=0.1)
                ], style={'marginBottom': 6}),
                html.Div([
                    html.Label("Default p-value cutoff (default_pvalue_cutoff):"),
                    dcc.Input(id='default_pvalue_cutoff', type='number', value=0.05, min=0.0, step=1e-4)
                ], style={'marginBottom': 6}),
                html.Div([
                    html.Label("UMAP to display:"),
                    dcc.RadioItems(id='hvg_umap_select', options=[
                        {'label': 'User cutoff', 'value': 'user'},
                        {'label': 'Default cutoff', 'value': 'default'}
                    ], value='user')
                ], style={'marginBottom': 6}),
                html.Div([
                    html.Label("UMAP plot type:"),
                    dcc.Dropdown(
                        id='hvg_UMAP_dropdown',
                        options=[
                            {'label': 'leiden', 'value': 'leiden_n'},
                            {'label': '# UMIs', 'value': 'total_counts'},
                            {'label': "# UMIs [ln(1+UMIs)]", 'value': 'log1p_total_counts'},
                            {'label': '# unique genes', 'value': 'n_genes'},
                        ],
                        value='leiden_n',
                        multi=False,
                        searchable=True
                    )
                ], style={'marginBottom': 6}),
                html.Div([
                    html.Label("Dimensions:"),
                    dcc.RadioItems(id='hvg_n_dims_radio', options=[{'label': '2D', 'value': 2}, {'label': '3D', 'value': 3}], value=2)
                ], style={'marginBottom': 10}),
                html.Div([
                    html.Label("Max points to plot (0 = no downsampling):"),
                    dcc.Input(id='hvg_max_points', type='number', value=10000, min=0, step=1000)
                ], style={'marginBottom': 10}),
                html.Div([
                    html.Label("p-value cutoff (pvalue_cutoff):"),
                    dcc.Input(id='pvalue_cutoff', type='number', value=0.01, min=0.0, step=1e-4)
                ], style={'marginBottom': 10}),
                dbc.Button("Recalculate HVGs", id='recalc_hvgs_button', color='primary'),
                html.Div(id='hvg_status', style={'marginTop': 10}),
            ], width=6),
            dbc.Col(children=[
                html.H4("Top HVGs"),
                dcc.Graph(id='hvg_umap_graph', config={'responsive': False}, style={'height': '600px', 'width': '100%'}),
                html.Div(id='hvg_list', style={'whiteSpace': 'pre-wrap', 'maxHeight': '300px', 'overflowY': 'auto'})
            ], width=6)
        ])
    ])
    return m
