"""
Dash layouts for correlation network analysis interface.
"""

from dash import dcc, html
import dash_bootstrap_components as dbc

from correlations.correlation_components import (
    create_correlation_controls,
    create_correlation_upload_controls,
    correlation_network_components
)


def correlation_tab_layout():
    """Create the main layout for correlation network analysis tab."""
    return html.Div([
        correlation_network_components,
        
        dbc.Container([
            dbc.Row([
                dbc.Col([
                    html.H2("Correlation Network Analysis", className="mb-3"),
                    html.P(
                        "Analyze gene-gene correlation networks using modified Pearson correlation "
                        "coefficients (mcPCCs). Communities are detected using the Leiden algorithm, "
                        "and inter-community relationships are visualized.",
                        className="text-muted mb-4"
                    )
                ])
            ]),
            
            dbc.Row([
                # Left sidebar - controls
                dbc.Col([
                    create_correlation_upload_controls(),
                    create_correlation_controls(),
                ], width=3),
                
                # Main visualization area
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader("Community Centroid Graph"),
                        dbc.CardBody([
                            dcc.Loading(
                                id="loading_centroid_plot",
                                type="default",
                                children=dcc.Graph(
                                    id='centroid_graph',
                                    config={
                                        'displayModeBar': True,
                                        'displaylogo': False,
                                        'modeBarButtonsToRemove': ['lasso2d', 'select2d']
                                    },
                                    style={'height': '500px'}
                                )
                            )
                        ])
                    ], className="mb-3"),
                    
                    dbc.Card([
                        dbc.CardHeader("Gene Network"),
                        dbc.CardBody([
                            dcc.Loading(
                                id="loading_network_plot",
                                type="default",
                                children=dcc.Graph(
                                    id='gene_network_graph',
                                    config={
                                        'displayModeBar': True,
                                        'displaylogo': False,
                                        'modeBarButtonsToRemove': ['lasso2d', 'select2d']
                                    },
                                    style={'height': '600px'}
                                )
                            )
                        ])
                    ], className="mb-3"),
                    
                    dbc.Card([
                        dbc.CardHeader("Community Statistics"),
                        dbc.CardBody([
                            html.Div(id='correlation_statistics_table')
                        ])
                    ])
                ], width=9)
            ])
        ], fluid=True)
    ])


def create_correlation_statistics_table(var_with_rank, communities_dict):
    """Create a summary table of community statistics."""
    if var_with_rank is None or communities_dict is None:
        return html.P("No data loaded", className="text-muted")
    
    # Calculate statistics per community
    stats_rows = []
    for community_id in sorted(communities_dict.keys()):
        community_genes = communities_dict[community_id]
        community_df = var_with_rank[var_with_rank['Community'] == community_id]
        
        avg_centrality = community_df['Centrality'].mean() if 'Centrality' in community_df.columns else 0
        max_centrality = community_df['Centrality'].max() if 'Centrality' in community_df.columns else 0
        top_gene = community_df.nlargest(1, 'Centrality').index[0] if len(community_df) > 0 else 'N/A'
        
        stats_rows.append(
            html.Tr([
                html.Td(str(community_id)),
                html.Td(str(len(community_genes))),
                html.Td(f"{avg_centrality:.4f}"),
                html.Td(f"{max_centrality:.4f}"),
                html.Td(top_gene)
            ])
        )
    
    table = dbc.Table([
        html.Thead(
            html.Tr([
                html.Th("Community ID"),
                html.Th("# Genes"),
                html.Th("Avg Centrality"),
                html.Th("Max Centrality"),
                html.Th("Top Gene")
            ])
        ),
        html.Tbody(stats_rows)
    ], striped=True, bordered=True, hover=True, responsive=True, size="sm")
    
    return table
