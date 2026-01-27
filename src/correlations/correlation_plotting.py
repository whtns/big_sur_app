"""
Plotly-based visualization functions for correlation networks.
Converts matplotlib/networkx visualizations to interactive Plotly graphs.
"""

import plotly.graph_objs as go
import numpy as np
import pandas as pd
import seaborn as sns

from correlations.correlation_functions import *


def create_centroid_plot_data(G_centroid, centroid_locations, communities_dict):
    """
    Create Plotly traces for centroid graph visualization.
    
    Returns:
    --------
    traces : list
        List of Plotly graph objects for nodes and edges
    color_mapping : dict
        Dictionary mapping community ID to color
    """
    traces = []
    
    # Assign colors to communities
    n_communities = len(communities_dict)
    colors = sns.color_palette("hls", n_communities)
    color_mapping = {
        community: f'rgb({int(c[0]*255)},{int(c[1]*255)},{int(c[2]*255)})'
        for community, c in zip(communities_dict.keys(), colors)
    }
    
    # Calculate edge weights for scaling
    weights = [G_centroid[u][v]['pos_weight'] for u, v in G_centroid.edges()]
    weights.extend([G_centroid[u][v]['neg_weight'] for u, v in G_centroid.edges()])
    max_weight = max(weights) if weights else 1
    min_weight = min(weights) if weights else 0
    
    # Draw positive edges (red)
    pos_edge_traces = create_edge_traces(
        G_centroid, 
        centroid_locations, 
        'pos_weight', 
        'red',
        max_weight, 
        min_weight
    )
    traces.extend(pos_edge_traces)
    
    # Draw negative edges (blue)
    neg_edge_traces = create_edge_traces(
        G_centroid,
        centroid_locations,
        'neg_weight',
        'blue',
        max_weight,
        min_weight
    )
    traces.extend(neg_edge_traces)
    
    # Draw nodes with community labels
    node_x = []
    node_y = []
    node_text = []
    node_colors = []
    
    for node in G_centroid.nodes():
        x, y = centroid_locations[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(f'Community {node}<br>{len(communities_dict[node])} genes')
        node_colors.append(color_mapping[node])
    
    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode='markers+text',
        text=[str(node) for node in G_centroid.nodes()],
        textposition='middle center',
        textfont=dict(size=10, color='white'),
        hovertext=node_text,
        hoverinfo='text',
        marker=dict(
            size=30,
            color=node_colors,
            line=dict(width=2, color='black')
        ),
        showlegend=False
    )
    traces.append(node_trace)
    
    return traces, color_mapping


def create_edge_traces(G, positions, weight_key, color, max_weight, min_weight):
    """Create edge traces for Plotly with varying widths."""
    traces = []
    
    for u, v in G.edges():
        weight = G[u][v][weight_key]
        if weight > 0:
            x0, y0 = positions[u]
            x1, y1 = positions[v]
            
            # Scale width
            max_width = 5
            min_width = 0.5
            if max_weight > min_weight:
                width = min_width + (weight - min_weight) / (max_weight - min_weight) * (max_width - min_width)
            else:
                width = min_width
            
            edge_trace = go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode='lines',
                line=dict(width=width, color=color),
                hoverinfo='none',
                showlegend=False
            )
            traces.append(edge_trace)
    
    return traces


def create_network_plot_data(
    G_subset_edges,
    position,
    positive_edges_to_plot,
    negative_edges_to_plot,
    communities_dict,
    color_mapping,
    community_0,
    community_1
):
    """
    Create Plotly traces for gene network visualization.
    
    Returns:
    --------
    traces : list
        List of Plotly graph objects
    """
    traces = []
    
    # Draw positive edges (red)
    for u, v in positive_edges_to_plot:
        if u in position and v in position:
            x0, y0 = position[u]
            x1, y1 = position[v]
            edge_trace = go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode='lines',
                line=dict(width=1, color='red'),
                hoverinfo='none',
                showlegend=False
            )
            traces.append(edge_trace)
    
    # Draw negative edges (blue)
    for u, v in negative_edges_to_plot:
        if u in position and v in position:
            x0, y0 = position[u]
            x1, y1 = position[v]
            edge_trace = go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode='lines',
                line=dict(width=1, color='blue'),
                hoverinfo='none',
                showlegend=False
            )
            traces.append(edge_trace)
    
    # Draw nodes by community
    for community in [community_0, community_1]:
        genes_of_community = communities_dict[community]
        genes_in_network = [g for g in genes_of_community if g in G_subset_edges.nodes()]
        
        if genes_in_network:
            node_x = [position[gene][0] for gene in genes_in_network]
            node_y = [position[gene][1] for gene in genes_in_network]
            
            node_trace = go.Scatter(
                x=node_x,
                y=node_y,
                mode='markers+text',
                text=genes_in_network,
                textposition='top center',
                textfont=dict(size=8),
                hovertext=genes_in_network,
                hoverinfo='text',
                marker=dict(
                    size=8,
                    color=color_mapping[community],
                    line=dict(width=1, color='black')
                ),
                name=f'Community {community}',
                showlegend=True
            )
            traces.append(node_trace)
    
    return traces


def plot_correlation_network(
    G,
    var_with_rank,
    communities_dict,
    k_centroids=0.1,
    k_network=0.1,
    seed=0,
    n_top_genes=100,
    positive_edge_weight_cutoff=0.1,
    negative_edge_weight_cutoff=-0.1,
    community_0=0,
    community_1=0,
    calculate_spring=True,
    recalculate_position=True
):
    """
    Create both centroid and gene network plots as Plotly figures.
    
    Returns:
    --------
    centroid_fig : plotly.graph_objs.Figure
        Centroid graph figure
    network_fig : plotly.graph_objs.Figure
        Gene network figure
    """
    # Build centroid graph
    G_centroid, df = build_centroid_graph(G, communities_dict)
    centroid_locations = calculate_centroid_layout(
        G_centroid, 
        k_centroids=k_centroids, 
        seed=seed, 
        calculate_spring=calculate_spring
    )
    
    # Create centroid plot
    centroid_traces, color_mapping = create_centroid_plot_data(
        G_centroid, 
        centroid_locations, 
        communities_dict
    )
    
    centroid_fig = go.Figure(
        data=centroid_traces,
        layout=go.Layout(
            title=dict(text='Community Centroids', font=dict(size=16)),
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white'
        )
    )
    
    # Create gene network
    G_subset_edges, position, positive_edges_to_plot, negative_edges_to_plot = create_gene_network_subset(
        communities_dict,
        G,
        community_0=community_0,
        community_1=community_1,
        k_network=k_network,
        seed=seed,
        positive_edge_weight_cutoff=positive_edge_weight_cutoff,
        negative_edge_weight_cutoff=negative_edge_weight_cutoff,
        recalculate_position=recalculate_position
    )
    
    network_traces = create_network_plot_data(
        G_subset_edges,
        position,
        positive_edges_to_plot,
        negative_edges_to_plot,
        communities_dict,
        color_mapping,
        community_0,
        community_1
    )
    
    if community_0 == community_1:
        title_text = f'Community {community_0} Network'
    else:
        title_text = f'Community {community_0} and {community_1} Network'
    
    title_text += f'<br><sub>Positive edges > {positive_edge_weight_cutoff:.3f}, Negative edges < {negative_edge_weight_cutoff:.3f}</sub>'
    
    network_fig = go.Figure(
        data=network_traces,
        layout=go.Layout(
            title=dict(text=title_text, font=dict(size=14)),
            showlegend=True,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=60),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white',
            legend=dict(x=0, y=1)
        )
    )
    
    return centroid_fig, network_fig
