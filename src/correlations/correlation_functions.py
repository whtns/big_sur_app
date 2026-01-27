"""
Correlation network analysis functions for BigSuR.
Adapted from T-cell correlation analysis workflow.
"""

import numpy as np
import pandas as pd
import networkx
import leidenalg
import igraph as ig
from scipy.sparse import csr_matrix, load_npz
from scipy.stats import binomtest
from statsmodels.stats.multitest import fdrcorrection
from itertools import combinations
from networkx.algorithms.centrality import eigenvector_centrality

from helper_functions import *


def create_graph_from_mcPCCs(mcPCCs_thresholded, var):
    """Create graph from mcPCCs matrix."""
    mcPCCs_thresholded_dense = mcPCCs_thresholded.todense()
    
    G = networkx.from_numpy_array(mcPCCs_thresholded_dense)
    nodes_mapper = {i: var.index[i] for i in range(var.shape[0])}
    G = networkx.relabel_nodes(G, nodes_mapper)
    
    return G


def calculate_communities(G, resolution_parameter=1, edge_cutoff=0):
    """
    Run Leiden community detection algorithm on the graph.
    Only considers edges greater than edge_cutoff.
    """
    # Remove negative edges for community detection
    G_positive_only = G.copy()
    remove_edges = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] < edge_cutoff]
    G_positive_only.remove_edges_from(remove_edges)
    
    # Detect communities
    G_ig = ig.Graph.from_networkx(G_positive_only)
    
    # Apply Leiden algorithm
    partition = leidenalg.find_partition(
        G_ig, 
        leidenalg.RBConfigurationVertexPartition, 
        seed=0, 
        resolution_parameter=resolution_parameter
    )
    
    # Get the node-community mapping
    membership = partition.membership
    membership = np.array(membership)
    communities = [
        np.array(list(G_positive_only.nodes))[np.where(membership == i)] 
        for i in range(len(np.unique(membership)))
    ]
    
    return communities, G_positive_only


def calculate_eigenvector_centrality(G_positive_only, communities, var):
    """Calculate eigenvector centrality for each cluster in the graph."""
    rank_df_clusters = pd.DataFrame()
    
    for cluster_number in range(len(communities)):
        cluster = G_positive_only.subgraph(communities[cluster_number])
        if cluster.number_of_nodes() > 1:
            rank = eigenvector_centrality(cluster, max_iter=1000)
            rank_df = pd.DataFrame(data={'Gene': rank.keys(), 'Centrality': rank.values()})
            rank_df['Community'] = cluster_number
            rank_df = rank_df.sort_values('Centrality', ascending=False)
            rank_df['Rank within community'] = range(rank_df.shape[0])
            rank_df_clusters = pd.concat([rank_df_clusters, rank_df])
    
    rank_df_clusters = rank_df_clusters.set_index('Gene')
    var_with_rank = rank_df_clusters.join(var)
    
    return var_with_rank


def detect_significant_number_of_negatively_correlated_edges(G_centroid, communities_dict):
    """
    For each pair of clusters, calculate the number of negatively correlated edges
    and perform a binomial test.
    """
    master_list_positive = np.empty(0)
    master_list_negative = np.empty(0)
    
    # Convert keys to list to ensure type clarity
    community_keys = list(communities_dict.keys())
    generator = combinations(community_keys, 2)
    
    for cluster_number, cluster_number_2 in generator:
        cluster_1 = communities_dict[cluster_number]
        cluster_2 = communities_dict[cluster_number_2]
        number_of_possible_edges = len(cluster_1) * len(cluster_2)
        
        for weight_type in ['pos_weight', 'neg_weight']:
            actual_number_of_correlated_edges = G_centroid.edges[(cluster_number, cluster_number_2)][weight_type]
            actual_number_of_nodes = len(cluster_1) + len(cluster_2)
            
            if actual_number_of_correlated_edges > 0:
                expected_probability = 2 * actual_number_of_correlated_edges / (
                    actual_number_of_nodes * (actual_number_of_nodes - 1)
                )
                test_output = binomtest(
                    actual_number_of_correlated_edges, 
                    number_of_possible_edges, 
                    expected_probability, 
                    alternative='greater'
                )
                temp_dict = {
                    'Cluster 1': cluster_number, 
                    'Cluster 2': cluster_number_2,
                    "Number of nodes in cluster 1": cluster_1.shape[0],
                    "Number of nodes in cluster 2": cluster_2.shape[0],
                    'Number of correlated edges': test_output.k,
                    'Total possible edges': test_output.n,
                    'p_value': test_output.pvalue
                }
                
                if weight_type == 'pos_weight':
                    master_list_positive = np.append(master_list_positive, temp_dict)
                elif weight_type == 'neg_weight':
                    master_list_negative = np.append(master_list_negative, temp_dict)
    
    probability_of_positive_correlations_df = pd.DataFrame.from_records(master_list_positive)
    probability_of_negative_correlations_df = pd.DataFrame.from_records(master_list_negative)
    probability_of_positive_correlations_df['weight type'] = 'positive'
    probability_of_negative_correlations_df['weight type'] = 'negative'
    probability_of_correlations_df = pd.concat(
        (probability_of_positive_correlations_df, probability_of_negative_correlations_df)
    )
    probability_of_correlations_df['FDR corrected p_value'] = fdrcorrection(
        probability_of_correlations_df['p_value']
    )[1]
    probability_of_correlations_df = probability_of_correlations_df.reset_index(drop=True)
    
    return probability_of_correlations_df


def build_centroid_graph(G, communities_dict):
    """
    Build a centroid graph where nodes are communities and edges represent
    inter-community correlations.
    """
    G_centroid = networkx.Graph()
    G_centroid.add_nodes_from(communities_dict.keys())
    generator = combinations(communities_dict.keys(), 2)
    
    # Calculate edges of centroid graph by summing the number of positive and negative
    # cross-community edges
    for community_A, community_B in generator:
        positive_cross_edges = []
        negative_cross_edges = []
        
        for u in communities_dict[community_A]:
            for v in communities_dict[community_B]:
                if G.has_edge(u, v):
                    d = G[u][v]
                    if d['weight'] > 0:
                        positive_cross_edges.append((u, v))
                    elif d['weight'] < 0:
                        negative_cross_edges.append((u, v))
        
        num_positive_cross_edges = len(positive_cross_edges)
        num_negative_cross_edges = len(negative_cross_edges)
        G_centroid.add_edge(
            community_A, 
            community_B,
            pos_weight=num_positive_cross_edges,
            neg_weight=num_negative_cross_edges,
            weight=num_positive_cross_edges - num_negative_cross_edges
        )
    
    df = detect_significant_number_of_negatively_correlated_edges(G_centroid, communities_dict)
    df_of_edges_to_remove = df[(df['FDR corrected p_value'] > 0.05)]
    edges_to_remove_because_of_stats = [
        (df_of_edges_to_remove.loc[row, 'Cluster 1'], df_of_edges_to_remove.loc[row, 'Cluster 2'])
        for row in df_of_edges_to_remove.index
    ]
    G_centroid.remove_edges_from(edges_to_remove_because_of_stats)
    
    return G_centroid, df


def calculate_centroid_layout(G_centroid, k_centroids=0.1, seed=0, calculate_spring=True):
    """Calculate layout for centroid graph."""
    if calculate_spring:
        centroid_locations = networkx.spring_layout(
            G_centroid,
            center=(0, 0),
            k=k_centroids,
            seed=seed,
            weight='weight'
        )
    else:
        # Circular layout
        centroid_locations = points_on_circle(1, len(G_centroid.nodes()))
    
    min_x = min([centroid_locations[node][0] for node in centroid_locations.keys()])
    
    # If any isolated nodes, arrange them separately
    edges_to_remove = [
        (u, v) for u, v, d in G_centroid.edges(data=True) 
        if (d['pos_weight'] == 0) and (d['neg_weight'] == 0)
    ]
    G_centroid.remove_edges_from(edges_to_remove)
    isolated_nodes = list(networkx.isolates(G_centroid))
    new_pos = 1
    for node in isolated_nodes:
        centroid_locations[node] = np.array([min_x - 0.25, new_pos])
        new_pos -= 0.1
    
    return centroid_locations


def points_on_circle(radius, n_points):
    """Generate evenly spaced points on a circle."""
    angles = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
    x = radius * np.cos(angles + np.pi / 4)
    y = radius * np.sin(angles + np.pi / 4)
    pos = {
        i: np.array([x_coord, y_coord]) 
        for i, (x_coord, y_coord) in enumerate(zip(x, y))
    }
    return pos


def make_communities_dict(var_with_rank, n_top_genes=100):
    """Create dictionary of communities with top N genes by centrality."""
    communities_dict = {}
    for community in range(var_with_rank['Community'].nunique()):
        community_df = var_with_rank[var_with_rank['Community'] == community].sort_values(
            by='Centrality', ascending=False
        )
        if community_df.shape[0] < 5:
            break
        else:
            top_genes = community_df[:n_top_genes].index
            temp_dict = {community: top_genes}
            communities_dict.update(temp_dict)
    
    return communities_dict


def create_gene_network_subset(
    communities_dict, 
    G, 
    community_0=0, 
    community_1=1,
    k_network=0.1, 
    seed=0,
    positive_edge_weight_cutoff=0.1,
    negative_edge_weight_cutoff=-0.1,
    recalculate_position=True
):
    """Create a subset network for specified communities with edge filtering."""
    if community_0 == community_1:
        genes_to_plot = communities_dict[community_0]
    else:
        genes_to_plot = np.append(communities_dict[community_0], communities_dict[community_1])
    
    G_subset = G.subgraph(genes_to_plot)
    
    # Subset genes to plot based on edge weights
    positive_edges_to_plot = np.array([
        (u, v) for u, v, d in G_subset.edges(data=True) 
        if d['weight'] > positive_edge_weight_cutoff
    ])
    negative_edges_to_plot = np.array([
        (u, v) for u, v, d in G_subset.edges(data=True) 
        if d['weight'] < negative_edge_weight_cutoff
    ])
    
    exists_positive_edges_to_plot = positive_edges_to_plot.shape[0] > 0
    exists_negative_edges_to_plot = negative_edges_to_plot.shape[0] > 0
    
    if exists_positive_edges_to_plot and exists_negative_edges_to_plot:
        edge_list = np.concatenate((positive_edges_to_plot, negative_edges_to_plot)).tolist()
    elif exists_positive_edges_to_plot:
        edge_list = positive_edges_to_plot.tolist()
    elif exists_negative_edges_to_plot:
        edge_list = negative_edges_to_plot.tolist()
    else:
        edge_list = []
    
    edge_list = [(u, v) for u, v in edge_list]
    
    G_subset_edges = G_subset.edge_subgraph(edge_list).copy()
    isolated_nodes = list(networkx.isolates(G_subset_edges))
    G_subset_edges.remove_nodes_from(isolated_nodes)
    
    if recalculate_position:
        position = networkx.spring_layout(G_subset_edges, k=k_network, seed=seed)
    else:
        position = networkx.spring_layout(G_subset, k=k_network, seed=seed)
    
    return G_subset_edges, position, positive_edges_to_plot, negative_edges_to_plot


def load_correlation_data(session_ID, mcPCCs_path=None, pvalues_path=None, 
                          pvalue_threshold=1e-3, correlation_threshold=0.1):
    """
    Load mcPCCs and p-values, apply thresholds, and create correlation graph.
    
    Parameters:
    -----------
    session_ID : str
        Session identifier
    mcPCCs_path : str, optional
        Path to mcPCCs .npz file
    pvalues_path : str, optional
        Path to BH-corrected p-values .npz file
    pvalue_threshold : float
        P-value threshold for filtering
    correlation_threshold : float
        Absolute correlation threshold for filtering
    
    Returns:
    --------
    G : networkx.Graph
        Correlation graph
    var_with_rank : pd.DataFrame
        Gene rankings with community assignments
    communities_dict : dict
        Dictionary of communities
    """
    # Default paths if not provided
    if mcPCCs_path is None:
        mcPCCs_path = os.path.join(save_analysis_path, str(session_ID), "mcPCCs.npz")
    if pvalues_path is None:
        pvalues_path = os.path.join(save_analysis_path, str(session_ID), "BH_corrected_pvalues.npz")
    
    # Load data
    mcPCCs = load_npz(mcPCCs_path)
    BH_corrected_pvalues = load_npz(pvalues_path)
    
    # Apply thresholds
    mcPCCs_thresholded = mcPCCs.copy()
    mcPCCs_thresholded[BH_corrected_pvalues > pvalue_threshold] = 0
    mcPCCs_thresholded[np.abs(mcPCCs_thresholded) < correlation_threshold] = 0
    
    # Load adata to get var
    adata = cache_adata(session_ID)
    if adata is None:
        raise ValueError(f"Could not load AnnData for session {session_ID}")
    
    # Create graph from correlations
    G = create_graph_from_mcPCCs(mcPCCs_thresholded, adata.var)
    communities, G_positive_only = calculate_communities(G)
    var_with_rank = calculate_eigenvector_centrality(G_positive_only, communities, adata.var)
    var_with_rank['numerical index'] = range(var_with_rank.shape[0])
    
    # Create communities dictionary
    communities_dict = make_communities_dict(var_with_rank, n_top_genes=100)
    
    return G, var_with_rank, communities_dict
