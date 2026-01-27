# Correlation functions
def create_graph_from_mcPCCs(mcPCCs_thresholded, var):
    '''Create graph from mcPCCs matrix.'''
    mcPCCs_thresholded_dense = mcPCCs_thresholded.todense()
    #mcPCCs_thresholded_dense = np.where(np.abs(np.round(mcPCCs_thresholded_dense, 8) == 1), 0, mcPCCs_thresholded_dense)

    G = networkx.from_numpy_array(mcPCCs_thresholded_dense)
    nodes_mapper = {i: var.index[i] for i in range(var.shape[0])}
    G = networkx.relabel_nodes(G, nodes_mapper)

    return G

def calculate_communities(G, resolution_parameter = 1, edge_cutoff = 0):
    '''Run walktrap community detection algorithm on the graph, only considering edges greater than edge_cutoff.'''
    # Remove negative edges for community detection
    G_postive_only = G.copy()
    remove_edges = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] < edge_cutoff]
    G_postive_only.remove_edges_from(remove_edges)

    # Detect communities
    G_ig = ig.Graph.from_networkx(G_postive_only)

    # Apply Leiden algorithm
    partition = leidenalg.find_partition(G_ig, leidenalg.RBConfigurationVertexPartition, seed = 0, resolution_parameter = resolution_parameter)

    # Get the node-community mapping
    membership = partition.membership
    membership = np.array(membership)
    #coms = algorithms.walktrap(G_postive_only)
    #out = coms.communities
    communities = [np.array(list(G_postive_only.nodes))[np.where(membership == i)] for i in range(len(np.unique(membership)))]

    return communities, G_postive_only

def calculate_eigenvector_centrality(G_postive_only, communities, var):
    '''Calculate eigenvector centrality for each cluster in the graph.'''
    # Calculate eigenvector centrality for each cluster
    rank_df_clusters = pd.DataFrame()

    for cluster_number in range(len(communities)):
        cluster = G_postive_only.subgraph(communities[cluster_number])
        if cluster.number_of_nodes() > 1:
            rank = eigenvector_centrality(cluster, max_iter=1000)
            rank_df = pd.DataFrame(data = {'Gene':rank.keys(), 'Centrality':rank.values()})
            rank_df['Community'] = cluster_number
            rank_df = rank_df.sort_values('Centrality', ascending = False)
            rank_df['Rank within community'] = range(rank_df.shape[0])
            rank_df_clusters = pd.concat([rank_df_clusters, rank_df])

    rank_df_clusters = rank_df_clusters.set_index('Gene')
    var_with_rank = rank_df_clusters.join(var)

    return var_with_rank

# Visualization functions
def points_on_circle(radius, n_points):
        angles = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
        x = radius * np.cos(angles + np.pi/4)
        y = radius * np.sin(angles + np.pi/4)
        pos = {i: np.array([x_coord, y_coord]) for i, (x_coord, y_coord) in enumerate(zip(x, y))}
        return pos

def detect_significant_number_of_negatively_correlated_edges(G_centroid, communities_dict):
    '''Optimized version: For each pair of clusters, calculate the number of negatively correlated edges and perform a binomial test.'''
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
                expected_probability =  2*actual_number_of_correlated_edges / (actual_number_of_nodes*(actual_number_of_nodes-1))
                test_output = binomtest(actual_number_of_correlated_edges, number_of_possible_edges, expected_probability, alternative='greater')
                temp_dict = {'Cluster 1': cluster_number, 'Cluster 2': cluster_number_2, "Number of nodes in cluster 1": cluster_1.shape[0], "Number of nodes in cluster 2": cluster_2.shape[0], f'Number of correlated edges':test_output.k, 'Total possible edges':test_output.n, 'p_value': test_output.pvalue}
                if weight_type == 'pos_weight':
                    master_list_positive = np.append(master_list_positive, temp_dict)
                elif weight_type == 'neg_weight':
                        master_list_negative = np.append(master_list_negative, temp_dict)

    probability_of_positive_correlations_df = pd.DataFrame.from_records(master_list_positive)
    probability_of_negative_correlations_df = pd.DataFrame.from_records(master_list_negative)
    probability_of_positive_correlations_df['weight type'] = 'positive'
    probability_of_negative_correlations_df['weight type'] = 'negative'
    probability_of_correlations_df = pd.concat((probability_of_positive_correlations_df, probability_of_negative_correlations_df))
    probability_of_correlations_df['FDR corrected p_value'] = fdrcorrection(probability_of_correlations_df['p_value'])[1]
    probability_of_correlations_df = probability_of_correlations_df.reset_index(drop = True)

    return probability_of_correlations_df

def build_centroid_graph(G, communities_dict):
    G_centroid = networkx.Graph()
    G_centroid.add_nodes_from(communities_dict.keys())
    generator = combinations(communities_dict.keys(), 2)

    ## Calculate edges of centroid graph by summing the number of positive and negative cross-community edges
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
        G_centroid.add_edge(community_A, community_B, pos_weight=num_positive_cross_edges, neg_weight=num_negative_cross_edges, weight = num_positive_cross_edges - num_negative_cross_edges)

    df = detect_significant_number_of_negatively_correlated_edges(G_centroid, communities_dict)
    df_of_edges_to_remove = df[(df['FDR corrected p_value'] > 0.05)]
    edges_to_remove_because_of_stats = [(df_of_edges_to_remove.loc[row, 'Cluster 1'], df_of_edges_to_remove.loc[row, 'Cluster 2']) for row in df_of_edges_to_remove.index]
    G_centroid.remove_edges_from(edges_to_remove_because_of_stats)

    return G_centroid, df


def calculate_centroid_layout(G_centroid, k_centroids=0.1, seed = 0, calculate_spring = True):
## Calculate centroid layout
    if calculate_spring == True:
        centroid_locations = networkx.spring_layout(G_centroid, center = (0, 0), k = k_centroids, method = 'force', weight = 'weight', seed = seed)
    else:
        centroid_locations = points_on_circle(1, len(G_centroid.nodes()))

    min_x = min([centroid_locations[node][0] for node in centroid_locations.keys()])

    ## If any isolated nodes, arrange them separately
    edges_to_remove = [(u, v) for u, v, d in G_centroid.edges(data=True) if (d['pos_weight'] == 0) and (d['neg_weight'] == 0)]
    G_centroid.remove_edges_from(edges_to_remove)
    isolated_nodes = list(networkx.isolates(G_centroid))
    new_pos = 1
    for node in isolated_nodes:
        centroid_locations[node] = np.array([min_x - 0.25, new_pos])
        new_pos -= 0.1

    return centroid_locations

def draw_centroid_plot(G_centroid, centroid_locations, communities_dict, centroid_ax, color_of_comunities):
    networkx.draw_networkx_nodes(G_centroid, pos=centroid_locations, ax=centroid_ax, node_size = 0.01)

    ## Draw intercommunity edges
    weights = [G_centroid[u][v]['pos_weight'] for u,v in G_centroid.edges()]
    weights = np.append(weights, [G_centroid[u][v]['neg_weight'] for u,v in G_centroid.edges()])
    max_weight = max(weights)
    min_weight = min(weights)
    for weight_type, color in zip(['pos_weight', 'neg_weight'], ['red', 'blue']):
        weights = [G_centroid[u][v][weight_type] for u,v in G_centroid.edges()]
        max_width = 1
        min_width = 0.01
        edge_widths = [min_width + (w - min_weight) / (max_weight - min_weight) * (max_width - min_width) for w in weights]
        edge_widths = [weight*5 for weight in edge_widths]
        networkx.draw_networkx_edges(G_centroid, pos=centroid_locations, edge_color=color, width=edge_widths, ax=centroid_ax, connectionstyle="arc3,rad=0.1" if weight_type == 'pos_weight' else "arc3,rad=-0.1", arrows = True)

    ## Draw community labels
    for node in G_centroid.nodes():
        centroid_ax.annotate(str(node), xy=centroid_locations[node], textcoords='offset points', xytext=(0,0), ha='center', va='center', fontsize=7, color='white', bbox=dict(boxstyle="circle,pad=1", fc=color_of_comunities[node], ec='black', lw=0.5))

    ## Add title
    centroid_ax.set_title('Community centroids', fontsize=10)

def create_gene_network_subset(communities_dict, G, community_0 = 0, community_1 = 1, k_network = 0.1, seed = 0, positive_edge_weight_cutoff = 0.1, negative_edge_weight_cutoff = -0.1, recalculate_position = True):
    if community_0 == community_1:
        genes_to_plot = communities_dict[community_0]
    else:
        genes_to_plot = np.append(communities_dict[community_0], communities_dict[community_1])
    G_subset = G.subgraph(genes_to_plot)

    ## Subset genes to plot based on edge weights
    positive_edges_to_plot = np.array([(u,v) for u, v, d in G_subset.edges(data = True) if d['weight'] > positive_edge_weight_cutoff])
    negative_edges_to_plot = np.array([(u,v) for u, v, d in G_subset.edges(data = True) if d['weight'] < negative_edge_weight_cutoff])
    exists_positive_edges_to_plot = positive_edges_to_plot.shape[0] > 0
    exists_negative_edges_to_plot = negative_edges_to_plot.shape[0] > 0

    if exists_positive_edges_to_plot and exists_negative_edges_to_plot:
        edge_list = np.concatenate((positive_edges_to_plot, negative_edges_to_plot)).tolist()
    elif exists_positive_edges_to_plot:
        edge_list = positive_edges_to_plot.tolist()
    elif exists_negative_edges_to_plot:
        edge_list = negative_edges_to_plot.tolist()

    edge_list = [(u,v) for u,v in edge_list]

    G_subset_edges = G_subset.edge_subgraph(edge_list).copy()
    isolated_nodes = list(networkx.isolates(G_subset_edges))
    G_subset_edges.remove_nodes_from(isolated_nodes)

    if recalculate_position == True:
        position = networkx.spring_layout(G_subset_edges, k = k_network, seed = seed, method = 'force') 
    else:
        position = networkx.spring_layout(G_subset, k = k_network, seed = seed, method = 'force') 

    return G_subset_edges, position, positive_edges_to_plot, negative_edges_to_plot

def draw_network_plot(G_subset_edges, community_0, community_1, communities_dict, position, positive_edges_to_plot, negative_edges_to_plot, positive_edge_weight_cutoff, negative_edge_weight_cutoff, network_ax, color_of_comunities, font_size = 6):
    ## Subset genes to plot based on communities

    ## Draw edges
    networkx.draw_networkx_edges(G_subset_edges, pos=position, edgelist=positive_edges_to_plot, edge_color='red', width=0.25, ax=network_ax)
    networkx.draw_networkx_edges(G_subset_edges, pos=position, edgelist=negative_edges_to_plot, edge_color='blue', width=0.25, ax=network_ax)

    ## Draw nodes
    networkx.draw_networkx_nodes(G_subset_edges, pos=position, ax=network_ax, node_size=0.1)

    ## Draw labels
    for community in [community_0, community_1]:
        genes_of_community = communities_dict[community]
        labels = {i:i for i in genes_of_community if i in G_subset_edges.nodes()}
        networkx.draw_networkx_labels(G_subset_edges, pos=position, labels=labels, bbox=dict(boxstyle="round,pad=0.2", fc=color_of_comunities[community], ec='black', lw=0, alpha=0.5), ax = network_ax, font_size = font_size, font_color = 'black')

    ## Add title
    if community_0 == community_1:
        text = f'Community {community_0}'
    else:
        text = f'Community {community_0} and {community_1}'
    network_ax.set_title(f'{text} network \n (positive edges > {np.round(positive_edge_weight_cutoff, 5)}) \n (negative edges < {np.round(negative_edge_weight_cutoff, 5)})', fontsize=10)

def make_communities_dict(var_with_rank, n_top_genes = 100):
    communities_dict = {}
    for community in range(var_with_rank['Community'].nunique()):
        community_df = var_with_rank[var_with_rank['Community'] == community].sort_values(by = 'Centrality', ascending = False)
        if community_df.shape[0] < 5:
            break
        else:
            top_genes = community_df[:n_top_genes].index
            temp_dict = {community:top_genes}
            communities_dict.update(temp_dict)
    
    return communities_dict

def plot_network(G, var_with_rank, k_centroids=0.1, k_network=0.1, seed = 0, axes = None, n_top_genes = 100, positive_edge_weight_cutoff = 0.1, negative_edge_weight_cutoff = -0.1, community_0 = 0, community_1 = 0, calculate_spring = True, recalculate_position = True):

    centroid_ax, network_ax = axes

    # Make centroid plot
    ## Subset genes based on communities and centrality
    communities_dict = make_communities_dict(var_with_rank, n_top_genes = n_top_genes)

    ## Build centroid graph
    G_centroid, df = build_centroid_graph(G, communities_dict)
    centroid_locations = calculate_centroid_layout(G_centroid, k_centroids=k_centroids, seed = seed, calculate_spring = calculate_spring)

    ## Draw centroid plot
    ## Assign colors to communities
    color_of_comunities = {community:color for community, color in zip(communities_dict.keys(), sns.color_palette("hls", len(communities_dict)))}
    draw_centroid_plot(G_centroid, centroid_locations, communities_dict, centroid_ax, color_of_comunities)

    # Create network
    G_subset_edges, position, positive_edges_to_plot, negative_edges_to_plot = create_gene_network_subset(communities_dict, G, community_0 = community_0, community_1 = community_1, k_network = k_network, seed = seed, positive_edge_weight_cutoff = positive_edge_weight_cutoff, negative_edge_weight_cutoff = negative_edge_weight_cutoff, recalculate_position = recalculate_position)
    
    draw_network_plot(G_subset_edges, community_0, community_1, communities_dict, position, positive_edges_to_plot, negative_edges_to_plot, positive_edge_weight_cutoff, negative_edge_weight_cutoff, network_ax, color_of_comunities)

    return communities_dict

def update(G, var_with_rank, k_centroids=0.1, k_network=0.1, n_top_genes = 100, positive_edge_weight_cutoff = 0.1, negative_edge_weight_cutoff = -0.1, community_0 = 0, community_1 = 0, calculate_spring = True, recalculate_position = False):
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        ax.clear()
    communities_dict = plot_network(G, var_with_rank, k_centroids=k_centroids, k_network=k_network, seed = 0, axes = axes, n_top_genes = n_top_genes, positive_edge_weight_cutoff = positive_edge_weight_cutoff, community_0 = community_0, community_1 = community_1, negative_edge_weight_cutoff = negative_edge_weight_cutoff, calculate_spring = calculate_spring, recalculate_position = recalculate_position)
    fig.canvas.draw_idle()  # Efficient redraw