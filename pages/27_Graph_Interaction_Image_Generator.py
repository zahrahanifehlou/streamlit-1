import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from io import BytesIO

# WARNING : Traceback unable
st.set_option('deprecation.showPyplotGlobalUse', False)

# Load nodes and edges
@st.cache_data

def load_data():
    nodes_df = pd.read_csv('/mnt/shares/L/DB/Moussa/VS_Code/nodes.csv')
    edges_df = pd.read_csv('/mnt/shares/L/DB/Moussa/VS_Code/edges.csv')
    relations_df = pd.read_csv('/mnt/shares/L/DB/Moussa/VS_Code/kg.csv')
    return nodes_df, edges_df, relations_df

nodes_df, edges_df, relations_df = load_data()

st.title('Graph Interaction Image Generator')
st.write('Welcome to our graphical interaction image generator. Please enter gene/protein names, pathways, drugs, effects/phenotypes, anatomies, diseases, exposures, biological processes, molecular functions, or cellular components, separated by commas in the dialogue box below. Our software will then display all desired interactions related to the selected elements in an image that you can download and keep. The depth corresponds to the number of secondary interactions displayed for each neighbor of the selected elements. Feel free to also use the live Viewer as needed!')

def plot_graph(nodes_list, interaction_types_list, neighborhood_depth):
    # Create a graph
    G = nx.Graph()
    
    # Map node indices to names and types for faster lookup
    node_idx_to_name = nodes_df.set_index('node_index')['node_name'].to_dict()
    node_idx_to_type = nodes_df.set_index('node_index')['node_type'].to_dict()
    
    # Setup node colors and shapes based on types
    node_types = nodes_df['node_type'].unique()
    colors = plt.cm.get_cmap('tab20', len(node_types))
    node_color_map = {ntype: colors(i) for i, ntype in enumerate(node_types)}
    shapes = ['o', 's', 'D', '^', 'p', 'h', '*']  # Circle, square, diamond, triangle, pentagon, hexagon, star
    node_shape_map = {ntype: shapes[i % len(shapes)] for i, ntype in enumerate(node_types)}
    
    # Setup edge colors based on selected interaction types
    interaction_types = relations_df['display_relation'].unique()
    edge_colors = plt.cm.get_cmap('tab10', len(interaction_types))
    edge_color_map = {etype: edge_colors(i) for i, etype in enumerate(interaction_types)}
    
    # Filter nodes and relations based on input lists and depth
    relevant_nodes = set(nodes_list)
    current_layer = set(nodes_list)
    
    for _ in range(neighborhood_depth):
        next_layer = set()
        for node in current_layer:
            node_idx = nodes_df[nodes_df['node_name'] == node]['node_index'].iloc[0]
            
            # Find all edges that involve this node and are of the selected types
            edges = relations_df[((relations_df['x_name'] == node) | (relations_df['y_name'] == node)) &
                                 (relations_df['display_relation'].isin(interaction_types_list))]
            
            for _, row in edges.iterrows():
                x_name = row['x_name']
                y_name = row['y_name']
                relation = row['display_relation']
                
                # Add nodes and edges to the graph if the interaction type is selected
                if relation in interaction_types_list:
                    x_type = node_idx_to_type[row['x_index']]
                    y_type = node_idx_to_type[row['y_index']]
                    G.add_node(x_name, type=x_type)
                    G.add_node(y_name, type=y_type)
                    G.add_edge(x_name, y_name, interaction=relation)
                    
                    # Determine the next set of nodes to process
                    if x_name == node:
                        next_layer.add(y_name)
                    else:
                        next_layer.add(x_name)
                        
        current_layer = next_layer
        relevant_nodes.update(current_layer)
    
    # Drawing the graph with a more sophisticated layout
    pos = nx.kamada_kawai_layout(G)  # Kamada-Kawai layout for better visual separation
    plt.figure(figsize=(16, 16))  # Increase the figure size
    
    # Draw nodes with type-specific styles
    for ntype, shape in node_shape_map.items():
        nx.draw_networkx_nodes(G, pos, nodelist=[n for n in G.nodes if G.nodes[n]['type'] == ntype],
                               node_color=[node_color_map[ntype]],
                               node_shape=shape, node_size=1000)  # Reduced node size
    
    # Draw edges with interaction-specific colors
    nx.draw_networkx_edges(G, pos, edge_color=[edge_color_map[G[u][v]['interaction']] for u, v in G.edges() if G[u][v]['interaction'] in interaction_types_list])
    
    # Node labels in bold and even smaller size
    nx.draw_networkx_labels(G, pos, font_size=6)
    
    # Legend for node types and edge types
    legend_elements = [plt.Line2D([0], [0], marker=shape, color='w', label=f'{ntype} nodes',
                                  markerfacecolor=node_color_map[ntype], markersize=10)
                       for ntype, shape in node_shape_map.items()]
    legend_elements += [plt.Line2D([0], [0], color=edge_color_map[etype], lw=4, label=f'{etype}')
                        for etype in interaction_types_list if etype in edge_color_map]
    
    plt.legend(handles=legend_elements, loc='upper left')
    
    plt.show()
    st.pyplot()

# Step-by-step input
node_input = st.text_input('Enter nodes separated by commas (no space)')
nodes_list = node_input.split(',')

if st.button('Confirm Nodes'):
    st.session_state.nodes_confirmed = True

if 'nodes_confirmed' in st.session_state:
    interaction_types = relations_df['display_relation'].unique()
    selected_interactions = st.multiselect('Select interaction types', options=interaction_types)

    if selected_interactions:
        st.session_state.interactions_confirmed = True

if 'interactions_confirmed' in st.session_state:
    depth = st.slider('Select neighborhood depth', min_value=1, max_value=5, value=1)
    if st.button('Plot Graph'):
        plot_graph(nodes_list, selected_interactions, depth)
        st.session_state.plot_confirmed = True