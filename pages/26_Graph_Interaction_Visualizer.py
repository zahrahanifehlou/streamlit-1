import streamlit as st
import pandas as pd
import networkx as nx
from streamlit_agraph import agraph, Node, Edge, Config

@st.cache_data

def load_data():
    nodes_df = pd.read_csv('/mnt/shares/L/DB/Moussa/VS_Code/nodes.csv')
    relations_df = pd.read_csv('/mnt/shares/L/DB/Moussa/VS_Code/kg.csv')
    return nodes_df, relations_df

nodes_df, relations_df = load_data()

def plot_interactive_graph(nodes_list, interaction_types_list, neighborhood_depth, graph_height, graph_width, static_graph, show_labels):
    G = nx.Graph()
    node_idx_to_type = nodes_df.set_index('node_name')['node_type'].to_dict()
    colors = {
        'gene/protein': '#FF5733', 
        'drug': '#33FFCE', 
        'effect/phenotype': '#335BFF', 
        'disease': '#9BFF33', 
        'biological_process': '#F833FF', 
        'molecular_function': '#FFC300', 
        'cellular_component': '#DAF7A6', 
        'exposure': '#FF5733', 
        'pathway': '#C70039', 
        'anatomy': '#900C3F'
    }
    shapes = {
        'gene/protein': 'square',
        'pathway': 'square',
        'drug': 'star',
        'effect/phenotype': 'triangle',
        'anatomy': 'triangle',
        'disease': 'diamond',
        'exposure': 'diamond',
        'biological_process': 'hexagon',
        'molecular_function': 'hexagon',
        'cellular_component': 'hexagon'
    }
    interaction_colors = {
        'ppi': '#FF5733',
        'interacts with': '#577399',  # Slate blue gray
        'carrier': '#8D8741',  # Olive green
        'enzyme': '#659DBD',  # Earthy blue
        'target': '#FBEEC1',  # Light yellow
        'transporter': '#BC986A',  # Tan brown
        'contraindication': '#5D5C61',  # Grayish
        'indication': '#7395AE',  # Gentle blue
        'off-label use': '#B1A296',  # Dusty rose
        'synergistic interaction': '#6B5E62',  # Taupe
        'associated with': '#A4C3B2',  # Celadon
        'parent-child': '#B9F8D3',  # Soft mint
        'phenotype absent': '#E0BBE4',  # Lavender floral
        'phenotype present': '#957DAD',  # Pastel violet
        'side effect': '#D291BC',  # Pastel pink
        'linked to': '#FEC8D8',  # Tea rose
        'expression absent': '#6F9FD8',  # Sky blue
        'expression present': '#D6D4E0',  # Light periwinkle
}
    
    for node in nodes_list:
        if node in node_idx_to_type:
            G.add_node(node, type=node_idx_to_type[node])
            edges = relations_df[((relations_df['x_name'] == node) | (relations_df['y_name'] == node)) &
                                 (relations_df['display_relation'].isin(interaction_types_list))]
            for _, row in edges.iterrows():
                G.add_node(row['x_name'], type=node_idx_to_type.get(row['x_name'], 'Unknown'))
                G.add_node(row['y_name'], type=node_idx_to_type.get(row['y_name'], 'Unknown'))
                G.add_edge(row['x_name'], row['y_name'], interaction=row['display_relation'])

    nodes = [Node(id=n, label=n if show_labels else '', size=10, color=colors[G.nodes[n]['type']], shape=shapes[G.nodes[n]['type']], font={'size': 8}) for n in G]
    edges = [Edge(source=u, target=v, color=interaction_colors[G[u][v]['interaction']], label='') for u, v in G.edges()]
    
    config = Config(
        height=graph_height,
        width=graph_width,
        directed=False,
        nodeHighlightBehavior=True,
        highlightColor="#F7A7A6",
        node={'labelProperty': 'label' if show_labels else None},
        staticGraph=static_graph
    )
    agraph(nodes=nodes, edges=edges, config=config)
    
    st.markdown("### Legend")
    for ntype, color in colors.items():
        st.markdown(f"{ntype}: <span style='color:{color}'>■</span> {shapes[ntype]}", unsafe_allow_html=True)
    
    st.markdown("#### Interaction Types")
    for interaction_type, color in interaction_colors.items():
        st.markdown(f"{interaction_type}: <span style='color:{color}'>■</span>", unsafe_allow_html=True)
        
    st.subheader("Nodes in Graph")
    st.write("Here you will find the set of nodes displayed in the current graph.")
    st.table(G.nodes(data=True))

st.title('Graph Interaction Visualizer')
st.write('Welcome to our interaction viewer. Please enter gene/protein names, pathways, drugs, effects/phenotypes, anatomical features, diseases, exposures, biological processes, molecular functions, or cellular components, separated by commas in the dialogue box below. Our software will then display all desired interactions related to the elements that have been chosen. The depth corresponds to the number of secondary interactions displayed for each neighbor of the selected elements. Feel free to use the sidebar to make the graph more readable according to your needs. A image generator is also available! Dont hesitate to use it as needed!')

st.sidebar.title('Graph Configuration')
graph_height = st.sidebar.number_input('Height of the graph', min_value=400, max_value=1000, value=600)
graph_width = st.sidebar.number_input('Width of the graph', min_value=400, max_value=1000, value=800)
static_graph = st.sidebar.checkbox('Static Graph', value=False)
show_labels = st.sidebar.checkbox('Display Node Names', value=True)

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
        plot_interactive_graph(nodes_list, selected_interactions, depth, graph_height, graph_width, static_graph, show_labels)
