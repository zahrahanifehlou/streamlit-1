import networkx as nx
import igraph as ig
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from pyvis.network import Network
import streamlit.components.v1 as components

from streamlit_agraph import agraph, Node, Edge, Config
# from d3graph import d3graph, vec2adjmat

from streamlit_d3graph import d3graph, vec2adjmat

# from streamlit_visgraph import visgraph

st.set_page_config(layout="wide")

st.title("Load Metabolics Graph")
st.header("Reactomes, KEGG")
st.subheader(
    "Here you are",
    divider="rainbow",
)


uploaded_file = st.file_uploader(
    "Choose gml file",
    accept_multiple_files=False,
)
edge_df = pd.DataFrame()
if uploaded_file:
    RHSA1 = nx.read_gml(uploaded_file)
    edge_df = nx.to_pandas_edgelist(RHSA1)
    list_source = edge_df.source.unique()
    list_target = edge_df.target.unique()
    list_types = edge_df.type.unique()

    with st.container(border=True):
        col1, col2, col3 = st.columns(3)
        source = col1.multiselect("Select Sources", list_source)
        target = col2.multiselect("Select targets", list_target)
        types = col3.multiselect("Select Types", list_types)

if not edge_df.empty:
    G = RHSA1

    adjmat = nx.adjacency_matrix(G).todense()
    adjmat = pd.DataFrame(
        index=range(adjmat.shape[0]), data=adjmat, columns=range(adjmat.shape[0])
    )
    adjmat.columns = adjmat.columns.astype(str)
    adjmat.index = adjmat.index.astype(str)

    df = pd.DataFrame(index=adjmat.index)

    df["degree"] = np.array([*G.degree()])[:, 1]
    df["label"] = [i for i in (G.nodes)]
    df["relation"] = [G.nodes[i]["type"] for i in G.nodes]

    label = df["label"].values
    node_size = df["degree"].values

    adjmat = vec2adjmat(edge_df["source"], edge_df["target"])
    d3 = d3graph()

    d3.graph(adjmat)
    d3.set_node_properties(color=df["relation"].values, label=label)
    # d3.set_node_properties()
    d3.set_edge_properties(directed=True)

    for idx, row in edge_df.iterrows():
        d3.edge_properties[row["source"], row["target"]]["label"] = row["type"]

    d3.show(title="Metabolics Flux", figsize=(1200, 1200))


# # Initialize with clustering colors
# d3.graph(adjmat, color="cluster")
# GG = nx.DiGraph()

# st.write(KG)
# for s, t, r in edge_df[["source", "target", "type"]].values:
#     GG.add_edge(s, t, label=r)

# nodes = [Node(id=i, label=str(i), size=35, shape="dot") for i in GG.nodes]
# edges = [
#     Edge(source=i, target=j, label=k, type="CURVE_SMOOTH")
#     for (i, j, k) in edge_df[["source", "target", "type"]].values
# ]
# config = Config(
#     width=1200,
#     height=1200,
#     directed=True,
#     nodeHighlightBehavior=True,
#     highlightColor="#F7A7A6",
#     collapsible=True,
#     physics=True,
#     hierarchical=False,
#     node={"labelProperty": "label"},
#     link={"labelProperty": "label", "renderLabel": True},
# )
# return_value = agraph(nodes=nodes, edges=edges, config=config)
