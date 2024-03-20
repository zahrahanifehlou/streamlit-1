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


RHSA1 = nx.read_gml("../R-HSA-8949215.gml")
G = RHSA1
# lab = G.nodes
# for i in lab:
#     st.write(lab[i]["type"])
# data = nx.node_link_data(RHSA1)
# st.write(data)
# for node, nbrsdict in RHSA1.adj.items():
#     st.write(node, nbrsdict)
# G = nx.karate_club_graph()
# test = nx.to_dict_of_dicts(G)
# df_test = pd.DataFrame.from_dict(test)
# st.write(test)
adjmat = nx.adjacency_matrix(G).todense()
adjmat = pd.DataFrame(
    index=range(adjmat.shape[0]), data=adjmat, columns=range(adjmat.shape[0])
)
adjmat.columns = adjmat.columns.astype(str)
adjmat.index = adjmat.index.astype(str)
# adjmat.iloc[3, 4] = 5
# adjmat.iloc[4, 5] = 6
# adjmat.iloc[5, 6] = 7

df = pd.DataFrame(index=adjmat.index)

df["degree"] = np.array([*G.degree()])[:, 1]
df["label"] = [i for i in (G.nodes)]
df["relation"] = [G.nodes[i]["type"] for i in G.nodes]
# st.write(adjmat)
st.write(df)
# exit(0)
label = df["label"].values
node_size = df["degree"].values
# node_size
# adjmat = nx.to_pandas_adjacency(RHSA1)
edge_df = nx.to_pandas_edgelist(RHSA1)
st.write(edge_df)
st.write(adjmat)
adjmat = vec2adjmat(edge_df["source"], edge_df["target"])
d3 = d3graph()

d3.graph(adjmat)
d3.set_node_properties(color=df["relation"].values, label=label)
# d3.set_node_properties()
d3.set_edge_properties(directed=True)

for idx, row in edge_df.iterrows():
    d3.edge_properties[row["source"], row["target"]]["label"] = row["type"]

d3.show()


d3 = d3graph()
adjmat = d3.import_example("bigbang")
st.write(adjmat)
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
