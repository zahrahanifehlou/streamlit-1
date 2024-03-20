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


RHSA1 = nx.read_gml("../R-HSA-70171.gml")
# data = nx.node_link_data(RHSA1)
# st.write(data)
# for node, nbrsdict in RHSA1.adj.items():
#     st.write(node, nbrsdict)

adjmat = nx.to_pandas_adjacency(RHSA1)
edge_df = nx.to_pandas_edgelist(RHSA1)
# st.write(edge_df)
# st.write(adjmat)
d3 = d3graph()

d3.graph(adjmat)
# d3.set_node_properties(
#     label=df["label"].values,
#     color=df["label"].values,
#     size=df["degree"].values,
#     edge_size=df["degree"].values,
#     cmap="Set1",
# )
d3.set_edge_properties(directed=True)
d3.show()


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
