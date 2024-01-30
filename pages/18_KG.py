import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import streamlit.components.v1 as components
from pyvis.network import Network

st.set_page_config(
    layout="wide",
)


@st.cache_data
def load_data():
    primekg = pd.read_csv("../kg.csv", low_memory=False)
    return primekg


vert_size = st.sidebar.slider(
    "vertex size", min_value=0.2, max_value=20.0, step=0.1, value=1.0
)
lab_size = st.sidebar.slider(
    "label size", min_value=0.2, max_value=20.0, step=0.1, value=2.4
)

KG = load_data()

var_text = st.text_area(
    "Enter your list of entities", help="Name or ID separated by enter"
)
var_t = var_text.split("\n")
list_gene = [t.strip().upper() for t in var_t]


KG = KG.query("x_name == @list_gene | y_name==@list_gene")

graph_list = []

options = {
    "node_color": "blue",
    "node_size": vert_size,
    "width": 1,
    "font_size": lab_size,
    "edge_color": "green",
    "alpha": 0.6,
}
for i in KG.relation.unique():
    G = nx.Graph()
    net = Network(notebook=False)
    G.add_edges_from(KG[KG.relation == i][["x_name", "y_name"]].values, relation=i)

    net.from_nx(G)
    net.show(f"{i}.html", notebook=False)
    HtmlFile = open(f"{i}.html", "r", encoding="utf-8")
    source_code = HtmlFile.read()
    components.html(source_code, height=900, width=1200)

    G.clear()
