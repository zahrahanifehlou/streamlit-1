import streamlit as st
import pandas as pd
import networkx as nx
import streamlit.components.v1 as components
from pyvis.network import Network

st.set_page_config(
    layout="wide",
)
st.title("Knowledge Graph (PrimeKG)")
st.header("17080 diseases, 4M relationship,100K Nodes,20 Sources")
st.subheader(
    "Building a knowledge graph to enable precision medicine (Harvard MS)",
    divider="rainbow",
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
# st.write("Test", KG[(KG.x_name.str.contains("Tyrphostins"))])


# list_test = ["Risperidone", "autism spectrum disorder"]
# KGG = KG.query('x_name =="autism spectrum disorder" and y_name=="Risperidone"')
# st.write(KGG)

# tutu = KG.query('x_name == "autism spectrum disorder" ')
# st.write(tutu)

# x = (
#     KG.query('x_name == "autism spectrum disorder" ')
#     .get(["y_index"])
#     .drop_duplicates()
#     .values
# )
# x = x.squeeze()

# toto = KG.query('x_index in @x and y_name=="Risperidone"')
# st.write("toto", toto)
# titi = KG.query(
#     'x_name == "autism spectrum disorder" and y_index in [3199, 3417, 33913] '
# )
# st.write("titi", titi)
var_text = st.text_area(
    "Enter your list of entities", help="Name or ID separated by enter"
)
var_t = var_text.split("\n")
list_gene = [t.strip().upper() for t in var_t]
list_gene = var_t

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


if var_text:
    #     GG = nx.Graph()
    #     netgg = Network(notebook=False)
    #     GG.add_edges_from(KG[["x_name", "y_name"]].values)
    #     netgg.from_nx(GG)
    #     netgg.show_buttons(filter_=["physics"])
    #     netgg.show("total.html", notebook=False)

    #     HtmlFile2 = open("total.html", "r", encoding="utf-8")
    #     source_code2 = HtmlFile2.read()

    #     components.html(source_code2, height=1200, width=1200)
    # U = nx.Graph()
    for i in KG.relation.unique():
        G = nx.Graph()
        net = Network(notebook=False)
        G.add_edges_from(KG[KG.relation == i][["x_name", "y_name"]].values, relation=i)
        # if i == "indication":
        #     for node, degree in dict(G.degree()).items():
        #         st.write("node,deg", (node, degree))
        # remove = [node for node, degree in dict(G.degree()).items() if degree == 1]
        # G.remove_nodes_from(remove)
        # U.add_edges_from(G.edges(data=True))
        net.from_nx(G)
        net.show_buttons(filter_=["physics"])
        net.show(f"/mnt/shares/L/Temp/{i}.html", notebook=False)

        HtmlFile = open(f"/mnt/shares/L/Temp/{i}.html", "r", encoding="utf-8")
        source_code = HtmlFile.read()
        st.warning(i)
        components.html(source_code, height=1200, width=1200)

        G.clear()
# net = Network(notebook=False)
# net.from_nx(U)
# net.show_buttons(filter_=["physics"])
# net.show("total.html", notebook=False)

# HtmlFile = open("total.html", "r", encoding="utf-8")
# source_code = HtmlFile.read()

# components.html(source_code, height=1200, width=1200)
