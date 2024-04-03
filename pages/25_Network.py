import pandas as pd
import requests
from Bio.KEGG import REST
from bs4 import BeautifulSoup
import streamlit as st
from streamlit_d3graph import d3graph, vec2adjmat
import networkx as nx
import numpy as np
import graphistry

st.set_page_config(layout="wide")

st.title("Load Network Graph")
st.header("KEGG")
st.subheader(
    "Here you are",
    divider="rainbow",
)
url = "http://rest.genome.jp/link/hsa:6622/original"
response = requests.get(url)
soup = BeautifulSoup(response.content, "html.parser")
path_net = soup.get_text().split("\n")

st.write(path_net)
exit(0)


# url1 = "http://rest.kegg.jp/link/pathway/network"
# url2 = "http://rest.kegg.jp/link/disease/network"
# url3 = "http://rest.kegg.jp/link/hsa/network"


# def get_infos(url):
#     df = pd.DataFrame()
#     response = requests.get(url)
#     # Parse the response using BeautifulSoup
#     soup = BeautifulSoup(response.content, "html.parser")
#     path_net = soup.get_text().split("\n")
#     list1 = []
#     list2 = []
#     list3 = []
#     for t in path_net:
#         g = t.split("\t")
#         if len(g) == 2:
#             list1.append(g[0])
#             list2.append(g[1])
#             if "ds:" in g[1]:
#                 list3.append("disease")
#             if "hsa:" in g[1]:
#                 list3.append("gene")
#             if "path:" in g[1]:
#                 list3.append("pathway")

#     df["source"] = list1
#     df["target"] = list2
#     df["relation"] = list3
#     return df


# df_url1 = get_infos(url1)
# df_url2 = get_infos(url2)
# df_url3 = get_infos(url3)
# df_tot = pd.concat([df_url1, df_url2, df_url3])
df_tot = pd.read_csv("/mnt/shares/L/PROJECTS/JUMP-CRISPR/KeggNetworks/keggnet.csv")
st.write(df_tot)
# df_graph = pd.DataFrame()
# list_source = df_tot.source.unique()

# list_types = df_tot.relation.unique()

# with st.container(border=True):
#     col1, col2 = st.columns(2)
#     # source = col1.multiselect("Select Sources", list_source)
#     types = col1.selectbox("Select Types", list_types)
#     if types:
#         df_tot = df_tot[df_tot["relation"] == types]
#         list_target = df_tot.target.unique()
#         target = col2.selectbox("Select targets", list_target)
#         df_tot = df_tot[df_tot["target"] == target]
#         df_graph = df_tot
#     # types = col2.multiselect("Select Types", list_types)
#
if not df_tot.empty:
    G = nx.DiGraph()

    # # st.write(KG)
    # for s, t, r in df_graph[["source", "target", "relation"]].values:
    #     G.add_edge(s, t, label=r)

    # # G = nx.from_pandas_adjacency(df_tot)
    # st.write(G)

    # adjmat = nx.adjacency_matrix(G).todense()
    # adjmat = pd.DataFrame(
    #     index=range(adjmat.shape[0]), data=adjmat, columns=range(adjmat.shape[0])
    # )
    # adjmat.columns = adjmat.columns.astype(str)
    # adjmat.index = adjmat.index.astype(str)

    # df = pd.DataFrame(index=adjmat.index)

    # df["degree"] = np.array([*G.degree()])[:, 1]
    # df["label"] = [i for i in (G.nodes)]
    # # df["relation"] = [G.nodes[i]["relation"] for i in G.nodes]
    # st.write(df)
    # label = df["label"].values
    # node_size = df["degree"].values
    df_tot = df_tot.sample(1000)
    adjmat = vec2adjmat(df_tot["source"], df_tot["target"])
    d3 = d3graph()

    d3.graph(adjmat)
    # d3.set_node_properties(color=df["relation"].values, label=label)
    # d3.set_node_properties()
    d3.set_edge_properties(directed=True)
    d3.show(title="Metabolics Flux", figsize=(1200, 1200))
