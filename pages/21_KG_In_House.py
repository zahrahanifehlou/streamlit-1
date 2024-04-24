import streamlit as st
import pandas as pd
import networkx as nx
import streamlit.components.v1 as components
from pyvis.network import Network
import igraph as ig
from streamlit_agraph import agraph, Node, Edge, Config
import matplotlib
from sklearn import preprocessing
from streamlib import sql_df
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np


def draw_graph(KG: pd.DataFrame, deg: int, html: str):
    # GG = nx.Graph()
    GG = nx.from_pandas_edgelist(KG)
    st.write(GG.nodes())

    get_nodes = GG.nodes()

    keep = [node for node, degree in dict(GG.degree()).items() if degree > deg]

    remove = get_nodes - keep

    GG.remove_nodes_from(remove)

    net = Network(notebook=True)
    net = net.from_nx(GG)
    st.write(net)
    # net.show_buttons(filter_="physics")
    net.show(f"/mnt/shares/L/Temp/{html}.html", notebook=False)

    HtmlFile = open(f"/mnt/shares/L/Temp/{html}.html", "r", encoding="utf-8")
    source_code = HtmlFile.read()
    components.html(source_code, height=1200, width=1200)


choix_source = "Ksilink_25"
profile_conn = "postgres://arno:12345@192.168.2.131:5432/ksilink_cpds"
conn = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"
sql_profile = f"select * from aggprofile where metasource='{choix_source}'"
df_source = sql_df(sql_profile, profile_conn)
df_source = df_source.rename(columns={"metabatchid": "batchid"})
sql_batch = f"select * from platemap where source='{choix_source}'"
df_batchs = sql_df(sql_batch, conn)
df_batchs = df_batchs.drop_duplicates()
df_merge = df_source.merge(df_batchs, on="batchid")
st.write(df_merge.sample(5))
df_merge = df_merge.rename(
    columns={"batchid": "meta_batchid", "pubchemid": "meta_pubchemid"}
)
df_merge = df_merge.dropna(subset="meta_pubchemid")
# st.write(df_merge)
filter_col1 = [col for col in df_merge.columns if not col.startswith("meta")]

simi = cosine_similarity(df_merge[filter_col1])
df_merge["source"] = df_merge["meta_pubchemid"]
df_merge["target"] = df_merge["meta_pubchemid"]

simi_df_alt = pd.DataFrame(simi, index=df_merge["source"], columns=df_merge["target"])
simi_df_alt = simi_df_alt.where(np.triu(np.ones(simi_df_alt.shape), k=1).astype(bool))
simi_df_alt = simi_df_alt.stack().reset_index()
simi_df_alt.columns = ["source", "target", "similarity"]


# simi_high = simi_df_alt[simi_df_alt["similarity"] > 0.9]
# st.write(simi_high.reset_index(drop=True))

G = nx.from_pandas_edgelist(simi_df_alt)
st.write(G)
# node_and_degree = G.degree()
# most_connected_node = sorted(G.degree, key=lambda x: x[1], reverse=True)[0]
# degree = G.degree(most_connected_node)

# # Create egonet for the focal node
# hub_ego = nx.ego_graph(G, most_connected_node[0])

# Create the equivalent Node and Edge lists
nodes = [Node(id=i, label=str(i), size=20) for i in (G.nodes)]
edges = [
    Edge(source=i, target=j, type="CURVE_SMOOTH")
    for (i, j) in G.edges
    # if i in hub_ego.nodes and j in hub_ego.nodes
]


config = Config(
    width=1200,
    height=1200,
    directed=False,
    nodeHighlightBehavior=True,
    highlightColor="#F7A7A6",
    collapsible=True,
    node={"labelProperty": "label"},
    link={"labelProperty": "label", "renderLabel": True},
)

return_value = agraph(nodes=nodes, edges=edges, config=config)
# # st.write(GG.nodes())

# # get_nodes = GG.nodes()

# # keep = [node for node, degree in dict(GG.degree()).items() if degree > 0]

# # remove = get_nodes - keep

# # GG.remove_nodes_from(remove)

# net = Network(notebook=False)
# net = net.from_nx(GG)

# html = "tests"
# # net.show_buttons(filter_="physics")
# net.show(f"/mnt/shares/L/Temp/{html}.html", notebook=False)

# HtmlFile = open(f"/mnt/shares/L/Temp/{html}.html", "r", encoding="utf-8")
# source_code = HtmlFile.read()
# components.html(source_code, height=1200, width=1200)

# # draw_graph(simi_high, 0, "cpds")
