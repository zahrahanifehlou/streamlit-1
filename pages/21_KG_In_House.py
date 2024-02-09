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
    GG = nx.Graph()

    for s, t in KG[["x_name", "y_name"]].values:
        GG.add_edge(s, t)
    # GG.add_edges_from()
    # deg = st.slider("Select Degree:", 0, 10, 1)
    get_nodes = GG.nodes()
    keep = [node for node, degree in dict(GG.degree()).items() if degree > deg]

    remove = get_nodes - keep

    GG.remove_nodes_from(remove)
    df = nx.to_pandas_edgelist(GG)

    if not df.empty:
        net = Network(notebook=False)
        net = net.from_nx(GG)

        net.show_buttons(filter_="physics")
        net.show(f"/mnt/shares/L/Temp/{html}.html", notebook=False)

        HtmlFile = open(f"/mnt/shares/L/Temp/{html}.html", "r", encoding="utf-8")
        source_code = HtmlFile.read()
        components.html(source_code, height=1200, width=1200)


choix_source = "Ksilink_25"
profile_conn = "postgres://arno:12345@192.168.2.131:5432/ksilink_cpds"
conn = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"
sql_profile = f"select * from aggprofile where metasource='{choix_source}'"
df_source = sql_df(sql_profile, profile_conn)
sql_batch = f"select * from cpdbatchs where metasource='{choix_source}'"
df_batchs = sql_df(sql_batch, conn)
df_batchs = df_batchs.drop_duplicates()
df_merge = df_source.merge(df_batchs, left_on="metabatchid", right_on="batchid")
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
simi_df_alt.columns = ["x_name", "y_name", "relation"]


simi_high = simi_df_alt[simi_df_alt["relation"] > 0.95]
st.write(simi_high.reset_index(drop=True))
draw_graph(simi_df_alt, 3, "cpds")
