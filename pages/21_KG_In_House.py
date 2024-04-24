import streamlit as st
import pandas as pd
import networkx as nx
import streamlit.components.v1 as components
from streamlib import sql_df
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np



choix_source = "Ksilink_25"
conn = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"
sql_query="SELECT aggprofile.*, cpdbatchs.pubchemid , cpd.cpdname , cpd.smile , cpd.keggid  FROM aggprofile\
    INNER JOIN cpdbatchs ON cpdbatchs.batchid = aggprofile.batchid left JOIN cpd ON cpdbatchs.pubchemid = cpd.pubchemid WHERE aggprofile.source = 'Ksilink_25'"

df_merge = sql_df(sql_query, conn)
filter_col1 = df_merge.select_dtypes(include=[int, float]).columns.tolist()
st.write(len(df_merge))



simi = cosine_similarity(df_merge[filter_col1])

simi_df_alt = pd.DataFrame(simi, index=df_merge["source"], columns=df_merge["batchid"])
simi_df_alt = simi_df_alt.where(np.triu(np.ones(simi_df_alt.shape), k=1).astype(bool))
simi_df_alt = simi_df_alt.stack().reset_index()
simi_df_alt.columns = ["source", "batchid", "similarity"]
st.write(len(simi_df_alt))

# # simi_high = simi_df_alt[simi_df_alt["similarity"] > 0.9]
# # st.write(simi_high.reset_index(drop=True))

# G = nx.from_pandas_edgelist(simi_df_alt)
# st.write(G)
# # node_and_degree = G.degree()
# # most_connected_node = sorted(G.degree, key=lambda x: x[1], reverse=True)[0]
# # degree = G.degree(most_connected_node)

# # # Create egonet for the focal node
# # hub_ego = nx.ego_graph(G, most_connected_node[0])

# # Create the equivalent Node and Edge lists
# nodes = [Node(id=i, label=str(i), size=20) for i in (G.nodes)]
# edges = [
#     Edge(source=i, target=j, type="CURVE_SMOOTH")
#     for (i, j) in G.edges
#     # if i in hub_ego.nodes and j in hub_ego.nodes
# ]


# config = Config(
#     width=1200,
#     height=1200,
#     directed=False,
#     nodeHighlightBehavior=True,
#     highlightColor="#F7A7A6",
#     collapsible=True,
#     node={"labelProperty": "label"},
#     link={"labelProperty": "label", "renderLabel": True},
# )

# return_value = agraph(nodes=nodes, edges=edges, config=config)
# # # st.write(GG.nodes())

# # # get_nodes = GG.nodes()

# # # keep = [node for node, degree in dict(GG.degree()).items() if degree > 0]

# # # remove = get_nodes - keep

# # # GG.remove_nodes_from(remove)

# # net = Network(notebook=False)
# # net = net.from_nx(GG)

# # html = "tests"
# # # net.show_buttons(filter_="physics")
# # net.show(f"/mnt/shares/L/Temp/{html}.html", notebook=False)

# # HtmlFile = open(f"/mnt/shares/L/Temp/{html}.html", "r", encoding="utf-8")
# # source_code = HtmlFile.read()
# # components.html(source_code, height=1200, width=1200)

# # # draw_graph(simi_high, 0, "cpds")
