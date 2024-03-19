import networkx as nx
import igraph as ig
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from pyvis.network import Network
import streamlit.components.v1 as components

from streamlit_agraph import agraph, Node, Edge, Config

# from streamlit_visgraph import visgraph

st.set_page_config(layout="wide")


RHSA1 = nx.read_gml("../R-HSA-1430728.gml")
G = ig.Graph.from_networkx(RHSA1)
A = [edge.tuple for edge in G.es]
# In case your graph is directed
H = nx.DiGraph(A)
# H = G.to_nx()
st.write(RHSA1)

ig.config["plotting.backend"] = "matplotlib"
subax1 = ig.plot(
    G,
    vertex_label=G.vs["_nx_name"],
    vertex_label_size=2,
    vertex_size=10,
    # margin=10,
    # bbox=(0, 0, box_size, box_size),
)
# nx.draw(GG, with_labels=True,node_size=10,font_size=4)
st.pyplot(subax1.figure, use_container_width=False)


# st.title("Streamlit VisGraph - Game of Thrones example")

# got_data = pd.read_csv("https://www.macalester.edu/~abeverid/data/stormofswords.csv")

# sources = got_data["Source"]
# targets = got_data["Target"]
# weights = got_data["Weight"]
# nodes = []
# edges = []
# # node_config = NodeConfig(shape="dot")
# # edge_config = EdgeConfig()
# options = Config()
# edge_data = zip(sources, targets, weights)
# nodes_all = sources.tolist() + targets.tolist()
# node_data = list(set(nodes_all))
# for i in range(0, len(node_data)):
#     nodes.append(
#         Node(
#             id=i,
#             label=node_data[i],
#             title=node_data[i],
#             value=nodes_all.count(node_data[i]),
#             url="http://example/" + node_data[i],
#             # node_config=node_config,
#         )
#     )

# for e in edge_data:
#     src = e[0]
#     dst = e[1]
#     w = e[2]
#     edges.append(
#         Edge(
#             source=node_data.index(src),
#             target=node_data.index(dst),
#             # edge_config=edge_config,
#         )
#     )

# v = visgraph(nodes=nodes, edges=edges, config=options)


nodes = []
# thres_db = col_b.slider("cluster Thresholds", 2, 100,9,1)
net = Network(notebook=False, directed=True)

net.from_nx(RHSA1)
net.show_buttons(filter_="physics")
net.show("/mnt/shares/L/Temp/arno2.html", notebook=False)
HtmlFile = open("/mnt/shares/L/Temp/arno2.html", "r", encoding="utf-8")
source_code = HtmlFile.read()


components.html(source_code, height=1200, width=1200)


# for u in G.vs["_nx_name"]:
#     # st.write(u)
#     nodes.append(
#         Node(
#             id=u.index,
#             label=u,
#             shape="circle",
#             size=200,
#         )
#     )


# A = G.get_edgelist()

# edges = [Edge(source=i, target=j) for (i, j) in RHSA1.edges()]
# st.write(edges[0])
# config = Config(
#     width=1200,
#     height=1200,
#     directed=True,
#     # nodeHighlightBehavior=True,
#     # highlightColor="#F7A7A6",
#     # collapsible=True,
#     # node={"labelProperty": "label"},
#     hierarchical=False,
#     # physics=True,
#     # link={"labelProperty": "label", "renderLabel": True},
# )
# # tedges = RHSA1.edges()
# # st.write(tedges)
# return_value = agraph(nodes=nodes, edges=edges, config=config)

# nodes.append(
#     Node(
#         id="Spiderman",
#         label="Peter Parker",
#         size=25,
#         shape="circularImage",
#         image="http://marvel-force-chart.surge.sh/marvel_force_chart_img/top_spiderman.png",
#     )
# )  # includes **kwargs
# nodes.append(
#     Node(
#         id="Captain_Marvel",
#         size=25,
#         shape="circularImage",
#         image="http://marvel-force-chart.surge.sh/marvel_force_chart_img/top_captainmarvel.png",
#     )
# )
# edges.append(
#     Edge(
#         source="Captain_Marvel",
#         label="friend_of",
#         target="Spiderman",
#         # **kwargs
#     )
# )

# config = Config(
#     width=750,
#     height=950,
#     directed=True,
#     physics=True,
#     hierarchical=False,
#     # **kwargs
# )

# return_value = agraph(nodes=nodes, edges=edges, config=config)
