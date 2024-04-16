import streamlit as st
import pandas as pd
import networkx as nx
import streamlit.components.v1 as components
from pyvis.network import Network
import igraph as ig
from streamlit_agraph import agraph, Node, Edge, Config
import matplotlib
from sklearn import preprocessing

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

KG2 = load_data()
st.write("Test", KG2.sample(22).reset_index())


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
list_gene = [t.strip().upper() for t in var_t if t != ""]
# list_gene = var_t
# st.write(list_gene)
# exit(0)
pattern = "|".join(list_gene)
# st.write(pattern)
# KG = KG.query("x_name == @list_gene | y_name==@list_gene")
# KG = KG.query("x_name in @list_gene")
if var_text:
    # KG = KG[KG["x_name"].str.contains(pattern, case=False, regex=True)]
    KG = KG2.query("x_name == @list_gene | y_name==@list_gene")
    # for i in range()
    # KG2 = KG[KG["y_name"].str.contains(pattern, case=False, regex=True)]
    # KG = pd.concat([KG1, KG2])
    # KG = KG[KG["x_name"].isin(list_gene)]
    # st.write(KG)
    # list_y = KG["y_name"].to_list()
    # KG = KG3.query("relation!='anatomy_protein_present'")

    # KG = KG.query("relation!='anatomy_protein_present'")
    st.write(len(KG.loc[KG["relation"] == "protein_protein", "x_name"].unique()))
    # st.write(KG.loc[KG["relation"] == "protein_protein"])

else:
    exit(0)
# exit(0)
graph_list = []
# st.write(KG3.reset_index(drop=True))
# st.write(KG.reset_index(drop=True))
# exit(0)
options = {
    "node_color": "blue",
    "node_size": vert_size,
    "width": 1,
    "font_size": lab_size,
    "edge_color": "green",
    "alpha": 0.6,
}

nodes = []
if var_text:
    GG = nx.DiGraph()
    list_rel = KG.relation.unique()

    sel_rel = st.multiselect("Chose relations for the graph", list_rel)
    KG = KG[KG["relation"].isin(sel_rel)]
    st.write("KG1", len(KG))

    # merge = list(set(KG["x_name"].unique().tolist() + KG["y_name"].unique().tolist()))
    # KG = KG2.query("x_name in @merge | y_name in @merge")
    # st.write(len(KG.loc[KG["relation"]=="protein_protein","x_name"].unique()))
    # sel_rel2 = st.multiselect("Chose second relations for the graph", list_rel)
    # KG = KG[KG["relation"].isin(sel_rel2)]
    # st.write("KG",len(KG))
    for s, t, r in KG[["x_name", "y_name", "relation"]].values:
        GG.add_edge(s, t, label=r)
    # GG.add_edges_from(KG[["x_name", "y_name"]].values, label=KG["relation"].values)
    deg = st.slider("Select Degree:", 0, 10, 1)
    get_nodes = GG.nodes()
    keep = [node for node, degree in dict(GG.degree()).items() if degree > deg]
    # st.write(get_nodes)
    remove = get_nodes - keep
    # st.write(remove)
    GG.remove_nodes_from(remove)
    df = nx.to_pandas_edgelist(GG)
    # st.write(df)
    # G = ig.Graph.from_networkx(GG)
    # communities = G.community_edge_betweenness()
    # communities = communities.as_clustering()
    num_communities = len(KG.relation.unique())
    communities = KG.relation.unique()

    if not df.empty:
        # st.warning(f"{sel_rel}")
        if sel_rel:
            palette = ig.RainbowPalette(n=num_communities)
            # df = nx.to_pandas_edgelist(GG)
            # st.write(df)
            le = preprocessing.LabelEncoder()
            le.fit(df.label)
            df["categorical_label"] = le.transform(df.label)
            #
            # st.write("DF", df)
            net = Network(notebook=False, heading=sel_rel)
            list_node = list(df["source"].unique())
            list_node2 = list(df["target"].unique())
            list_nodes = list_node + list_node2
            # st.write(list_nodes)
            # st.write(list_gene)
            for i, name in enumerate(list_nodes):
                # st.write(i)
                # st.write(name.upper())
                # if name.upper() in list(set(merge) - set(list_gene)):
                #     # st.write(name)
                #     net.add_node(name, label=name, color="blue")
                if name.upper() in list_gene:
                    net.add_node(name, label=name, color="green")
                else:
                    net.add_node(name, label=name, color="blue")
            # net.add_nodes(list_nodes)
            # st.write(net)
            for i, j, k, m in df[
                ["source", "target", "categorical_label", "label"]
            ].values:
                # st.write(i)
                net.add_edge(
                    i,
                    j,
                    color=matplotlib.colors.to_hex(palette.get(k), keep_alpha=True),
                    title=m,
                )
            # couleur = matplotlib.colors.to_hex(palette.get(i), keep_alpha=True)
            # st.write(couleur)
            # st.write(u)
            # st.write(data)

            #     nodes.append(
            #         Node(
            #             id=i,
            #             label=u,r
            #             # color=couleur,
            #             title=u,
            #             # shape="circle",
            #             size=20,
            #         )
            #     )
            # net.from_nx(GG)
            net.show_buttons(filter_="physics")
            net.show("/mnt/shares/L/Temp/test4.html", notebook=False)

            HtmlFile = open("/mnt/shares/L/Temp/test4.html", "r", encoding="utf-8")
            source_code = HtmlFile.read()
            components.html(source_code, height=1200, width=1200)

        # A = GG.edges()
        # st.write(A)
        # st.write(nodes)
        # edges = [Edge(source=i, target=j, type="CURVE_SMOOTH") for (i, j) in A]
        # st.write(edges[19])
        # config = Config(
        #     width=1200,
        #     height=1200,
        #     directed=False,
        #     nodeHighlightBehavior=True,
        #     highlightColor="#F7A7A6",
        #     collapsible=True,
        #     node={"labelProperty": "label"},
        #     physics=True,
        #     # link={"labelProperty": "label", "renderLabel": True},
        # )
        # st.write("### Displaying the full network/graph clustered with Leiden approach")
        # return_value = agraph(nodes=nodes, edges=edges, config=config)
    #     netgg.from_nx(GG)
    #     netgg.show_buttons(filter_=["physics"])
    #     netgg.show("total.html", notebook=False)

    #     HtmlFile2 = open("total.html", "r", encoding="utf-8")
    #     source_code2 = HtmlFile2.read()

    #     components.html(source_code2, height=1200, width=1200)
    # U = nx.Graph()
    # G.clear()

    # for r in KG.relation.unique():
    #     G = nx.Graph()
    #     nodes = []
    #     net = Network(notebook=False)
    #     G.add_edges_from(KG[KG.relation == r][["x_name", "y_name"]].values, relation=r)
    #     # if i == "indication":
    #     #     for node, degree in dict(G.degree()).items():
    #     #         st.write("node,deg", (node, degree))
    #     remove = [node for node, degree in dict(G.degree()).items() if degree == 1]
    #     G.remove_nodes_from(remove)
    #     G = ig.Graph.from_networkx(G)
    #     communities = G.community_edge_betweenness()
    #     communities = communities.as_clustering()
    #     num_communities = len(communities)
    #     if num_communities > 0:
    #         palette = ig.RainbowPalette(n=num_communities)

    #         for i, g in enumerate(communities):
    #             # U.add_edges_from(G.edges(data=True))
    #             couleur = matplotlib.colors.to_hex(palette.get(i), keep_alpha=True)
    #             for u in G.vs[g]:
    #                 # st.write(couleur)
    #                 nodes.append(
    #                     Node(
    #                         id=u.index,
    #                         label=u["_nx_name"],
    #                         color=couleur,
    #                         title=u["_nx_name"],
    #                         # shape="circle",
    #                         size=20,
    #                     )
    #                 )
    #         # net.from_nx(G)
    #         # net.show_buttons(filter_=["physics"])
    #         # net.show(f"/mnt/shares/L/Temp/{i}.html", notebook=False)

    #         # HtmlFile = open(f"/mnt/shares/L/Temp/{i}.html", "r", encoding="utf-8")
    #         # source_code = HtmlFile.read()
    #         st.warning(r)
    #         A = G.get_edgelist()
    #         edges = [Edge(source=i, target=j, type="CURVE_SMOOTH") for (i, j) in A]
    #         config = Config(
    #             width=1200,
    #             height=1200,
    #             directed=False,
    #             nodeHighlightBehavior=True,
    #             highlightColor="#F7A7A6",
    #             collapsible=True,
    #             node={"labelProperty": "label"},
    #             physics=True,
    #             # link={"labelProperty": "label", "renderLabel": True},
    #         )
    #         st.write(
    #             "### Displaying the full network/graph clustered with Leiden approach"
    #         )
    #         return_value = agraph(nodes=nodes, edges=edges, config=config)
    #         # components.html(source_code, height=1200, width=1200)

    #         G.clear()
# net = Network(notebook=False)
# net.from_nx(U)
# net.show_buttons(filter_=["physics"])
# net.show("total.html", notebook=False)

# HtmlFile = open("total.html", "r", encoding="utf-8")
# source_code = HtmlFile.read()

# components.html(source_code, height=1200, width=1200)
