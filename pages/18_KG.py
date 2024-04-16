import streamlit as st
import pandas as pd
import networkx as nx
import streamlit.components.v1 as components
from pyvis.network import Network
import igraph as ig
from streamlit_agraph import agraph, Node, Edge, Config, ConfigBuilder
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


KG2 = load_data()
st.write("Test", KG2.sample(22).reset_index())


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
    st.write("Selected KG", KG)

else:
    exit(0)
# exit(0)
graph_list = []


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
            list_nodes = list(set(list_node + list_node2))
            # st.write(list_nodes)
            # st.write(list_gene)
            cpt_sel = 0
            cpt_other = 0
            for i, name in enumerate(list_nodes):
                # st.write(i)
                # st.write(name.upper())
                # if name.upper() in list(set(merge) - set(list_gene)):
                #     # st.write(name)
                #     net.add_node(name, label=name, color="blue")
                if name.upper() in list_gene:
                    net.add_node(name, label=name, color="grey")
                    cpt_sel = cpt_sel + 1
                else:
                    k = int(df.loc[df["source"] == name, "categorical_label"].unique())
                    # st.write(tuple(k))
                    cpt_other = cpt_other + 1
                    net.add_node(
                        name,
                        label=name,
                        color=matplotlib.colors.to_hex(palette.get(k), keep_alpha=True),
                    )
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
            nodes = [
                Node(
                    id=u["id"],
                    label=u["label"],
                    color=u["color"],
                    # title=list_desc[i],
                    shape=u["shape"],
                    font="10px arial grey",
                    size=15,
                )
                for u in net.nodes
            ]

            edges = [
                Edge(
                    source=v["from"],
                    target=v["to"],
                    title=v["title"],
                    color=v["color"],
                )
                for v in net.edges
            ]

            st.warning(f"Retrieved {cpt_sel} members; Added {cpt_other} interactions ")
            c = ConfigBuilder(nodes, edges)
            d = c.build()

            return_value = agraph(nodes, edges, config=d)
