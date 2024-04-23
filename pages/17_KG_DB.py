import streamlit as st
import age
from pyvis.network import Network
from streamlit_agraph import agraph, Node, Edge, ConfigBuilder
import psycopg2
from age import networkx
import networkx as nx


st.set_page_config(
    layout="wide",
)


# connect DB
graphName = "kg1"
conn = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    dbname="ksi_cpds",
    user="postgres",
    # password="123456"
)

age.setUpAge(conn, graphName)


var_text = st.text_area(
    "Enter your list of entities", help="Name or ID separated by enter"
)
var_t = var_text.split("\n")
deg = st.slider("Select Degree:", 0, 10, 2)
list_gene = [t.strip().upper() for t in var_t if t != ""]
if len(list_gene) > 0:
    gsql = (
        f"""SELECT * from cypher('%s', $$ MATCH p=(n )<-[r*{deg}]-() 
        where n.__id__ in {list_gene} 
        RETURN  p   $$) as (v agtype)"""
        % graphName
    )
    st.write(gsql)
    G = networkx.age_to_networkx(conn, graphName, query=gsql)
    label_dict = dict(G.nodes(data="properties", default=1))
    converted_dict = {key: value["__id__"] for key, value in label_dict.items()}
    H = nx.relabel_nodes(G, converted_dict)
    # H = nx.relabel_nodes(G, converted_dict)
    pos = nx.spring_layout(H)
    nx.draw_networkx_labels(H, pos)
    nx.draw_networkx_edges(H, pos, edge_color="r", arrows=False)
    # nx.draw_networkx_nodes(H, pos)

    net = Network(notebook=False)
    net.from_nx(H)
    # st.write(net.nodes)
    nodes = [
        Node(
            id=u["id"],
            label=u["id"],
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
            # title=v["title"],
            # color=v["color"],
        )
        for v in net.edges
    ]

    c = ConfigBuilder(
        nodes,
        edges,
    )
    d = c.build(dictify=False)
    return_value = agraph(nodes, edges, config=d)
