import streamlit as st
import age
from pyvis.network import Network
from streamlit_agraph import agraph, Node, Edge, ConfigBuilder,Config
import psycopg2
from age import networkx
import networkx as nx


st.set_page_config(
    layout="wide",
)


# connect DB
@st.cache_data
def get_graph(list_gene, depth):
    graphName = "kg2"
    conn = psycopg2.connect(
        host="192.168.2.131",
        port="5432",
        dbname="ksi_cpds",
        user="postgres",
        # password="123456"
        )
    gsql = (
        f"""SELECT * from cypher('%s', $$ MATCH p=(n )-[r*{depth}]-() 
        where n.__id__ in {list_gene} 
      
      
        RETURN  p   $$) as (v agtype)"""
        % graphName
    )
    age.setUpAge(conn, graphName)

    G = networkx.age_to_networkx(conn, graphName, query=gsql)
    return G




col1, col2 = st.columns(2)
with col1:
    depth = st.slider("Select depth:", 0, 4, 1)
    deg = st.slider("Select Degree:", 0, 10, 0)
    list_rel=['protein_protein', 'drug_protein', 'drug_drug','disease_protein', 'disease_disease','pathway_pathway', 'pathway_protein']
    sel_rel = st.multiselect("Chose relations for the graph", list_rel)
with col2:
    
    var_text = st.text_area(
        "Enter your list of entities", help="Name or ID separated by enter"
    )
    var_t = var_text.split("\n")
    list_gene = [t.strip().upper() for t in var_t if t != ""]



if len(list_gene) > 0:
    G=get_graph(list_gene, depth)
    
    edge_remove = set([d for d in G.edges(data="label", default=1) if  d[-1] not in (sel_rel)])
    G.remove_edges_from(edge_remove)
    
    remove = [x for x in G.nodes() if G.degree(x) < deg]
    G.remove_nodes_from(remove)
    
    
    
   
 
    
    label_dict = dict(G.nodes(data="properties", default=1))
    converted_dict = {key: value["__id__"] for key, value in label_dict.items()}
    H = nx.relabel_nodes(G, converted_dict)

    net = Network(notebook=False)
    net.from_nx(H)

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
    config = Config(height=1000,
                    width=1000,
                    nodeHighlightBehavior=True,
                    highlightColor="#F7A7A6",
                    directed=True,
                    collapsible=True,
                    physics=False, 
                    staticGraphWithDragAndDrop=True,
                    link={'labelProperty': 'label', 'renderLabel': True}
                    )

   
    return_value = agraph(nodes, edges, config=config)
