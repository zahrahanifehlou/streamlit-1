import streamlit as st
import age
from pyvis.network import Network
from streamlit_agraph import agraph, Node, Edge, ConfigBuilder
from streamlit_agraph import  Config as Config2
import psycopg2
from age import networkx
import networkx as nx
from time import process_time

from nebula3.Config import Config
from nebula3.gclient.net import ConnectionPool
st.set_page_config(
    layout="wide",
)




@st.cache_data
def get_graph_postgres(list_gene, depth,sel_rel):
    graphName = "graph_kg"
    conn = psycopg2.connect(
        host="192.168.2.131",
        port="5432",
        dbname="ksi_cpds",
        user="postgres",
        # password="123456"
        )
    if len(sel_rel)==0:
        sel_rel=['protein_protein' ]
    gsql = (
        f"""SELECT * from cypher('%s', $$ MATCH p=(n )-[r*{depth}]-() 
        where n.__id__ in {list_gene} 
        and type(relationships(p)[0] ) in {sel_rel}
        RETURN  p   $$) as (v agtype)"""
        % graphName
    )
    st.write(gsql)
    age.setUpAge(conn, graphName)

    G = networkx.age_to_networkx(conn, graphName, query=gsql)
    return G

# @st.cache_data
def get_graph_nebula(list_gene, depth,sel_rel):
    
    client = None
    space_name ="KG"
    config = Config()
    config.max_connection_pool_size = 2
    connection_pool = ConnectionPool()
    assert connection_pool.init([("192.168.2.130", 9669)], config) 
    # Get one gclient
    client = connection_pool.get_session("root", "nebula")
    assert client is not None
    resp = client.execute(
        "CREATE SPACE IF NOT EXISTS {} (vid_type=FIXED_STRING(30)); USE {};".format(
            space_name, space_name
        )
    )
    #resp =client.execute(space_name)
    assert resp.is_succeeded(), resp.error_msg()


    if len(sel_rel)==0:
        sel_rel=['protein_protein' ]
    gsql = (
        f"""MATCH p=(n )-[r*{depth}]-() 
        where n.protein.name in {list_gene} 
        and type(relationships(p)[0] ) in {sel_rel}
        RETURN  p"""
    
    )
    st.write(gsql)
    resp = client.execute(gsql)
    data_for_vis = resp.dict_for_vis()
  
    return data_for_vis


col1, col2 = st.columns(2)
with col1:
    depth = st.slider("Select depth:", 0, 4, 1)
    deg = st.slider("Minimum Degree:", 0, 10,1)
    
with col2:
    list_rel=['protein_protein', 'drug_protein','bioprocess_bioprocess',  'bioprocess_protein',  'molfunc_molfunc', 'molfunc_protein', 
      'drug_drug','disease_protein', 'disease_disease','pathway_pathway', 'pathway_protein']
    sel_rel = st.multiselect("Chose relations for the graph", list_rel,default ='protein_protein')
    
    var_text = st.text_area(
        "Enter your list of entities", help="Name or ID separated by enter"
    )
    var_t = var_text.split("\n")
    list_gene = [t.strip().upper() for t in var_t if t != ""]

postgres_graph, nebula_graph = st.tabs(
        ["postgres graph", "nebula graph  "]
    )
with postgres_graph:
    st.write("postgres graph")
    # if len(list_gene) > 0:
    #     G=get_graph_postgres(list_gene, depth,sel_rel)
        
        
        
    #     remove = [x for x in G.nodes() if G.degree(x) < deg]
    #     G.remove_nodes_from(remove)
        
    #     label_dict = dict(G.nodes(data="properties", default=1))
    #     converted_dict = {key: value["__id__"] for key, value in label_dict.items()}
    #     H = nx.relabel_nodes(G, converted_dict)

    #     net = Network(notebook=False)
    #     net.from_nx(H)

    #     nodes = [
    #         Node(
    #             id=u["id"],
    #             label=u["id"],
    #             color=u["color"],
    #             # title=list_desc[i],
    #             shape=u["shape"],
    #             font="10px arial grey",
    #             size=15,
    #         )
    #         for u in net.nodes
    #     ]

    #     edges = [
    #         Edge(
    #             source=v["from"],
    #             target=v["to"],
    #             # title=v["title"],
    #             # color=v["color"],
    #         )
    #         for v in net.edges
    #     ]
    #     config2 = Config2(height=600,
    #                     width=1000,
    #                     nodeHighlightBehavior=True,
    #                     highlightColor="#F7A7A6",
    #                     directed=True,
    #                     collapsible=True,
    #                     physics=False, 
    #                     staticGraphWithDragAndDrop=True,
    #                     link={'labelProperty': 'label', 'renderLabel': True}
    #                     )

    
    #     return_value = agraph(nodes, edges, config=config2)
    #     cols=st.columns(2)
    #     with(cols[0]):
    #         st.write("node list\n", H.nodes)
    #     with (cols[1]):
    #         st.write("edge list \n",H.edges)
        

with nebula_graph:
    if len(list_gene) > 0:
        t1_start = process_time()  
        data_for_vis=get_graph_nebula(list_gene, depth,sel_rel)
        t1_stop = process_time() 
        
      ##  "GO 2 STEP FROM '10013' OVER protein_protein REVERSELY  YIELD src(edge), dst(edge)"
        st.write("Elapsed time during the whole program in seconds:", 
                                                t1_stop-t1_start) 
        nodes=data_for_vis.get("nodes_dict")
        edges=data_for_vis.get("edges_dict")
        new_edge_list=[]
        for e in edges:
            e_list = e.strip("()").split(", ")
            first_num = e_list[0].strip("'")
            second_num = e_list[1].strip("'")
            result = (f'{first_num}', f'{second_num}')
            new_edge_list.append(result)
        
        color_dict={'protein':'green' , 'drug':'blue' ,'pathway':'yellow' ,'disease':'red' }
        node_color = [node['labels'][0] for node in nodes.values()]
        node_color = list(map(color_dict.get, node_color))

        G = nx.from_edgelist(new_edge_list)
        nx.draw(G,node_color=node_color)

        
        remove = [x for x in G.nodes() if G.degree(x) < deg]
        G.remove_nodes_from(remove)
    

        
        node_names = [node['props']['name'] for node in nodes.values()]
        node_id = [node['id'] for node in nodes.values()]
        node_dict = dict(zip( node_id,node_names))
        node_color = [node['labels'][0] for node in nodes.values()]
        

        H = nx.relabel_nodes(G, node_dict)

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
        

        # nodes = [
        #     Node(
        #         id=u["id"],
        #         label=u['props']['name'],
        #         #color=u["color"],
        #         # title=list_desc[i],
        #     # shape=u["shape"],
        #         font="10px arial grey",
        #         size=15,
        #     )
        #     for u in nodes.values()

        # ]

        # edges = [
        #     Edge(
        #         source=v["src"],
        #         target=v["dst"],
        #         # title=v["title"],
        #         # color=v["color"],
        #     )
      
        #      for v in edges.values()
        # ]
        config2 = Config2(height=600,
                        width=1000,
                        nodeHighlightBehavior=True,
                        highlightColor="#F7A7A6",
                        directed=True,
                        collapsible=True,
                        physics=False, 
                        staticGraphWithDragAndDrop=True,
                        link={'labelProperty': 'label', 'renderLabel': True}
                        )

    
        return_value = agraph(nodes, edges, config=config2)
        cols=st.columns(2)
        with(cols[0]):
            st.write("node list\n", H.nodes)
        with (cols[1]):
            st.write("edge list \n",H.edges)
