import streamlit as st
import age
from pyvis.network import Network
from streamlit_agraph import agraph, Node, Edge, ConfigBuilder
from streamlit_agraph import  Config as Config2
import psycopg2
from age import networkx
import networkx as nx
from time import process_time
import pandas as pd
from nebula3.Config import Config
from nebula3.gclient.net import ConnectionPool
st.set_page_config(
    layout="wide",
)

from nebula3.data.DataObject import Value, ValueWrapper
from nebula3.data.ResultSet import ResultSet




def cast(val: ValueWrapper):
    cast_as = {
    Value.NVAL: "as_null",
    Value.BVAL: "as_bool",
    Value.IVAL: "as_int",
    Value.FVAL: "as_double",
    Value.SVAL: "as_string",
    Value.LVAL: "as_list",
    Value.UVAL: "as_set",
    Value.MVAL: "as_map",
    Value.TVAL: "as_time",
    Value.DVAL: "as_date",
    Value.DTVAL: "as_datetime",
    Value.VVAL: "as_node",
    Value.EVAL: "as_relationship",
    Value.PVAL: "as_path",
    Value.GGVAL: "as_geography",
    Value.DUVAL: "as_duration",
        }
    _type = val._value.getType()
    if _type == Value.__EMPTY__:
        return None
    if _type in cast_as:
        return getattr(val, cast_as[_type])()
    if _type == Value.LVAL:
        return [x.cast() for x in val.as_list()]
    if _type == Value.UVAL:
        return {x.cast() for x in val.as_set()}
    if _type == Value.MVAL:
        return {k: v.cast() for k, v in val.as_map().items()}

def resp_to_dataframe(resp: ResultSet):
    assert resp.is_succeeded()
    output_table=[]
    for recode in resp:
        value_list = []
        for col in recode:
            val = cast(col)
            value_list.append(val)
        output_table.append(value_list)
    output_table=pd.DataFrame(output_table,columns=resp.keys())
    return(output_table)



@st.cache_data
def get_graph_nebula_match(list_gene, depth,sel_rel):
    
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
  
    assert resp.is_succeeded(), resp.error_msg()

    
    if len(sel_rel)==0:
        sel_rel=['protein_protein' ]
    gsql = (
        f"""MATCH (n )
        where n.protein.name in {list_gene} or n.drug.name in {list_gene} or n.pathway.name in {list_gene} or n.disease.name in {list_gene} or n.cellcomp.name in {list_gene}
        RETURN   id(n)"""
    )
    
    #st.write(gsql)
    
    t1_start = process_time()
    resp = client.execute(gsql)
    t1_stop = process_time()
    st.write(" time to get nodes in seconds:",  t1_stop-t1_start)  
    
    id_list=[val for val in resp]
    output=pd.DataFrame()
    if len(id_list) >0:
        id_list = ', '.join(f'{item}' for item in id_list)
        #sel_rel = ', '.join(f'{item}' for item in sel_rel)
        
    
    
    
        gsql = (
            f"""MATCH p=(n )-[r*{depth}]-(m) 
            where id(n) in [{id_list} ]
            and   id(m) in [{id_list} ] and  type(relationships(p)[0] ) in {sel_rel} RETURN  p"""
        )
        st.write(gsql)

        resp = client.execute(gsql)
        
        # gen name of nodes and edges
        data_for_vis = resp.dict_for_vis()
        nodes=data_for_vis.get("nodes_dict")
        edges=data_for_vis.get("edges_dict")
        new_edge_list=[]
        for e in edges:
            e_list = e.strip("()").split(", ")
            first_num = e_list[0].strip("'")
            second_num = e_list[1].strip("'")
            type = e_list[-1].strip("'")
            
            new_edge_list.append([first_num, second_num,type ])
        output_table=pd.DataFrame(new_edge_list,columns=["source","target","type"])
        node_names = [node['props']['name'] for node in nodes.values() if 'name' in node['props']]
        node_id = [node['id'] for node in nodes.values()]

        node_dict = dict(zip( node_id,node_names))

        output_table["source"]=output_table["source"].map(node_dict)
        output_table["target"]=output_table["target"].map(node_dict)

    return output_table


# @st.cache_data
# def get_graph_nebula_Go(list_gene, depth,sel_rel):
   
#     client = None
#     space_name ="KG"
#     config = Config()
#     config.max_connection_pool_size = 2
#     connection_pool = ConnectionPool()
#     assert connection_pool.init([("192.168.2.130", 9669)], config) 
#     # Get one gclient
#     client = connection_pool.get_session("root", "nebula")
#     assert client is not None
#     resp = client.execute(
#         "CREATE SPACE IF NOT EXISTS {} (vid_type=FIXED_STRING(30)); USE {};".format(
#             space_name, space_name
#         )
#     )
#     #resp =client.execute(space_name)
#     assert resp.is_succeeded(), resp.error_msg()

    
#     if len(sel_rel)==0:
#         sel_rel=['protein_protein' ]
#     gsql = (
#         f"""MATCH (n )
#         where n.protein.name in {list_gene} or n.drug.name in {list_gene} or n.pathway.name in {list_gene} or n.disease.name in {list_gene} or n.cellcomp.name in {list_gene}
#         RETURN   id(n)"""
#     )
    
#     #st.write(gsql)
    
#     t1_start = process_time()
#     resp = client.execute(gsql)
#     t1_stop = process_time()
#     st.write(" time to get nodes in seconds:",  t1_stop-t1_start)  
    
#     id_list=[val for val in resp]
#     output=pd.DataFrame()
#     if len(id_list) >0:
#         id_list = ', '.join(f'{item}' for item in id_list)
#         sel_rel = ', '.join(f'{item}' for item in sel_rel)
        
    
#         gsql2=f'''GO {depth} STEP FROM  {id_list} OVER {sel_rel}   YIELD   properties($^).name as source , properties($$).name as target , type(edge) as type;'''
#         st.write(gsql2)
        
#         t2_start = process_time()
#         resp = client.execute(gsql2)
#         t2_stop = process_time()
#         st.write(" time to get edges in seconds:",  t2_stop-t2_start) 
        
#         assert resp.is_succeeded(), resp.error_msg()
#         output=resp_to_dataframe(resp)
#     #st.write(output)
#     return output

col1, col2 = st.columns(2)
with col1:
    depth = st.slider("Select depth:", 0, 4, 1)
    deg = st.slider("Minimum Degree:", 0, 10,1)
    
with col2:
    list_rel=['protein_protein', 'drug_protein',
              #'bioprocess_bioprocess',  'bioprocess_protein',  'molfunc_molfunc', 'molfunc_protein', 
      'drug_drug','disease_protein', 'disease_disease','pathway_pathway', 'pathway_protein', 'cellcomp_cellcomp', 'cellcomp_protein']
    sel_rel = st.multiselect("Chose relations for the graph", list_rel,default ='protein_protein')
    
    var_text = st.text_area(
        "Enter your list of entities", help="Name or ID separated by enter"
    )
    var_t = var_text.split("\n")
    list_search = [t.strip().lower() for t in var_t if t != ""]



if len(list_search) > 0:
  
    graph_df=get_graph_nebula_match(list_search, depth,sel_rel)
    st.write(graph_df)
    
    if len(graph_df)>0:
        graph_df.dropna(inplace=True)
        graph_df['source'] = graph_df['source'].str.split(' ').str[0]
        graph_df['target'] = graph_df['target'].str.split(' ').str[0]
    
 
        t2_start = process_time() 
        H = nx.from_pandas_edgelist(graph_df)
        t2_stop = process_time()
        st.write(" time to draw graph in seconds:", t2_stop-t2_start)  
        
        
        remove = [x for x in H.nodes() if H.degree(x) < deg]
        H.remove_nodes_from(remove)

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






#----------------------------------------------------------------

# @st.cache_data
# def get_graph_postgres(list_gene, depth,sel_rel):
#     graphName = "graph_kg"
#     conn = psycopg2.connect(
#         host="192.168.2.131",
#         port="5432",
#         dbname="ksi_cpds",
#         user="postgres",
#         # password="123456"
#         )
#     if len(sel_rel)==0:
#         sel_rel=['protein_protein' ]
#     gsql = (
#         f"""SELECT * from cypher('%s', $$ MATCH p=(n )-[r*{depth}]-() 
#         where n.__id__ in {list_gene} 
#         and type(relationships(p)[0] ) in {sel_rel}
#         RETURN  p   $$) as (v agtype)"""
#         % graphName
#     )
#     st.write(gsql)
#     age.setUpAge(conn, graphName)

#     G = networkx.age_to_networkx(conn, graphName, query=gsql)
#     return G

# postgres_graph, nebula_graph = st.tabs(
#         ["postgres graph", "nebula graph  "]
#     )
# with postgres_graph:
#     st.write("postgres graph")
#     # if len(list_gene) > 0:
#     #     G=get_graph_postgres(list_gene, depth,sel_rel)
        
        
        
#     #     remove = [x for x in G.nodes() if G.degree(x) < deg]
#     #     G.remove_nodes_from(remove)
        
#     #     label_dict = dict(G.nodes(data="properties", default=1))
#     #     converted_dict = {key: value["__id__"] for key, value in label_dict.items()}
#     #     H = nx.relabel_nodes(G, converted_dict)

#     #     net = Network(notebook=False)
#     #     net.from_nx(H)

#     #     nodes = [
#     #         Node(
#     #             id=u["id"],
#     #             label=u["id"],
#     #             color=u["color"],
#     #             # title=list_desc[i],
#     #             shape=u["shape"],
#     #             font="10px arial grey",
#     #             size=15,
#     #         )
#     #         for u in net.nodes
#     #     ]

#     #     edges = [
#     #         Edge(
#     #             source=v["from"],
#     #             target=v["to"],
#     #             # title=v["title"],
#     #             # color=v["color"],
#     #         )
#     #         for v in net.edges
#     #     ]
#     #     config2 = Config2(height=600,
#     #                     width=1000,
#     #                     nodeHighlightBehavior=True,
#     #                     highlightColor="#F7A7A6",
#     #                     directed=True,
#     #                     collapsible=True,
#     #                     physics=False, 
#     #                     staticGraphWithDragAndDrop=True,
#     #                     link={'labelProperty': 'label', 'renderLabel': True}
#     #                     )

    
#     #     return_value = agraph(nodes, edges, config=config2)
#     #     cols=st.columns(2)
#     #     with(cols[0]):
#     #         st.write("node list\n", H.nodes)
#     #     with (cols[1]):
#     #         st.write("edge list \n",H.edges)
        


