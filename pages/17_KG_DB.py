import streamlit as st
import pandas as pd
from age import *
import streamlit.components.v1 as components
from pyvis.network import Network
import igraph as ig
from streamlit_agraph import agraph, Node, Edge, Config, ConfigBuilder
import matplotlib
from sklearn import preprocessing
import psycopg2
from age.networkx import *
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
list_gene = [t.strip().upper() for t in var_t if t != ""]
if len(list_gene) >0:
    gsql = (
        f"""SELECT * from cypher('%s', $$ MATCH p=(n )<-[r*2]-() 
        where n.__id__ in {list_gene} 
        RETURN  p   $$) as (v agtype)"""
        % graphName
    )
    st.write(gsql)
    G = age_to_networkx(conn, graphName, 
                        query=gsql )
    label_dict=dict(G.nodes(data="properties", default=1))
    converted_dict = {key: value['__id__'] for key, value in label_dict.items()}
    H = nx.relabel_nodes(G, converted_dict)
    H = nx.relabel_nodes(G, converted_dict)
    pos = nx.spring_layout(H)
    nx.draw_networkx_labels(H, pos)
    nx.draw_networkx_edges(H, pos, edge_color='r', arrows = True)
    nx.draw_networkx_nodes(H,pos)
 
    # nt = Network(notebook=False)
    # nt.from_nx(H)
    # nt.show('nx.html')