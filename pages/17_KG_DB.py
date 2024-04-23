import streamlit as st
import pandas as pd
import networkx as nx
import streamlit.components.v1 as components
from pyvis.network import Network
import igraph as ig
from streamlit_agraph import agraph, Node, Edge, Config, ConfigBuilder
import matplotlib
from sklearn import preprocessing
import psycopg2
from age.networkx import *
from age import *
import networkx as nx
import argparse

st.set_page_config(
    layout="wide",
)
st.title("Knowledge Graph (PrimeKG)")
st.header("17080 diseases, 4M relationship,100K Nodes,20 Sources")
st.subheader(
    "Building a knowledge graph to enable precision medicine (Harvard MS)",
    divider="rainbow",
)

# connect DB
graphName = 'kg1'
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

gsql=f"""SELECT * from cypher('%s', $$ MATCH p=(n )<-[r*2]-() 
    where n.__id__ in {list_gene}
 
    RETURN  p   $$) as (v agtype)"""% graphName
st.write(gsql)  