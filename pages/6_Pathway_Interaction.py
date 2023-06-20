
import xml.etree.ElementTree as ET

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psycopg2
import requests
import streamlit as st
from Bio.KEGG import REST
from bioservices import KEGG
from bs4 import BeautifulSoup

# st.write("# In Progress...")


@st.cache_resource
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])

def convert_hexa(dfa,col):
    cmap = plt.cm.get_cmap('bwr')
    hexcolors=[]
    # normalized_data = (dfa[col] - np.min(dfa[col])) / (np.max(dfa[col]) - np.min(dfa[col]))
    normalized_data=dfa[col]
    # colors = cmap(np.arange(-1.0,1.0,2.0/(len(normalized_data))))
    # colors=normalized_data
    colors=cmap(normalized_data)
    hex_colors = [mcolors.rgb2hex(color) for color in colors]
    hexcolors2=[]
    for cols in hex_colors:
        hexcolors2.append(cols.split('#')[1])
    dfa['hex_colors']=hex_colors
    dfa['num_colors']=normalized_data
    return dfa


conn = init_connection()

uploaded_files = st.file_uploader("Choose files", accept_multiple_files=True)
# st.write(uploaded_file)
list_df = []
for uploaded_file in uploaded_files:
    if ".csv" in uploaded_file.name:
        list_df.append(pd.read_csv(uploaded_file))

    if ".fth" in uploaded_file.name:
        list_df.append(pd.read_feather(uploaded_file))

if len(list_df) > 0:
    df = pd.concat(list_df)
    st.write(df)
    batch_list = df["batchid"].tolist()

    b_list=[]
    for b in batch_list:
        if "JCP2022_80" in b:
            b_list.append("'" + b + "'")
            # b_list.append(b.upper())

    sql_crispr = "select * from crispermeta where batchid  in (" + ",".join(b_list) + ")"
    df_int = pd.read_sql(sql_crispr, conn)

    df_int_int = df.merge(df_int, left_on='batchid',right_on='batchid')

    b_list2 = df_int["geneid"].to_list()
    if len(b_list2)>0:
        bq = []
        for bs in b_list2:
            bq.append("'" + bs + "'")

        sql_batch = "select * from genepath where geneid in (" + ",".join(bq) + ")"
        df_cpd = pd.read_sql(sql_batch, conn)
    df_a = df_int_int.merge(df_cpd,left_on='geneid',right_on='geneid')
    df_b= df_a[~df_a.pathid.str.contains('REACT')]

    b_path=df_b['pathid'].tolist()
    b_pathq=[]
    for bs in b_path:
            b_pathq.append("'" + bs + "'")

    sql_path = "select * from pathway where pathid in (" + ",".join(b_pathq) + ")"
    df_path = pd.read_sql(sql_path, conn)

    df_bb=df_b.merge(df_path,left_on='pathid',right_on='pathid')

    # st.write('df_bb')
    # st.write(df_bb)
    path_names=df_bb['name'].unique().tolist()
    path_name =st.selectbox('Chose Pathway',path_names)

    df_cancer=df_bb[df_bb['name']==path_name]
    st.write('path_name')
    # st.write(df_cancer)

    df_c=convert_hexa(df_cancer,'sim')

    st.write('df: '+path_name)
    st.write(df_c)
    # pathway_id='hsa05200'



    # url = f"http://rest.kegg.jp/get/{pathway_id}/kgml"
    # response = requests.get(url)
    # list_g = []
    # tree = ET.fromstring(response.content)
    # for entry in tree.iter('entry'):
    #     if entry.attrib["type"] == "gene":
    #         graphics=entry.find('graphics')
    #         list_g.append(entry.attrib["name"].split(' ')[0])

    # string_path = "http://www.kegg.jp/kegg-bin/show_pathway?" +'hsa05200'+'/'
    # list_gene=df_c['geneid'].tolist()
    # list_color=df_c['hex_colors'].tolist()
    # for g, c in zip(list_gene[:50], list_color[:50]):
    #     if g in list_g:
    #         string_path += g + "%09,"+"%23"+c +'/'

    # # string_path += "/default%3dpink/"
    # # string_path += "%09,blue"
    # st.write(string_path, unsafe_allow_html=True)
    df_save=pd.DataFrame()
    df_save['geneid']=df_c['geneid']
    df_save['hexcolor']=df_c['hex_colors']
    df_save.to_csv(path_name+'.txt',sep='\t',index=None)

    st.write('https://www.genome.jp/kegg/mapper/color.html', unsafe_allow_html=True)
