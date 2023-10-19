
# import xml.etree.ElementTree as ET

import sys

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
# import numpy as np
import pandas as pd
import psycopg2
# import requests
import streamlit as st

sys.path.append('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/')


from Bio import SeqIO
from Bio.Graphics.ColorSpiral import ColorSpiral
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG.REST import *
from streamlib import sql_df


def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])

conn_meta = init_connection()

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



# conn = init_connection()

uploaded_files = st.file_uploader("Choose files with geneid column", accept_multiple_files=True)
# st.write(uploaded_file)
list_df = []
for uploaded_file in uploaded_files:
    if ".csv" in uploaded_file.name:
        list_df.append(pd.read_csv(uploaded_file))

    if ".fth" in uploaded_file.name:
        list_df.append(pd.read_feather(uploaded_file))

if len(list_df) > 0:
    # df = pd.concat(list_df)
    # st.write(df)
    # batch_list = df["batchid"].tolist()

    # b_list=[]
    # for b in batch_list:
    #     if "JCP2022_80" in b:
    #         b_list.append("'" + b + "'")
    #         # b_list.append(b.upper())

    # sql_crispr = "select * from crispermeta where batchid  in (" + ",".join(b_list) + ")"
    # df_int = sql_df(sql_crispr, conn)

    # df_int_int = df.merge(df_int, left_on='batchid',right_on='batchid')
    df = pd.concat(list_df)
    if 'geneid' in df.columns:
        b_list2 = df["geneid"].to_list()
        if len(b_list2)>0:
            # bq = []
            # for bs in b_list2:
            #     bq.append("'" + bs + "'")
            list_geneid = [f"'{t}'" for t in df["geneid"]]

            sql_batch = f"select genepath.pathid,gene.* from genepath\
                        inner join gene on gene.geneid=genepath.geneid\
                        where gene.geneid in  ({','.join(list_geneid)})"
            # st.write(sql_batch)
            df_cpd = sql_df(sql_batch, conn_meta)
            df_a = df.merge(df_cpd,left_on='geneid',right_on='geneid')
            df_b= df_a[~df_a.pathid.str.contains('REACT')]
            df_group=df_b.groupby('pathid').count().reset_index()
            st.write(df_b)
            b_path=df_b['pathid'].tolist()
            b_pathq=[]
            for bs in b_path:
                    b_pathq.append("'" + bs + "'")

            sql_path = "select * from pathway where pathid in (" + ",".join(b_pathq) + ")"
            df_path = sql_df(sql_path, conn_meta)

            df_bb=df_b.merge(df_path,left_on='pathid',right_on='pathid')
            # st.write(df_bb)

            # st.write('df_bb')
            # st.write(df_bb)
            path_names=df_bb['pathname'].unique().tolist()
            path_name =st.selectbox('Chose Pathway',path_names)

            df_cancer=df_bb[df_bb['pathname']==path_name]
            df_cancer['sim']=1
            # st.write(path_name)
            # st.write(df_cancer)

            df_c=convert_hexa(df_cancer,'sim')

            # st.write('df: '+str(path_name))
            # st.write(df_c)
            # # pathway_id='hsa05200'
            # st.write(df_c['pathid'].values)
            img = kegg_get(df_c['pathid'].values[0], "image").read()
            # pathway = KGML_parser.read(kegg_get(df_c['pathid'].values[0], "kgml"))
            # canvas = KGMLCanvas(pathway)
            # canvas.import_imagemap = True
            # file = "fab_map_with_image.pdf"
            # canvas.draw(file)

            st.image(img)

    else:
        st.warning("NO GENEID in your file!!!, please use Convert Data to get it")
