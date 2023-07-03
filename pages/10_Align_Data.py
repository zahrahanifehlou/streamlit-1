
import warnings

import bbknn
import numpy as np
import pandas as pd
import plotly.express as px
import psycopg2
import streamlit as st
import umap

conn_profileDB = psycopg2.connect(
        host="192.168.2.131",
        port="5432",
        user="arno",
        database="ksilink_cpds",
        password="12345",
    )

col1, col2 = st.columns(2)
list_sources = [
    "source_1",
    "source_2",
    "source_3",
    "source_4",
    "source_5",
    "source_6",
    "source_7_25",
    "source_7_625",
    "source_8",
    "source_9",
    "source_10",
    "source_11",
    "source_13",
]
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])
conn = init_connection()

#  where batchid like 'BR%'
# choix_source = st.selectbox("Select the Source", list_sources)

with col1:
    option = st.selectbox("Pick one table", list_sources,key='src1')
    sql_src1 = f"SELECT * FROM platemap where project='{option}'"
    df_src1 = pd.read_sql(sql_src1,conn)
    st.write("You selected the source: ", df_src1)



with col2:
    option2 = st.selectbox("Pick one table", list_sources,key='src2')
    sql_src2 = f"SELECT * FROM platemap where project='{option2}'"
    df_src2 = pd.read_sql(sql_src2,conn)
    st.write("You selected the source: ", df_src2)

# df_join=df_src1.merge
res = set(df_src1['batchid']).intersection(df_src2['batchid'])
st.write(len(res))

# df_conc= pd.concat([df_src1,df_src2]).reset_index()
# num_col=[col for col in df_conc if 'meta' not in col]
# model = umap.UMAP(random_state=42, verbose=False).fit(df_conc[num_col])
# emb = model.transform(df_conc[num_col])

# df_all_umap = pd.DataFrame()
# df_all_umap["X_umap"] = emb[:, 0]
# df_all_umap["Y_umap"] = emb[:, 1]
# df_all_umap["source"] = df_conc["metasource"]
# # df_all_umap["location"] = df_prof_drug_meta_genes["meta_mainlocation"]

# fig3 = px.scatter(
#     df_all_umap,
#     x="X_umap",
#     y="Y_umap",
#     title="umap",
#     hover_data=["source"],
#     color="source",
# )
# st.plotly_chart(fig3, theme="streamlit", use_container_width=True)#
