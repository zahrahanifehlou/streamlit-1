
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
    "Amgen" ,
    "Astrazeneca",
     "Bayer",
     "Broad",
     "Eisai",
    "Janssen",
     "Ksilink_25",
     "Ksilink_625",
     "Merck",
    "Pfizer",
     "Servier" ,
    "Takeda" ,
    "CRISPER"
]
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])
conn = init_connection()

#  where batchid like 'BR%'
# choix_source = st.selectbox("Select the Source", list_sources)

def clean_pubchemid(pubchemid):
    return str(pubchemid).split('.')[0]

def sel_code(option, conn):
    #sql = f"SELECT cpdbatchs.pubchemid, cpdbatchs.batchid, platemap.plate, platemap.ctrltype, platemap.well FROM platemap INNER JOIN cpdbatchs ON cpdbatchs.batchid = platemap.batchid WHERE platemap.source = '{option}' "
    sql = f"SELECT cpdbatchs.pubchemid, cpdbatchs.batchid, platemap.plate, platemap.ctrltype, platemap.well FROM platemap \
    INNER JOIN cpdbatchs ON cpdbatchs.batchid = platemap.batchid and platemap.ctrltype= 'poscon_diverse' "
    df = pd.read_sql(sql, conn)
    df['pubchemid'] = df['pubchemid'].apply(clean_pubchemid)
    df=df.drop_duplicates(subset=[ "plate", "well"]).reset_index(drop=True)
    return df

# with col1:
#     option1 = st.selectbox("Pick one table", list_sources, key='src1')
#     df_src1 = sel_code(option1, conn)
#     st.write("You selected the source:", df_src1)

# with col2:
#     option2 = st.selectbox("Pick one table", list_sources, key='src2')
#     df_src2 = sel_code(option2, conn)
#     st.write("You selected the source:", df_src2)

# res = set(df_src1['pubchemid']).intersection(df_src2['pubchemid'])
# st.write(len(res))
# st.write(df_src1[df_src1["ctrltype"]=="poscon_cp"])
# st.write(df_src2[df_src2["ctrltype"]=="poscon_cp"])

# df=pd.concat([df_src2,df_src1])
df = sel_code('titi',conn)
# df=df[df["ctrltype"]=="poscon_cp"]


batch_list = [f"'{batchid}'" for batchid in df["batchid"].unique()]
# src_list = [f"'{s}'" for s in [option2, option1]]
#sql_profile = f"SELECT * FROM aggcombatprofile WHERE metabatchid IN ({','.join(batch_list)}) and metasource IN ({','.join(src_list)})"
sql_profile = f"SELECT * FROM aggcombatprofile WHERE metabatchid IN ({','.join(batch_list)})"
df_prof = pd.read_sql(sql_profile, conn_profileDB)
df_prof=df_prof.drop_duplicates(subset=[ "metabatchid","metasource"]).reset_index(drop=True)
df_prof.reset_index(drop=True, inplace=True)
st.write("profile")
st.write(df_prof)



#df_conc= pd.concat([df_src1,df_src2]).reset_index()
num_col=[col for col in df_prof if 'meta' not in col]
model = umap.UMAP(random_state=42, verbose=False).fit(df_prof[num_col])
emb = model.transform(df_prof[num_col])

df_all_umap = pd.DataFrame()
df_all_umap["X_umap"] = emb[:, 0]
df_all_umap["Y_umap"] = emb[:, 1]
df_all_umap["source"] = df_prof["metasource"]
# df_all_umap["location"] = df_prof_drug_meta_genes["meta_mainlocation"]
st.write(df_all_umap)
st.write(df_all_umap["source"].value_counts())
fig3 = px.scatter(
    df_all_umap,
    x="X_umap",
    y="Y_umap",
    title="umap",
    hover_data=["source"],
    color="source",
)
st.plotly_chart(fig3, theme="streamlit", use_container_width=True)#
conn_profileDB.close()