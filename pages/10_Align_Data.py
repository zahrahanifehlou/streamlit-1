
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
    # sql = f"SELECT cpdbatchs.pubchemid, cpdbatchs.batchid, platemap.plate, platemap.ctrltype, platemap.well \
    # FROM platemap INNER JOIN cpdbatchs ON cpdbatchs.batchid = platemap.batchid " # WHERE platemap.source = '{option}' "
    sql = f"SELECT cpdbatchs.pubchemid, cpdbatchs.batchid, platemap.plate, platemap.ctrltype, platemap.well FROM platemap \
    INNER JOIN cpdbatchs ON cpdbatchs.batchid = platemap.batchid and platemap.ctrltype= '{option}' "
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
list_cont=['poscon_diverse','poscon_cp','trt','poscon_orf']
option1 = st.selectbox("Pick one ctrl", list_cont, key='src1')
df = sel_code(option1,conn)
# df=df[df["ctrltype"]=="poscon_cp"]
# st.write(df)
# exit(0)


batch_list = [f"'{batchid}'" for batchid in df["batchid"].unique()]
# src_list = [f"'{s}'" for s in [option2, option1]]
#sql_profile = f"SELECT * FROM aggcombatprofile WHERE metabatchid IN ({','.join(batch_list)}) and metasource IN ({','.join(src_list)})"
sql_profile = f"SELECT * FROM aggcombatprofile WHERE metabatchid IN ({','.join(batch_list)})"
df_prof = pd.read_sql(sql_profile, conn_profileDB)
df_prof=df_prof.drop_duplicates(subset=[ "metabatchid","metasource"]).reset_index(drop=True)
st.write(df_prof["metasource"].value_counts())
if option1=='trt':


    samp = st.number_input('sample per source :snail:',min_value=100,max_value=1000,value=200,step=100)
    if samp>200:
        df_prof=df_prof[df_prof["metasource"]!='Broad']
    df_prof=df_prof.groupby('metasource').sample(n=samp).reset_index()
df_prof.reset_index(drop=True, inplace=True)
# st.write("profile")




#df_conc= pd.concat([df_src1,df_src2]).reset_index()
num_col=[col for col in df_prof if 'meta' not in col]
model = umap.UMAP(random_state=42, verbose=False).fit(df_prof[num_col])
emb = model.transform(df_prof[num_col])

df_all_umap = pd.DataFrame()
df_all_umap["X_umap"] = emb[:, 0]
df_all_umap["Y_umap"] = emb[:, 1]
df_all_umap["source"] = df_prof["metasource"]
# df_all_umap["location"] = df_prof_drug_meta_genes["meta_mainlocation"]
# st.write(df_all_umap)
# st.write(df_all_umap["source"].value_counts())
fig3 = px.scatter(
    df_all_umap,
    x="X_umap",
    y="Y_umap",
    title="umap",
    hover_data=["source"],
    color="source",
)
st.plotly_chart(fig3, theme="streamlit", use_container_width=True)#


from sklearn.manifold import TSNE

emb = TSNE(n_components=2, learning_rate='auto',init='random', perplexity=10).fit_transform(df_prof[num_col])
df_tsne = pd.DataFrame()
df_tsne["X_Tsne"] = emb[:, 0]
df_tsne["Y_Tsne"] = emb[:, 1]
df_tsne["source"] = df_prof["metasource"]

fig4 = px.scatter(
    df_tsne,
    x="X_Tsne",
    y="Y_Tsne",
    title="T-SNE",
    hover_data=["source"],
    color="source",
)
st.plotly_chart(fig4, theme="streamlit", use_container_width=True)#

import pacmap

embedding = pacmap.PaCMAP(n_components=2, n_neighbors=None, MN_ratio=0.5, FP_ratio=2.0)
X_transformed = embedding.fit_transform(df_prof[num_col], init="pca")

df_pacmap = pd.DataFrame()
df_pacmap["X_PacMap"] = X_transformed[:, 0]
df_pacmap["Y_PacMap"] = X_transformed[:, 1]
df_pacmap["source"] = df_prof["metasource"]

fig5 = px.scatter(
    df_pacmap,
    x="X_PacMap",
    y="Y_PacMap",
    title="PacMap",
    hover_data=["source"],
    color="source",
)
st.plotly_chart(fig5, theme="streamlit", use_container_width=True)#

conn_profileDB.close()