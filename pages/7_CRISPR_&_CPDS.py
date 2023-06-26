import warnings
from functools import partial
from pydoc import describe

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import scipy.cluster.hierarchy as sch
import streamlit as st
from plotly.figure_factory import create_dendrogram
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics.pairwise import cosine_similarity

pd.set_option("display.max_rows", 100)
pd.set_option("display.max_columns", 4000)
pio.templates.default = "plotly_dark"
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.filterwarnings("ignore")
import psycopg2
import umap

st.write("## Based on Kegg Annotations....")


conn = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksi_cpds",
    password="12345",
)
sql_kegg = "select * from keggcpdgene inner join cpd on cpd.keggid=keggcpdgene.keggid \
            inner join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid "

df_drug_meta = pd.read_sql(sql_kegg, conn)
df_drug_meta = df_drug_meta.loc[:,~df_drug_meta.columns.duplicated()].copy()


df_drug_meta.dropna(subset="geneid", axis=0, inplace=True)
# st.write(df_drug_meta.describe())


b_list2 = df_drug_meta["batchid"].to_list()
bq = []
for bs in b_list2:
    bq.append("'" + bs + "'")

conn2 = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksilink_cpds",
    password="12345",
)
sql_profile = (
    "select * from aggcombatprofile where metabatchid  in ("
    + ",".join(bq)
    + ")"
)

df_cpd_prof=pd.read_sql(sql_profile,conn2)
st.write(df_cpd_prof.sample(5))
list_sources = df_cpd_prof['metasource'].unique().tolist()
choix_source = st.selectbox("Select the Source", list_sources)

df_cpd_src=df_cpd_prof[df_cpd_prof['metasource']==choix_source].reset_index(drop=True)
# st.write(df_cpd_src.sample(5))

df_cpds_merge = df_cpd_src.merge(df_drug_meta,left_on='metabatchid',right_on='batchid').reset_index(drop=True)

# st.write('Merged', df_cpds_merge.sample(5))
b_list2 = df_cpds_merge["keggid"].to_list()
bq = []
for bs in b_list2:
    bq.append("'" + bs + "'")

sql_cpd_meta = "select * from keggcpd where keggid  in (" + ",".join(bq) + ")"
df_cpd_meta = pd.read_sql(sql_cpd_meta, conn)
df_cpds_merge_meta = df_cpds_merge.merge(df_cpd_meta,left_on='keggid',right_on='keggid').reset_index(drop=True)
st.write('Merged with Meta', df_cpds_merge_meta)

# df_cpds_merge_meta.to_csv('cpd_merge.csv',index=None)
################################# CRISPR ######################################

st.write("## Loading CRispr")
b_list2 = df_cpds_merge["geneid"].to_list()
bq = []
for bs in b_list2:
    bq.append("'" + bs + "'")

sql_crispr = "select * from crisperbatchs where geneid  in (" + ",".join(bq) + ")"
df_crispr = pd.read_sql(sql_crispr, conn)

b_list2 = df_crispr["batchid"].to_list()
bq = []
for bs in b_list2:
    bq.append("'" + bs + "'")

sql_crispr_prof = "select * from aggcombatprofile where metabatchid  in (" + ",".join(bq) + ")"
df_crispr_prof = pd.read_sql(sql_crispr_prof, conn2)
df_crispr_merge = df_crispr_prof.merge(df_crispr,left_on='metabatchid',right_on='batchid').reset_index(drop=True)

# st.write('CRISPR', df_crispr_merge.sample(5))


################################# COMPUTATION ######################################
st.write('## Computation part')

numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
df_num_crispr = df_crispr_merge.select_dtypes(include=numerics)
df_num_cpd = df_cpds_merge_meta.select_dtypes(include=numerics)
list_common_cols = list(set(df_num_crispr.columns).intersection(df_num_cpd.columns))


df_num_cpd["meta_efficacy"] = df_cpds_merge_meta["efficacy"].apply(
    lambda x: x.split(",")[0]
)
df_num_cpd["meta_target"] = df_cpds_merge_meta["geneid"]+'_'+df_cpds_merge_meta['keggid']
df_num_crispr["meta_target"] = df_crispr_merge["geneid"]
df_num_cpd["meta_Type"] = "CPDS"
df_num_crispr["meta_Type"] = "CRISPR"
df_num_crispr["meta_efficacy"] = "Unknown"


# st.write('join',df_join_kegg_geneid_prof_filt)
list_common_cols2 = list(set(df_num_crispr.columns).intersection(df_num_cpd.columns))
# st.write(list_common_cols2)
df_cpd_gene = pd.concat(
    [df_num_crispr[list_common_cols2], df_num_cpd[list_common_cols2]]
).reset_index(drop=True)
# st.write('CPD_GENES',df_cpd_gene)
list_df = ["cpds only", "crispr only", "cpds vs crispr"]

sel_data = st.selectbox("Select the data", list_df)

if sel_data == "cpds only":
    df = df_num_cpd
elif sel_data == "crispr only":
    df = df_num_crispr
elif sel_data == "cpds vs crispr":
    df = df_cpd_gene

st.write('df',len(df),len(df.columns))

# df=df.dropna(subset='Dr_Gene')
# st.write('df',df['target'])
cols = [c for c in list_common_cols2 if not c.startswith("meta")]
# df_num_cpd.to_csv('625.csv')
# df.drop(columns=["Unnamed: 0"],inplace=True)
# df = df[cols]
st.write("## UMAP" )

reducer = umap.UMAP(densmap=True, random_state=42, verbose=True)
embedding = reducer.fit_transform(df.select_dtypes(include=numerics))
df_emb = pd.DataFrame()
df_emb["x"] = embedding[:, 0]
df_emb["y"] = embedding[:, 1]
df_emb["target"] = df["meta_target"]
df_emb["Type"] = df["meta_Type"]
df_emb["efficacy"] = df["meta_efficacy"]
fig1 = px.scatter(df_emb, x="x", y="y", hover_data=["target"], color="efficacy")
st.plotly_chart(fig1, theme="streamlit", use_container_width=True)



st.write("## Similarity")
simi = cosine_similarity(df[cols])

sim_df = pd.DataFrame(simi)
sim_df['target']=df['meta_target']
st.write('sim_df',sim_df)
fig = px.imshow(sim_df, title="Profiles Similarities")
st.plotly_chart(fig)

# pw_func = partial(pdist, metric="euclidean")
# fig = create_dendrogram(
#     df_emb.select_dtypes(include=numerics),
#     orientation="left",
#     distfun=pw_func,
#     linkagefun=lambda x: sch.linkage(x, "ward"),
#     labels=df_emb.target.values,
# )
# # fig.update_layout(width=2000, height=1600,template='plotly_dark')
# st.plotly_chart(fig, theme="streamlit", use_container_width=True)

