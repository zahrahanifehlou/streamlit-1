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

pd.set_option("display.max_rows", 100)
pd.set_option("display.max_columns", 4000)
pio.templates.default = "plotly_dark"
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.filterwarnings("ignore")
import psycopg2
import umap

st.write('## Incoming....')

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
list_sources = [
        "src_7_25_combat",
        "src_7_625_combat",
    ]
choix_source = st.selectbox("Select the Source", list_sources)


conn = psycopg2.connect(host="192.168.2.131", port="5432", user="arno", database="ksi_cpds", password="12345")
sql_kegg = "select * from allkegg"

df_drug_meta = pd.read_sql(sql_kegg,conn)
df_drug_meta.dropna(subset='geneid',axis=0,inplace=True)

b_list2 = df_drug_meta["keggid"].to_list()
bq = []
for bs in b_list2:
    bq.append("'" + bs + "'")

sql_batch = "select batchs.batchid,cpd.keggid from cpd  inner join batchs on batchs.pubchemid=cpd.pubchemid where keggid  in (" + ",".join(bq) + ")"
df_ksi_pubchem2 = pd.read_sql(sql_batch, conn)
st.write(df_ksi_pubchem2.describe())

df_ksi_prof=pd.read_feather(f'/mnt/shares/L/DB/profiles/{choix_source}.fth')
# st.write("prof",df_ksi_prof.columns)

df_ksi_prof.drop('Metadata_pubchem_cid',inplace=True,axis=1)

df_inter=df_ksi_prof.merge(df_ksi_pubchem2, left_on='Metadata_JCP2022', right_on='batchid').reset_index(drop=True)
# df_inter.drop('morganfp',inplace=True,axis=1)
df_join_kegg_geneid_prof = df_inter.merge(df_drug_meta, left_on='keggid', right_on='keggid').reset_index(drop=True)

col_filt=[c for c in df_join_kegg_geneid_prof.columns if 'Metadata' not in c]

df_join_kegg_geneid_prof=df_join_kegg_geneid_prof[col_filt]
df_join_kegg_geneid_prof.drop('batchid',inplace=True,axis=1)
# st.write(df_join_kegg_geneid_prof)
df_join_kegg_geneid_prof_filt = df_join_kegg_geneid_prof


df_join_kegg_geneid_prof_filt['target']=df_join_kegg_geneid_prof_filt['keggid']+'_'+df_join_kegg_geneid_prof_filt['geneid']
st.write('Cpds/Drugs with Kegg Annotations',df_join_kegg_geneid_prof_filt)
df_join_kegg_geneid_prof_filt=df_join_kegg_geneid_prof_filt.drop('keggid',axis=1)

df_join_kegg_geneid_prof_filt = df_join_kegg_geneid_prof_filt.groupby(['target','geneid','symbol','name','efficacy']).median().reset_index()
# st.write('Cpds/Drugs with Kegg Annotations',df_join_kegg_geneid_prof_filt)


st.write('## Loading CRispr')

batch_list = df_join_kegg_geneid_prof_filt["geneid"].tolist()

b_list=[]
for b in batch_list:
        b_list.append("'" + b.lower() + "'")
        # b_list.append(b.upper())

sql_crispr = "select * from crispermeta where geneid  in (" + ",".join(b_list) + ")"
df_crispr = pd.read_sql(sql_crispr, conn)

tab1,tab2=st.tabs(["Metadata", "Profiles"])


df_crispr_data  = pd.read_feather('/mnt/shares/L/DB/profiles/src_13_combat.fth')
col_filt2=[c for c in df_crispr_data.columns if 'Metadata' not in c]
common_cols = list(set(col_filt).intersection(col_filt2))
df_crispr_data = df_crispr_data[common_cols]
st.write(df_crispr_data.head())

df_crispr_data_j = df_crispr.merge(df_crispr_data,left_on='batchid',right_on='Metadata_JCP2022').reset_index(drop=True)
df_crispr_data_j_agg=df_crispr_data_j.groupby(['batchid','geneid']).median().reset_index()
df_crispr_data_j_agg['target']=df_crispr_data_j_agg["geneid"]

gene_list = df_crispr_data_j_agg["geneid"].tolist()

b_list=[]
for b in gene_list:
        b_list.append("'" + b.lower() + "'")
        # b_list.append(b.upper())

sql_path = "select * from genepath where geneid  in (" + ",".join(b_list) + ")"
df_path = pd.read_sql(sql_path, conn)
df_path= df_path[~df_path.pathid.str.contains('REACT')].reset_index(drop=True)
tab1,tab2,tab3=st.tabs(["Metadata", "Profiles","Pathways"])

tab1.write(df_crispr)
tab2.write(df_crispr_data_j_agg)
tab3.write(df_path)

numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
df_num_crispr = df_crispr_data_j_agg.select_dtypes(include=numerics)
df_num_cpd = df_join_kegg_geneid_prof_filt.select_dtypes(include=numerics)
list_common_cols = list(set(df_num_crispr.columns).intersection(df_num_cpd.columns))



df_num_cpd["efficacy"] = df_join_kegg_geneid_prof_filt["efficacy"].apply(lambda x: x.split(",")[0])
df_num_cpd["target"] = df_join_kegg_geneid_prof_filt["target"]
df_num_crispr["target"]=df_crispr_data_j_agg["target"]
df_num_cpd["Type"]='CPDS'
df_num_crispr["Type"]='CRISPR'
df_num_crispr["efficacy"]='Unknown'

#st.write('join',df_join_kegg_geneid_prof_filt)
list_common_cols2 = list(set(df_num_crispr.columns).intersection(df_num_cpd.columns))
# st.write(list_common_cols2)
df_cpd_gene  = pd.concat([df_num_crispr[list_common_cols2],df_num_cpd[list_common_cols2]]).reset_index(drop=True)
# st.write('CPD_GENES',df_cpd_gene)
list_df=['cpds only', 'crispr only', 'cpds vs crispr']

sel_data= st.selectbox('Select the data',list_df)

if sel_data=='cpds only':
    df=df_num_cpd
elif sel_data=='crispr only':
    df=df_num_crispr
elif sel_data=='cpds vs crispr':
    df=df_cpd_gene

# df=df.dropna(subset='Dr_Gene')
# st.write('df',df['target'])
cols=[c for c in list_common_cols2 if not c.startswith("Meta")]
# df_num_cpd.to_csv('625.csv')
# df.drop(columns=["Unnamed: 0"],inplace=True)
df=df[cols]
st.write("## UMAP")

reducer = umap.UMAP(densmap=True, random_state=42,verbose=True)
embedding = reducer.fit_transform(df.select_dtypes(include=numerics))
df_emb=pd.DataFrame()
df_emb['x']=embedding[:, 0]
df_emb['y']=embedding[:, 1]
df_emb['target']=df['target']
df_emb['Type']=df['Type']
df_emb['efficacy']=df['efficacy']
fig1 = px.scatter(df_emb,x='x',y='y',hover_data=['target'],color='Type')
st.plotly_chart(fig1, theme="streamlit", use_container_width=True)

st.write('## Dendrogram')


pw_func = partial(pdist, metric='euclidean')
fig = create_dendrogram(df_emb.select_dtypes(include=numerics), orientation='left',distfun=pw_func, linkagefun=lambda x: sch.linkage(x, "ward"),labels=df_emb.target.values)
# fig.update_layout(width=2000, height=1600,template='plotly_dark')
st.plotly_chart(fig, theme="streamlit", use_container_width=True)
