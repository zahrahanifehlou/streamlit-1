

import sys

import pandas as pd
import plotly.express as px
import psycopg2
import streamlit as st
import umap
from sklearn.cluster import KMeans

sys.path.append('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/')
from streamlib import sql_df

st.set_page_config(

    layout="wide",

)

# def init_connection():
#     return psycopg2.connect(**st.secrets["postgres"])
# conn = init_connection()
# cols=st.columns(2)
# with cols[0]:

#     sql = st.text_area("write your SQL")
#     st.write(" example")
#     st.write("SELECT cpdgene.geneid, cpdgene.pubchemid, gene.symbol FROM cpdgene\
#              INNER JOIN gene ON  gene.geneid=cpdgene.geneid WHERE UPPER(gene.symbol) LIKE UPPER( 'H%')")

# with cols[1]:
#      st.image("DB.png")
# if len(sql)>0:
#     st.write(f" you sql code is : {sql}")

#     df_resutls = sql_df(sql, conn).reset_index(drop=True)
#     st.write(df_resutls)
# conn.close()


# def cluster_features(df: pd.DataFrame, fraction: float):
#     """The list of parameters that defines a cluster.

#     Parameters:
#     ===========
#     df: pd.DataFrame
#         The dataframe to select the features from.
#     fraction: float
#         The fraction of feature values that need to point in the same direction
#         in order to be selected for the cluster.

#     Returns: a list of selected feature names.
#     """
#     df_len = len(df)
#     result = []
#     for feat in df.columns:
#         count_plus = int((df[feat] >= 0.0).sum())
#         count_minus = int((df[feat] < 0.0).sum())
#         value = max(count_plus, count_minus) / df_len
#         if value >= fraction:
#             result.append(feat)
#     return result

# warnings.filterwarnings("ignore")

st.write('## Clustering')
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])
conn = init_connection()

conn_profileDB = psycopg2.connect(
        host="192.168.2.131",
        port="5432",
        user="arno",
        database="ksilink_cpds",
        password="12345",
    )
sql_source='select platemap.source from platemap'
df_source=sql_df(sql_source,conn)

sel_source=st.selectbox('Select the source',df_source['source'].unique())
if 'CRISPER' not in sel_source: 
    sql_data = f"select cpd.pubchemid,cpd.keggid,cpd.cpdname,cpdgene.geneid,cpdbatchs.batchid from cpd \
                        inner join cpdbatchs on cpdbatchs. pubchemid=cpd.pubchemid \
                        inner join platemap on platemap.batchid= cpdbatchs.batchid \
                        left join cpdgene on cpdbatchs.pubchemid= cpdgene.pubchemid \
                        where platemap.source='{sel_source}' and cpdgene.server='KEGG' \
                        group by  cpdbatchs.batchid,cpd.pubchemid,cpd.keggid,cpd.cpdname,cpdgene.geneid"
# sql_data="select * from platemap"
else:
    sql_data='select * from crisperbatchs'
data = sql_df(sql_data, conn)
if 'CRISPER' not in sel_source:
    data.dropna(subset='keggid',inplace=True)
    data.drop_duplicates(subset='pubchemid',inplace=True)
    data = data.loc[:,~data.columns.duplicated()].copy()
for col in data.columns:
    data=data.rename(columns={col:'meta_'+col})
# st.write(data.describe())
data=data.reset_index(drop=True)
st.write(data.sample(2))

list_batch = [f"'{meta_batchid}'" for meta_batchid in data["meta_batchid"]]

# st.write(len(data))
sql_drug_profile = f"select * from aggcombatprofile where metabatchid  in (" + ",".join(list_batch) + ")"
# sql_drug_profile = f"SELECT * FROM aggcombatprofile \
# WHERE metabatchid IN ({','.join(list_batch)})"
df_prof_drug = sql_df(sql_drug_profile, conn_profileDB)
df_prof_drug_meta=df_prof_drug.merge(data,left_on='metabatchid',right_on='meta_batchid').reset_index(drop=True)
num_col=[col for col in df_prof_drug_meta if 'meta' not in col]

# st.write(df_prof_drug_meta)

list_geneid = [f"'{geneid}'" for geneid in df_prof_drug_meta["meta_geneid"]]
sql_gene=f"select * from gene where gene.geneid in ({','.join(list_geneid)})"
df_genes = sql_df(sql_gene,conn)
# df_genes.dropna(subset='mainlocation',inplace=True)
for col in df_genes.columns:
    df_genes=df_genes.rename(columns={col:'meta_'+col})


# st.write(df_genes)

df_prof_drug_meta_genes =df_prof_drug_meta.merge(df_genes,left_on='meta_geneid',right_on='meta_geneid').reset_index(drop=True)
df_prof_drug_meta_genes=df_prof_drug_meta_genes[df_prof_drug_meta_genes['metasource']==sel_source].reset_index(drop=True)
st.write(df_prof_drug_meta_genes.sample(2))


model = umap.UMAP(random_state=42, verbose=False).fit(df_prof_drug_meta_genes[num_col])
emb = model.transform(df_prof_drug_meta_genes[num_col])

df_all_umap = pd.DataFrame()
df_all_umap["X_umap"] = emb[:, 0]
df_all_umap["Y_umap"] = emb[:, 1]
df_all_umap["symbol"] = df_prof_drug_meta_genes["meta_symbol"]
# df_all_umap["location"] = df_prof_drug_meta_genes["meta_mainlocation"]

fig3 = px.scatter(
    df_all_umap,
    x="X_umap",
    y="Y_umap",
    title="umap",
    hover_data=["symbol"],
    # color="symbol",
)
st.plotly_chart(fig3, theme="streamlit", use_container_width=True)#

numerics = ["float16", "float32", "float64"]
X=df_all_umap.select_dtypes(include=numerics)
nb_cluster=st.slider('Number of clusters',min_value=2,max_value=30,value=2,step=1)
kmeans = KMeans(n_clusters=nb_cluster, random_state=0, n_init="auto").fit(X)
df_all_umap['Cluster']=kmeans.labels_
list_chr=[]
for g in df_all_umap['Cluster']:
    list_chr.append(chr(g+97))
df_all_umap['Clusterchr']=list_chr
df_all_umap =  df_all_umap.drop('Cluster',axis=1)
fig4 = px.scatter(
    df_all_umap,
    x="X_umap",
    y="Y_umap",
    title="umap",
    hover_data=["symbol"],
    color="Clusterchr",
)
st.plotly_chart(fig4, theme="streamlit", use_container_width=True)#
st.write('Unique Genes: ',len(df_prof_drug_meta_genes['meta_geneid'].unique()))
st.write('Unique Symbols: ',len(df_prof_drug_meta_genes['meta_symbol'].unique()))
st.write('Unique Drugs: ',len(df_prof_drug_meta_genes['metabatchid'].unique()))
st.write(df_all_umap)

from sklearn.neighbors import NearestNeighbors


list_symb=df_all_umap['symbol'].to_list()
sel_symb = st.selectbox('Select a gene:', list_symb)
sel_card = st.slider('card of neighbours', min_value=2,max_value=30,value=5,step=1)
neigh = NearestNeighbors(n_neighbors=sel_card,n_jobs=-1)
neigh.fit(df_all_umap.select_dtypes(include=numerics))
neib =neigh.kneighbors(df_all_umap[df_all_umap['symbol']==sel_symb].select_dtypes(include=numerics))[1].tolist()
df_neib = df_all_umap.iloc[neib[0]]

st.write(df_neib)

df_all_umap['color']='others'
df_all_umap.loc[df_all_umap["symbol"].isin(df_neib['symbol'].tolist()), "color"] = "similar profile"

fig5 = px.scatter(
    df_all_umap,
    x="X_umap",
    y="Y_umap",    
    hover_data=["symbol"],
    color_discrete_sequence=["green", "red","blue" ],
    title=f"similar CRISPER profiles to {sel_symb}   ",
    color="color",
    opacity=0.2
)
for trace in fig5.data:
    if trace.name=='similar profile':
        # st.write(trace)
        trace.marker.opacity=0.9
        trace.marker.size=15
    # else:
    #     trace.marker.line.color = 'rgba(0, 0, 0, 1.0)'
    #     trace.marker.line.width = 2
        
st.plotly_chart(fig5, theme="streamlit", use_container_width=True)#
# st.error('you failed!!!')