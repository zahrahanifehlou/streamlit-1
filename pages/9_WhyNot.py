
import warnings

import numpy as np
import pandas as pd
import plotly.express as px
import psycopg2
import streamlit as st
import umap


def cluster_features(df: pd.DataFrame, fraction: float):
    """The list of parameters that defines a cluster.

    Parameters:
    ===========
    df: pd.DataFrame
        The dataframe to select the features from.
    fraction: float
        The fraction of feature values that need to point in the same direction
        in order to be selected for the cluster.

    Returns: a list of selected feature names.
    """
    df_len = len(df)
    result = []
    for feat in df.columns:
        count_plus = int((df[feat] >= 0.0).sum())
        count_minus = int((df[feat] < 0.0).sum())
        value = max(count_plus, count_minus) / df_len
        if value >= fraction:
            result.append(feat)
    return result

warnings.filterwarnings("ignore")

st.write('## why not?')
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


data = pd.read_sql("select cpd.pubchemid,cpd.keggid,cpd.name,cpdgene.geneid,cpdbatchs.batchid from cpd \
                    inner join cpdbatchs on cpdbatchs. pubchemid=cpd.pubchemid \
                    inner join platemap on platemap.batchid= cpdbatchs.batchid \
                    left join cpdgene on cpdbatchs.pubchemid= cpdgene.pubchemid \
                    where platemap.project='source_7_625' and cpdgene.server='KEGG' \
                    group by  cpdbatchs.batchid,cpd.pubchemid,cpd.keggid,cpd.name,cpdgene.geneid", conn)

data.dropna(subset='keggid',inplace=True)
data.drop_duplicates(subset='pubchemid',inplace=True)
data = data.loc[:,~data.columns.duplicated()].copy()
for col in data.columns:
    data=data.rename(columns={col:'meta_'+col})
# st.write(data.describe())


list_batch = [f"'{meta_batchid}'" for meta_batchid in data["meta_batchid"]]

# st.write(len(data))

sql_drug_profile = f"SELECT * FROM aggcombatprofile \
WHERE aggcombatprofile.metabatchid IN ({','.join(list_batch)}) AND aggcombatprofile.metasource='source_7_625'"
df_prof_drug = pd.read_sql(sql_drug_profile, conn_profileDB)
df_prof_drug_meta=df_prof_drug.merge(data,left_on='metabatchid',right_on='meta_batchid').reset_index(drop=True)
num_col=[col for col in df_prof_drug_meta if 'meta' not in col]

# st.write(df_prof_drug_meta)

list_geneid = [f"'{geneid}'" for geneid in df_prof_drug_meta["meta_geneid"]]
sql_gene=f"select * from gene where gene.geneid in ({','.join(list_geneid)})"
df_genes = pd.read_sql(sql_gene,conn)
# df_genes.dropna(subset='mainlocation',inplace=True)
for col in df_genes.columns:
    df_genes=df_genes.rename(columns={col:'meta_'+col})


# st.write(df_genes)

df_prof_drug_meta_genes =df_prof_drug_meta.merge(df_genes,left_on='meta_geneid',right_on='meta_geneid').reset_index(drop=True)
# st.write(df_prof_drug_meta_genes)


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
    color="symbol",
)
st.plotly_chart(fig3, theme="streamlit", use_container_width=True)#

st.write('Unique Genes: ',len(df_prof_drug_meta_genes['meta_geneid'].unique()))
st.write('Unique Symbols: ',len(df_prof_drug_meta_genes['meta_symbol'].unique()))
st.write('Unique Drugs: ',len(df_prof_drug_meta_genes['metabatchid'].unique()))

st.error('you failed!!!')