
import warnings

import numpy as np
import pandas as pd
import psycopg2
import streamlit as st


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

# data = pd.read_sql("select * from cpd  inner join cpdbatchs \
#                    on cpdbatchs.pubchemid=cpd.pubchemid \
#                    inner join platemap on platemap.batchid=cpdbatchs.batchid \
#                    where platemap.project='source_7_625'", conn)
data = pd.read_sql("select cpd.pubchemid,cpd.keggid,cpd.name,cpdgene.geneid,cpdbatchs.batchid from cpd \
                    inner join cpdbatchs on cpdbatchs. pubchemid=cpd.pubchemid \
                    inner join platemap on platemap.batchid= cpdbatchs.batchid \
                    left join cpdgene on cpdbatchs.pubchemid= cpdgene.pubchemid \
                    where platemap.project='source_7_625' and cpdgene.server='KEGG' \
                    group by  cpdbatchs.batchid,cpd.pubchemid,cpd.keggid,cpd.name,cpdgene.geneid", conn)

data.dropna(subset='keggid',inplace=True)
data.drop_duplicates(subset='pubchemid',inplace=True)
data = data.loc[:,~data.columns.duplicated()].copy()
st.write(data.describe())

# list_keggid = [f"'{keggid}'" for keggid in data["keggid"]]
# sql = f"SELECT * from  keggcpdgene WHERE keggcpdgene.keggid \
#                 IN ({','.join(list_keggid)})"
# df_cpdGene = pd.read_sql(sql, conn).reset_index(drop=True)



# df_join = df_cpdGene.merge(data,left_on='keggid',right_on='keggid')


# st.write(len(df_join['keggid'].unique()),len(df_join['geneid'].unique()))
# df_join.drop_duplicates(subset='pubchemid',inplace=True)
list_batch = [f"'{batchid}'" for batchid in data["batchid"]]

st.write(len(data))

sql_drug_profile = f"SELECT * FROM aggcombatprofile \
WHERE aggcombatprofile.metabatchid IN ({','.join(list_batch)}) AND aggcombatprofile.metasource='source_7_625'"
df_prof_drug = pd.read_sql(sql_drug_profile, conn_profileDB)

st.write(df_prof_drug.describe())
st.error('you failed!!!')