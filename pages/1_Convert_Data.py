
import pandas as pd
import streamlit as st
import psycopg2


st.set_page_config(layout="wide")
import sys

sys.path.append('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/')
from streamlib import sql_df

def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])

conn_meta = init_connection()
# conn_meta.close()

def convert_df(df):
       return df.to_csv(index=False).encode('utf-8')

# def init_connection():
#     return psycopg2.connect(**st.secrets["postgres"])
# conn = init_connection()
sql_query = 'select * from gene'
df_meta = sql_df(sql_query, conn_meta)
st.write("Example in Genes", df_meta.sample(5))
len_df = len(df_meta['geneid'].unique())

st.write(f'{len_df} Unique Genes Registered')

# sql_cpd = "select cpd.*, cpdbatchs.batchid, keggcpdgne.geneid from cpd \
#           INNER join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid\
#           INNER join keggcpdgene on cpd.keggid=keggcpdgene.keggid"
        
sql_cpd="select cpd.*, cpdbatchs.batchid, keggcpdgene.geneid, gene.* from cpd \
        INNER join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid \
        INNER join keggcpdgene on cpd.keggid=keggcpdgene.keggid \
        INNER join gene on gene.geneid=keggcpdgene.geneid"

# st.write(sql_cpd)         
df_cpd_meta = sql_df(sql_cpd, conn_meta)
df_cpd_meta = df_cpd_meta.loc[:, ~df_cpd_meta.columns.duplicated()]
# df_cpd_meta.drop_duplicates()
st.write("Example in Cpds with geneid", df_cpd_meta.sample(5))
len_df = len(df_cpd_meta['pubchemid'].unique())
st.write(f'{len_df} Unique Cpds in Jump with geneid Registered')

sql_kegg="select keggcpdgene.*, gene.*, keggcpd.* from keggcpdgene \
  INNER join gene on keggcpdgene.geneid=gene.geneid \
  INNER join keggcpd on keggcpdgene.keggid=keggcpd.keggid"
# sql_kegg = "select * from keggcpd"
df_kegg = sql_df(sql_kegg,conn_meta)
df_kegg = df_kegg.loc[:, ~df_kegg.columns.duplicated()]
st.write("Kegg Data", len(df_kegg['keggid'].unique()))

####
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
    st.write('Your Data:', st.data_editor(df))
   
#test
# test test
    col1,col2,col3=st.columns(3)
    sel_col = col1.selectbox("Choose your column", df.columns)
    sel_data = col2.selectbox("Choose dataset", ['genes', 'cpds in Jump with geneinfos', 'cpds in Kegg', 'cpds in Jump'])

    b_list2 = df[sel_col].astype(str).to_list()
   
 
    bq = []
    for bs in b_list2:
        bq.append("'" + bs.strip().upper() + "'")

    if sel_data=='genes':
      sel_match = col3.selectbox("Choose matching col", df_meta.columns)
      sql_genes = f"select * from gene where {sel_match}  in  (" + ",".join(bq) + ")"
      df_genes = sql_df(sql_genes,conn_meta)
      st.write(df_genes)
      df_pub=df_genes

      
      st.download_button(
                label="Save",data=convert_df(df_genes),file_name=f"{sel_data}.csv",mime='csv')
    if sel_data=='cpds in Jump with geneinfos':
      sel_match = col3.selectbox("Choose matching col", df_cpd_meta.columns)
      df_cpds= df_cpd_meta[df_cpd_meta[sel_match].astype(str).str.contains('|'.join(b_list2))].drop_duplicates(subset='pubchemid').reset_index(drop=True)
      st.write(df_cpds)
      st.download_button(
                label="Save",data=convert_df(df_cpds),file_name=f"{sel_data}.csv",mime='csv')
      df_pub=df_cpds

    
    if sel_data=='cpds in Kegg':
      sel_match = col3.selectbox("Choose matching col", df_kegg.columns)
      df_kegg= df_kegg[df_kegg[sel_match].astype(str).str.contains('|'.join(b_list2))].drop_duplicates(subset=['keggid']).reset_index(drop=True)
      st.write(df_kegg)
      st.download_button(
                label="Save",data=convert_df(df_kegg),file_name=f"{sel_data}.csv",mime='csv')
      df_pub=df_kegg
    if sel_data=='cpds in Jump':
      sql_jump="select * from cpd"
      df_jump=sql_df(sql_jump,conn_meta)
      sel_match = col3.selectbox("Choose matching col", df_jump.columns)
      df_jump= df_jump[df_jump[sel_match].astype(str).str.contains('|'.join(b_list2))].drop_duplicates(subset=['pubchemid']).reset_index(drop=True)
      st.write(df_jump)
      st.download_button(
                label="Save",data=convert_df(df_jump),file_name=f"{sel_data}.csv",mime='csv')
      df_pub=df_jump
    st.session_state['Convert']=df_pub
conn_meta.close()