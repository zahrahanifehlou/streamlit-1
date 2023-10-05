
import sys

import pandas as pd
import psycopg2
import streamlit as st

from Bio.KEGG.REST import *
sys.path.append('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/')
from streamlib import sql_df
import urllib
conn_meta = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksi_cpds",
    password="12345",
)
# A bit of helper code to shorten long text
def head(text):
    """ Print the first lines lines of the passed text.
    """
    stre =  (text.split('\n'))
    return stre
def pub(sel_dis):
    keg_dis =  head(kegg_get(sel_dis).read())
    list_ref=[]
    list_auth=[]
    list_title=[]
    list_jour=[]
    for item in keg_dis:
        if 'REFERENCE' in item:
            list_ref.append('https://pubmed.ncbi.nlm.nih.gov/'+item.split('REFERENCE')[1].split(':')[1])
        if 'AUTHORS' in item:
            list_auth.append(item.split('AUTHORS')[1])
        if 'TITLE' in item:
            list_title.append(item.split('TITLE')[1])
        if 'JOURNAL' in item:
            list_jour.append(item.split('JOURNAL')[1])
    df_pub=pd.DataFrame()
    df_pub['REFERENCE']=list_ref
    df_pub['AUTHORS']=list_auth
    df_pub['TITLE']=list_title
    df_pub['JOURNAL']=list_jour
    return df_pub


st.title('Diseases Information')
st.subheader("From Disease to CPDS/Genes", divider='rainbow')
sql_dis = 'select * from disease'

df_dis = sql_df(sql_dis,conn_meta)

# st.write(df_dis)
sel_disease = st.selectbox('Chose the disease:', df_dis["disname"].unique())

code_dis =df_dis[df_dis["disname"]==sel_disease]['disid'].values[0]
st.write("---")
st.write('Kegg Code for selected disease: ',code_dis)
tog_ref = st.sidebar.toggle('Get References',help='retrieve main litterature')
if tog_ref:
    df_pubmed=pub(code_dis)
    st.write('REFERENCES:', df_pubmed)
st.write("---")

col1,col2=st.columns(2)

sql_known_drugs=f"select * from keggcpddis where disid='{code_dis}'"
df_known_drug=sql_df(sql_known_drugs,conn_meta)
# st.write(df_known_drug)


with col1:
    st.write("---")
    st.write('Known Drugs for this disease',df_known_drug)
    st.write("---")


#################################GET Genes ##################

sql_known_genes=f"select * from genedis where disid='{code_dis}'"
df_link_genes=sql_df(sql_known_genes,conn_meta)


with col2:
    st.write("---")
    st.write('Known genes for this disease',df_link_genes)
    st.write("---")

#

on = st.sidebar.toggle('Retrieve Cpds Info',help='Will retrieve all infos')
if on and not df_known_drug.empty:
    # @st.cache_data
    def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
        return df.to_csv().encode('utf-8')
    list_drug = [f"'{t}'" for t in df_known_drug['keggid']]
    # sql_drug="select * from keggcpd where keggid  in (" + ",".join(list_drug) + ")"
    sql_drug = f"SELECT keggcpd.*, keggcpdgene.geneid,gene.geneid, gene.symbol,cpdbatchs.batchid FROM keggcpd\
            LEFT JOIN keggcpdgene ON keggcpdgene.keggid=keggcpd.keggid\
            LEFT JOIN gene ON gene.geneid=keggcpdgene.geneid\
            LEFT JOIN cpdbatchs ON cpdbatchs.pubchemid=keggcpd.pubchemid\
            WHERE keggcpd.keggid IN ({','.join(list_drug)})"
    # sql_drug="select * from keggcpd"sql_df
    df_drug=sql_df(sql_drug,conn_meta)
    df_drug = df_drug.loc[:, ~df_drug.columns.duplicated()].copy()
    df_drug=df_drug.drop_duplicates(subset=["keggid", "batchid"]).reset_index(drop=True)
    st.write("Drug Infos",df_drug)
    on_save1 = st.sidebar.toggle('Save Drug Infos')
    if on_save1:
        st.download_button(label="Save",data=convert_df(df_drug),file_name="data.csv",mime='csv')


    list_drug_ksi = [f"'{t}'" for t in df_drug['batchid']]
    sql_ksi = f"select * from platemap where platemap.batchid in ({','.join(list_drug_ksi)})"
    df_ksi = sql_df(sql_ksi,conn_meta)
    list_source= df_ksi.source.unique()
    sel_source = st.selectbox('Chose the source',list_source )
    df_sel = df_ksi[df_ksi['source']==sel_source]
    st.write(f'Data from source: {sel_source}',df_sel)

    on_save = st.sidebar.toggle('Save Data from source')
    if on_save:
        st.download_button(label="Save",data=convert_df(df_sel),file_name="data.csv",mime='csv')


on_gene = st.sidebar.toggle('Retrieve Genes Info',help='Will retrieve all gene infos')
if on_gene and not df_link_genes.empty:
    list_genes = [f"'{t.upper()}'" for t in df_link_genes['geneid']]
    sql_genes="select * from gene where geneid  in (" + ",".join(list_genes) + ")"
    # sql_drug="select * from keggcpd"
    df_genes=sql_df(sql_genes,conn_meta)
    st.write("Genes Infos",df_genes)

st.subheader("From gene to diseases", divider='rainbow')
sql_help = "select gene.symbol from gene"
df_gene_help=sql_df(sql_help,conn_meta)
list_help_genes=df_gene_help['symbol']
gene_of_int = st.selectbox("chose one of them", list_help_genes)
# gene_of_int= st.text_input('enter your gene of interest').upper()
sql_known_dis=f"select * from gene where symbol='{gene_of_int}'"
df_gene_of_int=sql_df(sql_known_dis,conn_meta)
if not df_gene_of_int.empty:
    geneid=df_gene_of_int[df_gene_of_int["symbol"]==gene_of_int]['geneid'].values[0]
    # st.write(symbol)
    sql_known_gene_in_dis=f"select * from genedis where geneid='{geneid}'"
    df_gene_in_dis=sql_df(sql_known_gene_in_dis,conn_meta)
    if not df_gene_in_dis.empty:
        disid=df_gene_in_dis["disid"].values
        df_list_dis = df_dis[df_dis["disid"].isin(disid)].reset_index(drop=True)
        st.write("known diseases involving this gene",df_list_dis)
        sel_dis = st.selectbox('select the disease for references', df_list_dis['disid'].to_list())
        if tog_ref:
            df_pub = pub(sel_dis)
            st.write('REFERENCES',df_pub)
else:
    st.warning("Unknown gene")


# from tdc.multi_pred import GDA
# data = GDA(name = 'DisGeNET')
# split = data.get_split()

# st.write(split)
