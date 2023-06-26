import gseapy as gp
import numpy as np
import pandas as pd
import plotly.express as px
import psycopg2
import streamlit as st

# st.write("# In Progress...")


@st.cache_resource
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])


conn = init_connection()



def get_infos(df):
    list_cols=df.columns.tolist()
    col_sel = st.radio('Column to select in the df',list_cols)
    df_cpds_in_kegg = pd.DataFrame()
    # df_cpds_tot = get_compound_info()
    if col_sel=='keggid':
        df=df[df.keggid !='No result']
        list_gene=df['keggid'].unique().tolist()
        batch_list = []
        for batch in list_gene:
            batch_list.append("'" + batch + "'")
        sql_kegg ="select * from keggcpdgene inner join gene on gene.geneid=keggcpdgene.geneid where keggid in (" + ",".join(batch_list) + ")"

        df_cpds_in_kegg= pd.read_sql(sql_kegg,conn)
        df_cpds_in_kegg = df_cpds_in_kegg.loc[:,~df_cpds_in_kegg.columns.duplicated()].copy()
        st.write('DF',df_cpds_in_kegg)
        # st.write(list_gene)
    if col_sel=='geneid':
        list_gene=df['geneid'].unique().tolist()
        batch_list = []
        for batch in list_gene:
            batch_list.append("'" + batch + "'")
        sql_kegg ="select * from gene where geneid in (" + ",".join(batch_list) + ")"

        df_cpds_in_kegg= pd.read_sql(sql_kegg,conn)
        df_cpds_in_kegg = df_cpds_in_kegg.loc[:,~df_cpds_in_kegg.columns.duplicated()].copy()
        st.write('DF',df_cpds_in_kegg)





    dftiti=df_cpds_in_kegg.drop_duplicates(keep='first')
    dftiti.dropna(subset='symbol',axis=0,inplace=True)
    # st.write("Data from Kegg:", df_cpds_in_kegg)
    st.write("Data from Kegg:", dftiti)
    human = gp.get_library_name(organism='Human')
    col1,col2 = st.columns(2)
    with col1:
        sel_col = st.selectbox(':red[Choose column] :worried:',dftiti.columns.to_list())

    with col2:
        sel = st.selectbox(':green[Choose DB]',human)


    gene_list=dftiti['symbol'].squeeze().str.strip().to_list()

    enr = gp.enrichr(gene_list=gene_list, # or "./tests/data/gene_list.txt",
                 gene_sets=[sel],
                 organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )
    # st.write(gene_list)
    df_enr = enr.results

    df_enr['Log_10_Pv']=-np.log10(df_enr['Adjusted P-value'])
    df_enr.rename(columns={'Term': 'Pathways'}, inplace=True)
    st.write(df_enr)

    fig = px.bar(df_enr,x='Pathways',y='Log_10_Pv')
    st.plotly_chart(fig, theme="streamlit", use_container_width=True)



def upload_files():
    uploaded_files = st.file_uploader("## Choose files with keggid in cols, if df null", accept_multiple_files=True)
    # st.write(uploaded_file)
    list_df = []
    for uploaded_file in uploaded_files:
        if ".csv" in uploaded_file.name:
            list_df.append(pd.read_csv(uploaded_file))

        if ".fth" in uploaded_file.name:
            list_df.append(pd.read_feather(uploaded_file))
    return(list_df)

if "df_keggGene" not in st.session_state:

    # st.write(uploaded_file)
    list_df = upload_files()

    if len(list_df) > 0:
        df = pd.concat(list_df)
        st.write("Data from file:", df)
        get_infos(df)
else:
    df = st.session_state["df_keggGene"]
    if len(df)>0:
        st.write("Data from Get Structures:", df)
        get_infos(df)
    else:
        st.write('# No data in df: check output of Get Structure, load file(s)')
        list_df = upload_files()

        if len(list_df) > 0:
            df = pd.concat(list_df)
            st.write("Data from file:", df)
            get_infos(df)
