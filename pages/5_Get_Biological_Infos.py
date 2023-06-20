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


col_sel = st.radio('Column to select in the df',('name','keggid'))


def get_infos(df):



    df_cpds_in_kegg = pd.DataFrame()
    df_cpds_tot = get_compound_info()
    if col_sel=='keggid':
        df=df[df.keggid !='No result']
        list_gene=df['keggid'].unique().tolist()
        # st.write(list_gene)

        df_cpds_in_kegg= df_cpds_tot[df_cpds_tot['keggid'].isin(list_gene)].reset_index()
    if col_sel=='name':
        list_name=df['name'].unique().tolist()
        df_cpds_in_kegg= df_cpds_tot[df_cpds_tot['name'].isin(list_name)].reset_index()

    # df_cpds_in_kegg=df_cpds_in_kegg[df_cpds_in_kegg.geneid !='None']
    dftiti=df_cpds_in_kegg.drop_duplicates(keep='first')
    dftiti.dropna(subset='symbol',axis=0,inplace=True)
    # st.write("Data from Kegg:", df_cpds_in_kegg)
    st.write("Data from Kegg:", dftiti)
    human = gp.get_library_name(organism='Human')
    col1,col2 = st.columns(2)
    with col1:
        sel_col = st.selectbox(':red[Choose column] :worried:',df_cpds_tot.columns.to_list())

    with col2:
        sel = st.selectbox(':green[Choose DB]',human)
    # go_mf = gp.get_library(name=sel, organism='Human')

    # st.write("DataTypes of Selection:", go_mf.keys())

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
    # from gseapy.plot import barplot, dotplot
    # fig2 = barplot(enr.res2d,title='GO Biological Processes seroneg day 1 (up)',color = 'r')
    # st.pyplot(fig2)

    # path_names=df_enr['Pathways'].unique().tolist()
    # path_name =st.selectbox('Chose Pathway',path_names)
    #  # ### GET PATHWAYS_ID IN DATA
    # b_list=["'" + b + "'"for b in path_names]
    # sql_path = "select * from pathway where name  in (" + ",".join(b_list) + ")"
    # df_path_names = pd.read_sql(sql_path, conn)
    # df_path_names.drop('servername',inplace=True,axis=1)
    # df_conc=df_path_names.merge(df_enr,left_on='name',right_on='Pathways')
    # st.write("Data from pathway:", df_path_names)

    # string_path = "http://www.kegg.jp/kegg-bin/show_pathway?" +path_name+'/'
    # list_gene=df_cpds_in_kegg['geneid'].tolist()
    # # list_color=df_c['hex_colors'].tolist()
    # for g in list_gene:
    #     string_path += g.lower() + "%09,"
    # string_path += "/default%3dpink/"
    # string_path += "%09,blue"

    # st.write(string_path, unsafe_allow_html=True)

    # ### GET PATHWAYS IN DATA
    # b_list=["'" + b.lower() + "'"for b in df_cpds_in_kegg['geneid']]
    # sql_genepath = "select * from genepath where geneid  in (" + ",".join(b_list) + ")"
    # df_path = pd.read_sql(sql_genepath, conn)
    # df_path= df_path[~df_path.pathid.str.contains('REACT')]

    # st.write("Data from genepath:", df_path)

    # ### GET PATHWAYS_NAMES IN DATA
    # b_list=["'" + b.lower() + "'"for b in df_path['pathid']]
    # sql_path = "select * from pathway where pathid  in (" + ",".join(b_list) + ")"
    # df_path_names = pd.read_sql(sql_path, conn)
    # df_path_names.drop('servername',inplace=True,axis=1)
    # df_conc=df_path.merge(df_path_names,left_on='pathid',right_on='pathid')
    # st.write("Data from pathway:", df_conc)


    # ### GET ALL GENES IN PATHS
    # b_list=["'" + b.lower() + "'"for b in df_path['pathid'].unique()]
    # sql_genepath = "select * from genepath where pathid  in (" + ",".join(b_list) + ")"
    # df_path_tot = pd.read_sql(sql_genepath, conn)
    # df_path_tot= df_path_tot[~df_path_tot.pathid.str.contains('REACT')]



    # ### COMPUTE NUMBER OF GENES IN PATHs
    # df_path_tot['#Genes'] = df_path_tot.groupby(['pathid']).transform('count')

    # df_path_tot.drop('geneid',axis=1,inplace=True)
    # df_path_tot.drop_duplicates(subset='pathid',keep='first',inplace=True)
    # # st.write("# genes in pathways:", df_path_tot)

    # z = df_conc['pathid'].value_counts()
    # df_conc['#Genes'] = df_conc['pathid'].map(z)
    # df_conc.drop('geneid',axis=1,inplace=True)
    # df_conc.drop_duplicates(subset='pathid',keep='first',inplace=True)

    # # df_conc['#Genes'] = df_conc.groupby(['pathid']).transform('count')
    # # st.write(df_conc)

    # df_m = df_conc.merge(df_path_tot,left_on='pathid',right_on='pathid')
    # df_m['Percentage']=df_m['#Genes_x']/df_m['#Genes_y']
    # st.write(df_m)





def get_compound_info():
    df_kegg_d=pd.read_csv("/mnt/shares/L/DB/KEGG/Drugs_And_CPDS_inKEGG/drug_in_all_kegg.csv")
    df_kegg_c=pd.read_csv("/mnt/shares/L/DB/KEGG/Drugs_And_CPDS_inKEGG/drug_in_all_kegg.csv")
    df_kegg=pd.concat([df_kegg_d,df_kegg_c],axis=0)
    return df_kegg



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

if "df_kegg" not in st.session_state:

    # st.write(uploaded_file)
    list_df = upload_files()

    if len(list_df) > 0:
        df = pd.concat(list_df)
        st.write("Data from file:", df)
        get_infos(df)
else:
    df = st.session_state["df_kegg"]
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
