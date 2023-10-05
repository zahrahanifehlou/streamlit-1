

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
import umap

st.set_page_config(layout="wide")
import sys

sys.path.append('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/')
from streamlib import (conn_meta, conn_prof, get_cpds, get_list_category,
                       get_stringDB, get_stringDB_enr, int_to_str, sql_df,
                       str_to_float)


def get_relation(df_genes):
    if not df_genes.empty:
        thres=0
        list_edges=get_stringDB(df_genes,0,'symbol')
        st.write(f'retrieved : {len(list_edges)} Interaction with {thres} threshold')

        st.write('## Computing Network')

        import networkx as nx
        H = nx.Graph(list_edges)
        import igraph as ig
        import leidenalg as la
        G = ig.Graph.from_networkx(H)
        partition = la.find_partition(G, la.ModularityVertexPartition)
        subg = partition.subgraphs()
        list_gene=[]
        list_clust=[]
        cluster=96
        # thres_db = col_b.slider("cluster Thresholds", 2, 100,9,1)
        st.write(f'Total Number of clusters before Threshold: {len(subg)}')

        ################################### ENRICHMENT ##############################################
        list_cat=get_list_category(df_genes,'symbol')
        categ = st.selectbox("Select Category", list_cat)
        df_go_ento = get_stringDB_enr(df_genes,'symbol',categ)
        st.write(f'{categ} Enrichment',df_go_ento)
        df_go_ento['log_p_val']=np.log10(df_go_ento['p_val'].apply(str_to_float))
        fig_bar=px.bar(df_go_ento,x='Description',y='log_p_val')
        st.plotly_chart(fig_bar)

        import matplotlib.pyplot as plt

        # subax1 = plt.subplot(121)
        st.write('### Displaying the full network/graph clustered with Leiden approach')
        # col1, col2,col3 = st.columns(3)

        vert_size = st.sidebar.slider('vertex size',min_value=0.2,max_value=20.0,step=0.1,value=1.0)
        lab_size = st.sidebar.slider('label size',min_value=0.2,max_value=20.0,step=0.1,value=1.0)
        box_size = st.sidebar.slider('box size',min_value=400,max_value=1600,step=50,value=600)
        ig.config['plotting.backend'] = 'matplotlib'
        subax1=ig.plot(partition,vertex_label=partition.graph.vs['_nx_name'],vertex_label_size=lab_size,vertex_size=vert_size,
                    margin=10, bbox=(0,0,box_size, box_size))
        # nx.draw(GG, with_labels=True,node_size=10,font_size=4)
        st.pyplot(subax1.figure,use_container_width=False)



################################### LOAD GENE CSV DATA ##############################################
for key in st.session_state.keys():
        del st.session_state[key]
df_genes=get_cpds()
unique_gene=df_genes['symbol'].unique()
gene_col='symbol'
# st.write("TOTOT: ",len(df_genes['keggid']))
on = st.sidebar.toggle('Activate Data Upload')
if on:

    uploaded_file = st.file_uploader("Choose csv file with a list of gene symbols  otherwise the reference 300  GeneMOA will be loaded", accept_multiple_files=False)
    if uploaded_file:
        df_genes = pd.read_csv(uploaded_file)
        st.write('Loaded csv file:', df_genes.head())
        gene_col = st.text_input('Write the column name with gene symbols',value='symbol')
        if gene_col:
            df_genes.dropna(subset=gene_col,inplace=True)
            unique_gene=df_genes[gene_col].unique()
            st.write(f'you entered : {len(unique_gene)} genes')
            get_relation(df_genes)
        # st.write(df_genes)
    else:
        df_genes = pd.read_csv('Supp Data 3 - 300-Gene MoA sgRNA Library.csv')
        st.write('Loaded csv file:', df_genes.head())
        gene_col = st.text_input('Write the column name with gene symbols',value='symbol')
    if gene_col:
        unique_gene=df_genes[gene_col].unique()
        
st.write(f'### You entered : {len(unique_gene)} genes' )
disp= st.sidebar.toggle('Display Data')
if disp:
    st.write(df_genes)

################################### LOAD CPDS ##############################################
choice =st.sidebar.radio('Chose between cpds or crispr data',options=['Cpds','CRISPR'])
df_inter=pd.DataFrame()
if choice=='Cpds':
    st.write("## Loading CPDS")

    b_list2 = unique_gene
    bq = []
    for bs in b_list2:
        if bs!='':
            bq.append("'" + bs + "'")


    sql_kegg = "select cpdpath.pathid,keggcpdgene.geneid,gene.symbol, keggcpd.*, cpdbatchs.* from keggcpd\
                inner join keggcpdgene on keggcpd.keggid=keggcpdgene.keggid\
                inner join cpd on cpd.keggid=keggcpd.keggid \
                inner join cpdpath on cpdpath.pubchemid=cpd.pubchemid \
                inner join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid \
                inner join gene on keggcpdgene.geneid=gene.geneid"

    df_drug_meta = sql_df(sql_kegg, conn_meta)
    df_drug_meta = df_drug_meta.loc[:, ~df_drug_meta.columns.duplicated()].copy()
    df_drug_meta = df_drug_meta.drop_duplicates(subset=["keggid", "source"]).reset_index(drop=True)

    # df_drug_meta.dropna(subset="geneid", axis=0, inplace=True)
    # disp2= st.toggle('Display Data')
    if disp:
        st.write(df_drug_meta.sample(5))
    # # st.write(df_drug_meta.describe())

    list_sources = df_drug_meta["source"].unique().tolist()
    choix_source = st.selectbox("Select the Source", list_sources)

    df_cpd_src = df_drug_meta[df_drug_meta["source"] == choix_source].reset_index(drop=True)
    b_list2 = df_cpd_src["batchid"].to_list()
    bq = []
    for bs in b_list2:
        bq.append("'" + bs + "'")


    sql_profile = f"select * from aggcombatprofile where metabatchid  in (" + ",".join(bq) + ")"

    df_cpd_prof = sql_df(sql_profile, conn_prof)
    df_cpd_prof = df_cpd_prof[df_cpd_prof["metasource"] == choix_source].reset_index(drop=True)
    if disp:
        st.write(df_cpd_prof.sample(5))

    # # st.write(df_cpd_src.sample(5))


    df_cpds_merge = df_cpd_prof.merge(df_cpd_src, left_on="metabatchid", right_on="batchid").reset_index(drop=True)


    # st.write("TOTOT: ",unique_gene)
    df_inter=df_cpds_merge[df_cpds_merge[gene_col].isin(unique_gene)].reset_index(drop=True)
    if disp:
        st.write("Merged", df_cpds_merge)
        st.write("CDPS_GENES_COMMON", df_inter)
        st.write(" Unique GENES", df_inter['geneid'].unique())


################################### LOAD CRISPR ##############################################
if choice=='CRISPR':
    st.write("## Loading Crispr")
    b_list2 = unique_gene
    bq = []
    for bs in b_list2:
        bq.append("'" + bs + "'")

    sql_genes = "select synonyms,geneid,chromosome,locus,mainlocation,symbol from gene where symbol  in (" + ",".join(bq) + ")"
    df_genes = sql_df(sql_genes, conn_meta)
    df_genes['chromosome']=df_genes['chromosome'].apply(int_to_str)
    st.write("GeneInfos", df_genes.sample(5))

    b_list2 = df_genes["geneid"].to_list()
    bq = []
    for bs in b_list2:
        bq.append("'" + bs + "'")

    sql_crispr = "select * from crisperbatchs where geneid  in (" + ",".join(bq) + ")"
    df_crispr = sql_df(sql_crispr, conn_meta)

    b_list2 = df_crispr["batchid"].to_list()
    bq = []
    for bs in b_list2:
        bq.append("'" + bs + "'")

    sql_crispr_prof = "select * from aggcombatprofile where metabatchid  in (" + ",".join(bq) + ")"
    df_crispr_prof = sql_df(sql_crispr_prof, conn_prof)
    df_crispr_merge = df_crispr_prof.merge(df_crispr, left_on="metabatchid", right_on="batchid").reset_index(drop=True)

    df_inter=df_crispr_merge.merge(df_genes, left_on="geneid", right_on="geneid").reset_index(drop=True)
    if disp:
        st.write('Crispr with profiles',df_inter)


################################### NETWORK ##############################################
st.write("## Loading StringDB PPI")
# col_a,col_b=st.columns(2)
thres = st.sidebar.slider("Interaction Thresholds", 0.0, 1.0,0.4,0.02)
if not df_inter.empty:
    list_edges=get_stringDB(df_inter,thres,'symbol')
    st.write(f'retrieved : {len(list_edges)} Interaction with {thres} threshold')

    st.write('## Computing Network')

    import networkx as nx
    H = nx.Graph(list_edges)
    import igraph as ig
    import leidenalg as la
    G = ig.Graph.from_networkx(H)
    partition = la.find_partition(G, la.ModularityVertexPartition)
    subg = partition.subgraphs()
    list_gene=[]
    list_clust=[]
    cluster=96
    # thres_db = col_b.slider("cluster Thresholds", 2, 100,9,1)
    st.write(f'Total Number of clusters before Threshold: {len(subg)}')
    for g in subg:
        cluster=cluster+1
        # print(g.vs['_nx_name'])
        # if len(g.vs['_nx_name'])>thres_db:
        for name in g.vs['_nx_name']:
            list_clust.append(chr(cluster))
            list_gene.append(name)
    df_clust=pd.DataFrame()
    df_clust['symbol']=list_gene
    df_clust['cluster']=list_clust
    if df_clust.empty:
        st.warning("Please reduce the threshold to get more interaction")
        exit()
    if disp:
        st.write("DF_CLUST",df_clust)
    df_umap_cluster=df_inter.merge(df_clust,left_on='symbol',right_on='symbol').reset_index(drop=True)
    # df_umap_cluster['chromosome']=df_umap_cluster['chromosome'].apply(int_to_str)
    if disp:
        st.write(df_umap_cluster)

    ################################### ENRICHMENT ##############################################
    list_cat=get_list_category(df_umap_cluster,'symbol')
    categ = st.selectbox("Select Category", list_cat)
    df_go_ento = get_stringDB_enr(df_umap_cluster,'symbol',categ)
    # df_umap_cluster['chromosome']=df_umap_cluster['chromosome'].apply(int_to_str)
    if disp:
        st.write(df_umap_cluster)
    list_enr={}
    for grpName, rows in df_clust.groupby('cluster'):
        df_temp = get_stringDB_enr(rows['symbol'].unique(),cat=categ)
        # st.write(df_temp['Description'], grpName)
        
        # list_enr['cluster'].append(grpName)
        if len(df_temp['Description'])>0:
            list_enr.update({grpName:df_temp['Description'][0]})
            
        else:
             list_enr.update({grpName:'Null'})
    ################################### ENRICHMENT ##############################################
    

    st.write(f'{categ} Enrichment',df_go_ento)
    df_go_ento['log_p_val']=-np.log10(df_go_ento['p_val'].apply(str_to_float))
    fig_bar=px.bar(df_go_ento,x='Description',y='log_p_val')
    st.plotly_chart(fig_bar)

    import matplotlib.pyplot as plt

    # subax1 = plt.subplot(121)
    st.write('### Displaying the full network/graph clustered with Leiden approach')
    # col1, col2,col3 = st.columns(3)

    vert_size = st.sidebar.slider('vertex size',min_value=0.2,max_value=20.0,step=0.1,value=1.0)
    lab_size = st.sidebar.slider('label size',min_value=0.2,max_value=20.0,step=0.1,value=2.4)
    box_size = st.sidebar.slider('box size',min_value=400,max_value=1600,step=50,value=600)
    ig.config['plotting.backend'] = 'matplotlib'
    subax1=ig.plot(partition,vertex_label=partition.graph.vs['_nx_name'],vertex_label_size=lab_size,vertex_size=vert_size,
                margin=10, bbox=(0,0,box_size, box_size))
    # nx.draw(GG, with_labels=True,node_size=10,font_size=4)
    st.pyplot(subax1.figure,use_container_width=False)
    nbr_of_cluster = ord(df_clust['cluster'].max())-96

    #comment for nothing

    st.write(f'## Computing UMAP with the main {nbr_of_cluster} clusters')
    st.write('To increase or decrease number of clusters please change cluster threshold above')
    import umap
    numerics = ["float16", "float32", "float64"]
    emb = umap.UMAP(random_state=42, verbose=False).fit_transform(df_umap_cluster.select_dtypes(include=numerics))
    df_umap = pd.DataFrame()
    df_umap["X"] = emb[:, 0]
    df_umap["Y"] = emb[:, 1]
    df_umap["target"] = df_umap_cluster['symbol']
    df_umap['size']=5
    df_umap["cluster"] = df_umap_cluster['cluster']
    df_umap['cluster'] = df_umap['cluster'].replace(list_enr)
    if choice=='Cpds':
        df_umap["keggid"] = df_umap_cluster['keggid']
        fig1 = px.scatter(df_umap, x="X", y="Y",color='cluster',text='target',size='size',width=800,height=800,hover_data=['target','keggid'])
    else:
        fig1 = px.scatter(df_umap, x="X", y="Y",color='cluster',text='target',size='size',width=800,height=800,hover_data=['target'])
    st.plotly_chart(fig1, theme="streamlit", use_container_width=True)


    ################################### SIMILARITY ##############################################
    st.write("## Similarity")
    import seaborn as sns
    df_umap_cluster.set_index('symbol',inplace=True)
    lut = dict(zip(set(df_umap_cluster['cluster'].values), sns.hls_palette(len(set(df_umap_cluster['cluster'])), l=0.5, s=0.8)))
    row_colors = pd.DataFrame(df_umap_cluster['cluster'].values)[0].map(lut)

    sns.set(font_scale=0.6)

    figsize = [20, 20]
    fig_clusmap, ax1 = plt.subplots()
    fig_clusmap = sns.clustermap(
                df_umap_cluster.select_dtypes(include=numerics),
                # metric="cosine",
                figsize=figsize,
                row_colors=row_colors.values,
            method="ward",
                xticklabels=False,
                yticklabels=True,
                col_cluster=False,
                cmap="vlag",
                center=0,
                vmin=-5,
                vmax=5,
                # ax=ax1
            )

    st.pyplot(fig_clusmap)

