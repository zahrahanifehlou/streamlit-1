import umap
from sklearn.decomposition import PCA
from streamlib import sql_df, convert_df
import sys
import warnings
import numpy as np
import pandas as pd
import psycopg2
import streamlit as st
import plotly.express as px
from sklearn.cluster import KMeans

warnings.filterwarnings("ignore")
sys.path.append("lib/")


def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])


conn = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"
profile_conn = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksilink_cpds",
    password="12345",
)

####################################

# def int_to_str(ints):
#     return str(ints)


# st.header("Display known drugs and/or Crispr profiles",divider='rainbow')

# conn = psycopg2.connect(
#     host="192.168.2.131",
#     port="5432",
#     user="arno",
#     database="ksi_cpds",
#     password="12345",
# )
# sql_kegg = "select cpdpath.pathid,keggcpdgene.geneid, keggcpd.*, cpdbatchs.* from keggcpd\
#             inner join keggcpdgene on keggcpd.keggid=keggcpdgene.keggid\
#             inner join cpd on cpd.keggid=keggcpd.keggid \
#             inner join cpdpath on cpdpath.pubchemid=cpd.pubchemid \
#             inner join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid "

# df_drug_meta = sql_df(sql_kegg, conn)
# df_drug_meta = df_drug_meta.loc[:, ~df_drug_meta.columns.duplicated()].copy()
# df_drug_meta = df_drug_meta.drop_duplicates(subset=["keggid", "source"]).reset_index(drop=True)
# on = st.sidebar.toggle("Display Data")

# # df_drug_meta.dropna(subset="geneid", axis=0, inplace=True)
# if on:
#     st.write("Sample of Data", df_drug_meta.sample(5))
# # # st.write(df_drug_meta.describe())

# list_sources = df_drug_meta["source"].unique().tolist()
# choix_source = st.selectbox("Select the Source", list_sources)

# df_cpd_src = df_drug_meta[df_drug_meta["source"] == choix_source].reset_index(drop=True)
# b_list2 = df_cpd_src["batchid"].to_list()
# bq = []
# for bs in b_list2:
#     bq.append("'" + bs + "'")

# conn2 = psycopg2.connect(
#     host="192.168.2.131",
#     port="5432",
#     user="arno",
#     database="ksilink_cpds",
#     password="12345",
# )
# sql_profile = f"select * from aggprofile where metabatchid  in (" + ",".join(bq) + ")"

# df_cpd_prof = sql_df(sql_profile, conn2)
# df_cpd_prof = df_cpd_prof[df_cpd_prof["metasource"] == choix_source].reset_index(drop=True)
# if on:
#     st.write(df_cpd_prof.sample(5))

# # # st.write(df_cpd_src.sample(5))

# df_cpds_merge = df_cpd_prof.merge(df_cpd_src, left_on="metabatchid", right_on="batchid").reset_index(drop=True)

# if on:
#     st.write("Merged", df_cpds_merge.sample(5))

# # ################################# CRISPR ######################################


# dict1=df_cpds_merge.set_index("geneid").to_dict()["pathid"]
# dict_eff=df_cpds_merge.set_index("geneid").to_dict()["efficacy"]
# st.write("## Loading CRispr")
# b_list2 = df_cpds_merge["geneid"].to_list()
# bq = []
# for bs in b_list2:
#     bq.append("'" + bs + "'")

# sql_crispr = "select * from genebatchs where geneid  in (" + ",".join(bq) + ")"
# df_crispr = sql_df(sql_crispr, conn)

# df_crispr["pathid"]=df_crispr["geneid"].map(dict1)
# df_crispr["efficacy"]=df_crispr["geneid"].map(dict_eff)

# sql_genes = "select synonyms,geneid,chromosome,locus,mainlocation,symbol from gene where geneid  in (" + ",".join(bq) + ")"
# df_genes = sql_df(sql_genes, conn)
# df_genes['chromosome']=df_genes['chromosome'].apply(int_to_str)
# if on:
#     st.write("GeneInfos", df_genes.sample(5))

# # ######################################### DONE ##########################################

# sql_path = "select geneid,pathid from genepath where geneid  in (" + ",".join(bq) + ")"
# df_paths = sql_df(sql_path, conn)

# dict_path={k: [g['pathid'].tolist()] for k,g in df_paths.groupby('geneid')}
# df_paths2=pd.DataFrame(list(dict_path.items()),columns=['geneid','pathids'])
# # st.write(df_paths2)

# # dict_genes={k: [g['geneid'].tolist()] for k,g in df_paths.groupby('pathid')}
# # df_genes2=pd.DataFrame(list(dict_genes.items()),columns=['paths','geneids'])
# # st.write(df_genes2)

# # ######################################### UNTIL HERE ##########################################


# b_list2 = df_crispr["batchid"].to_list()
# bq = []
# for bs in b_list2:
#     bq.append("'" + bs + "'")

# sql_crispr_prof = "select * from aggprofile where metabatchid  in (" + ",".join(bq) + ")"
# df_crispr_prof = sql_df(sql_crispr_prof, conn2)
# df_crispr_merge = df_crispr_prof.merge(df_crispr, left_on="metabatchid", right_on="batchid").reset_index(drop=True)

# df_crispr_merge_infos=df_crispr_merge.merge(df_genes, left_on="geneid", right_on="geneid").reset_index(drop=True)

# # st.write("CRISPR", df_crispr_merge)
# # st.write("CRISPR INFOS", df_crispr_merge_infos)

# # ################################# MERGE CPDS INFOS ######################################

# df_cpds_merge_infos=df_cpds_merge.merge(df_genes, left_on="geneid", right_on="geneid").reset_index(drop=True)

# # ################################# MERGE with pathways ######################################

# df_crispr_merge_infos_all=df_crispr_merge_infos.merge(df_paths2, left_on="geneid", right_on="geneid").reset_index(drop=True)
# df_cpds_merge_infos_all=df_cpds_merge_infos.merge(df_paths2, left_on="geneid", right_on="geneid").reset_index(drop=True)
# # st.write("CRISPR INFOS", df_crispr_merge_infos_all)
# # st.write("CPDS INFOS", df_cpds_merge_infos_all)

# # ################################# COMPUTATION ######################################

# st.write("## Computation part")
# # st.write("CPDS INFOS", df_cpds_merge_infos)
# numerics = ["float16", "float32", "float64"]
# df_num_crispr = df_crispr_merge_infos_all.select_dtypes(include=numerics)
# df_num_cpd = df_cpds_merge.select_dtypes(include=numerics)
# list_common_cols = list(set(df_num_crispr.columns).intersection(df_num_cpd.columns))
# df_cpds_merge.dropna(subset='efficacy',inplace=True)

# list_meta = ["geneid","pathid", "chromosome", "synonyms","mainlocation",'pathids','symbol']
# df_num_cpd[list_meta]=df_cpds_merge_infos_all[list_meta]
# df_num_cpd["efficacy"] = df_cpds_merge["efficacy"].apply(
#     lambda x: x.split(",")[0]
# )
# df_num_cpd["Type"] = "CPDS"


# df_num_crispr["Type"] = "CRISPR"
# df_num_crispr[list_meta]=df_crispr_merge_infos_all[list_meta]
# df_num_crispr["efficacy"] = df_crispr_merge_infos_all["efficacy"].apply(
#     lambda x: x.split(",")[0])


# list_common_cols2 = list(set(df_num_crispr.columns).intersection(df_num_cpd.columns))
# df_cpd_gene = pd.concat([df_num_crispr[list_common_cols2], df_num_cpd[list_common_cols2]]).reset_index(drop=True)

# # combat = Combat()
# # Y_align = combat.fit_transform(Y=df_cpd_gene.select_dtypes(include=numerics).values, b=df_cpd_gene['Type'].values)
# # df_combat = pd.DataFrame(Y_align,columns=df_cpd_gene.select_dtypes(include=numerics).columns)

# # df_combat["Type"] = df_cpd_gene["Type"]

# # df_combat["efficacy"] = df_cpd_gene["efficacy"]
# # df_combat[list_meta]=df_crispr_merge_infos_all[list_meta]


# col1, col2 = st.columns(2)

# list_df = ["cpds only", "crispr only", "cpds vs crispr"]

# sel_data = col1.selectbox("Select the data", list_df)

# if sel_data == "cpds only":
#     df = df_num_cpd
# elif sel_data == "crispr only":
#     df = df_num_crispr
# elif sel_data == "cpds vs crispr":
#     df = df_cpd_gene
# # elif sel_data=="cpds vs crispr aligned":
# #     df=df_combat

# st.write(sel_data, len(df), len(df.columns))
# ###################################################################################################################
# ###################################################################################################################
# # st.write("## StringDB" )
# # st.write('DF',df)# .dropna(subset='symbol',inplace=True)

# # import networkx as nx
# # list_edges = get_stringDB(df,0.4,'symbol')
# # H = nx.Graph(list_edges) #Get networkx graph
# # import igraph as ig
# # import leidenalg as la
# # G = ig.Graph.from_networkx(H) # Convert to igraph
# # partition = la.find_partition(G, la.ModularityVertexPartition) # Leiden Algorithm
# # ################### Create cluster Names #######################
# # subg = partition.subgraphs()
# # list_gene=[]
# # list_clust=[]
# # cluster=96
# # for g in subg:
# #     cluster=cluster+1
# #     # print(g.vs['_nx_name'])
# #     # if len(g.vs['_nx_name'])>9:
# #     for name in g.vs['_nx_name']:
# #         list_clust.append(chr(cluster))
# #         list_gene.append(name)

# # df_clust=pd.DataFrame()
# # df_clust['symbol']=list_gene
# # df_clust['cluster']=list_clust


# # df_umap_cluster=df.merge(df_clust,left_on='symbol',right_on='symbol').reset_index(drop=True)

# # st.write(df['pathids'])
# # st.write(sel_data,df.select_dtypes(include=numerics).sample(5))
# st.write("## UMAP" )

# if on:
#     st.write("Data",df)
# list_opt = ["geneid","pathid", "chromosome", "synonyms","mainlocation",'Type','efficacy','symbol']

# sel_opt = col2.selectbox("Select the coloring", list_opt)
# reducer = umap.UMAP(random_state=42, verbose=False,learning_rate=0.1)
# embedding = reducer.fit_transform(df.select_dtypes(include=numerics))
# df_emb = pd.DataFrame()
# df_emb["x"] = embedding[:, 0]
# df_emb["y"] = embedding[:, 1]
# df_emb[list_opt]=df[list_opt]
# df_emb["size"] = 5

# fig1 = px.scatter(df_emb, x="x", y="y", hover_data=["synonyms","geneid","pathid","chromosome",'symbol'],
#                    color=sel_opt,text=df_emb["geneid"],size='size')#,symbol='meta_target')#,text=df_emb[sel_opt])
# st.plotly_chart(fig1, theme="streamlit", use_container_width=False)


# st.write("## Similarity")

# import matplotlib.pyplot as plt
# import seaborn as sns

# df.set_index(sel_opt,inplace=True)
# # plt_src,col_colors=get_col_colors(tmp)
# sns.set(font_scale=0.5)
# fig_clusmap, ax1 = plt.subplots()
# fig_clusmap = sns.clustermap(
#             df.select_dtypes(include=numerics),
#             metric="cosine",
#             # col_colors=col_colors,
#         # method="ward",
#             xticklabels=False,
#             yticklabels=True,
#             col_cluster=False,
#             cmap="vlag",
#             center=0,
#             vmin=-5,
#             vmax=5,
#         )

# st.pyplot(fig_clusmap)
