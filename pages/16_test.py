import streamlit as st
from streamlit_plotly_events import plotly_events
import numpy as np
import plotly.express as px
from streamlib import sql_df,find_sim_cpds
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import cosine_similarity


conn_meta = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"
conn_prof = "postgres://arno:12345@192.168.2.131:5432/ksilink_cpds"

st.title('Experimental: works only on 131')
sql_rep='select symbol1 from crisprcos'
df_rep= sql_df(sql_rep,conn_prof)
df_rep = df_rep.drop_duplicates().reset_index(drop=True)
st.write(df_rep)
bq = []
for bs in df_rep['symbol1']:
    bq.append("'" + bs + "'")
    
sql_umqpemd =  f"SELECT * FROM aggcombatprofile where metasource='CRISPER' and metabatchid  in (" + ",".join(bq) + ") "
df_src_emd = sql_df(sql_umqpemd, conn_prof)
df_sel = df_src_emd[df_src_emd["metabatchid"]=='DMSO']
sim_crispr = find_sim_cpds(df_src_emd, df_sel)
df_hist_cpd = pd.DataFrame(
    {"sim": sim_crispr.flatten().tolist(
    ), "metabatchid": df_src_emd["metabatchid"]}
)
df_hist_cpd['metabatchid']=df_hist_cpd['metabatchid'].str.split('_').str[0]

df_rna=pd.read_csv('Recursion_U2OS_expression_data.csv')
dict_rna = df_rna.set_index('gene').to_dict()['zfpkm']

df_hist_cpd['zfpkm']=df_hist_cpd['metabatchid'].map(dict_rna)
df_hist_cpd.dropna(subset='zfpkm', inplace=True)
df_hist_cpd=df_hist_cpd[df_hist_cpd['zfpkm']>-30]
df_hist_cpd=df_hist_cpd[df_hist_cpd['zfpkm']<-3]
st.write(df_hist_cpd)
fig = px.scatter(df_hist_cpd,x='sim',y='zfpkm',hover_data=['metabatchid'])
st.plotly_chart(fig)


# def toto():
#     sql_umqpemd =  f"SELECT * FROM aggcombatprofile where metasource='CRISPER'"
#     df_src_emd = sql_df(sql_umqpemd, conn_prof)
#     return df_src_emd


# df_dcp= pd.read_csv('/mnt/shares/L/Temp/Test/test.csv')
# # numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
# # cols_num = df_dcp.select_dtypes(include=numerics).columns
# # st.dataframe(df_dcp[cols_num].style.format('{:.4g}'))
# # st.write(df_dcp.dtypes)
# st.dataframe(df_dcp,column_config={"EC_50_Nuclei_Tot":st.column_config.NumberColumn("EC_50_Nuclei_Tot", format="%7.2g",)})
# numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
# cols_num = df_dcp.select_dtypes(include=numerics).columns
# for col in cols_num:
#     df_dcp[col] = df_dcp[col].apply(lambda x: "{:.2e}".format(x))



# st.dataframe(df_dcp)

# df=toto()
# all_cpds=["Staurosporine","Saccharin","Sorbitol",   "CA-074","BAYK8644",
#     "Lys05","Cucurbitacin","FCCP","Rapamycin","Cladribine","Cytarabine",
#     "Etoposide","Berberine","Fluphenazine","Docetaxel","Oxibendazole",
#     "Ethoxyquin","Rotenone","GSK2879552","BAY K8644","NMDA","Tetrandrine",
#     'Aloxistatin','Dexamethasone','Quinidine','LY2109761','AMG900','NVS-PAK1-1','Mirk-IN-1','Daporinad']
# df_dcps=pd.DataFrame()
# df_dcps['cpds']=all_cpds
# df_dcps.to_csv('/mnt/shares/L/Temp/cpds_tox.csv',index=None)
# cols = [c for c in df.columns if  not c.startswith("meta")]
# meta_cols = [c for c in df.columns if  c.startswith("meta")]
# X = df[cols]
# # KNN
# nb_cluster=st.slider('Number of clusters',min_value=2,max_value=30,value=15,step=1)
# #kmeans = KMeans(n_clusters=nb_cluster, random_state=0).fit(X)


# ## PCA
# # from sklearn.decomposition import PCA
# from cuml import PCA
# pca_2d = PCA(n_components=2)
# projection_2d = pca_2d.fit_transform(X)
# # st.write(projection_2d[0])
# emd_pca = pd.DataFrame()
# emd_pca["pca1"] = projection_2d[0]
# emd_pca["pca2"] = projection_2d[1]
# kmeans = KMeans(n_clusters=nb_cluster, random_state=0).fit(emd_pca)
# # emd_pca[meta_cols] = df[meta_cols]

# emd_pca['Cluster']=kmeans.labels_
# emd_pca['Cluster'] = emd_pca['Cluster'].astype(str)
# fig1 = px.scatter(
#                 emd_pca,
#                 x="pca1",
#                 y="pca2",
#                 color="Cluster",
            
               
#                 title=f" PCA ",
#                 # hover_data=["metabatchid"],
#             )
            
# st.plotly_chart(fig1, theme="streamlit",
#                                 use_container_width=True)

# ## UMAP
# # import umap
# from cuml import UMAP
# umap_2d = UMAP(n_components=2, n_neighbors=30,
#                 min_dist=0, verbose=True, random_state=42)
# mapper=umap_2d.fit(X)
# projection_2d = mapper.transform(X)
# emd_umap = pd.DataFrame()
# emd_umap["umap1"] = projection_2d[0]
# emd_umap["umap2"] = projection_2d[1]
# # emd_umap[meta_cols] = df[meta_cols]
# emd_umap['Cluster']=kmeans.labels_
# emd_umap['Cluster'] = emd_umap['Cluster'].astype(str)
# fig = px.scatter(
#                 emd_umap,
#                 x="umap1",
#                 y="umap2",
#                 color="Cluster",
              
                
#                 title=f" UMAP ",
#                 # hover_data=["metabatchid"],
#             )
            
# st.plotly_chart(fig, theme="streamlit",
#                                 use_container_width=True)
# from cuml.metrics import trustworthiness
# cu_score=trustworthiness(X,projection_2d)
# st.write(" cuml's trustworthiness score : ",cu_score)

# st.write(df)
# import matplotlib.pyplot as plt
# import plotly.figure_factory as ff
# df_g1=emd_pca.groupby("Cluster").mean()
# df_g2=emd_umap.groupby("Cluster").mean()
# st.write(df_g1)
# st.write(df_g2)
# fig3=ff.create_quiver( df_g1["pca1"],df_g1["pca2"],(df_g2["umap1"]-df_g1["pca1"]),df_g2["umap2"]-df_g1["pca2"],scale=1)
# st.plotly_chart(fig3)
# import umap.plot
# fig6 = umap.plot.diagnostic(mapper, diagnostic_type='pca')
# # st.write(fig6)
# st.pyplot(fig6.axes.figure)

# fig7=umap.plot.diagnostic(mapper, diagnostic_type='vq')
# st.pyplot(fig7.axes.figure)

# fig8=umap.plot.diagnostic(mapper, diagnostic_type='local_dim')
# st.pyplot(fig8.axes.figure)

# fig9=umap.plot.diagnostic(mapper, diagnostic_type='neighborhood')
# st.pyplot(fig9.axes.figure)

# fig10 = umap.plot.connectivity(mapper, edge_bundling='hammer')
# st.pyplot(fig10.axes.figure)
# X=df_g1.values

# grid_unit_size = .02
# margin = 0.25

# x_min = X[:, 0].min() - margin
# x_max = X[:, 0].max() + margin
# y_min = X[:, 1].min() - margin
# y_max = X[:, 1].max() + margin

# x_range = np.arange(start = x_min, stop = x_max, step = grid_unit_size )
# y_range = np.arange(start = y_min, stop = y_max, step = grid_unit_size )

# x_gridvalues, y_gridvalues = np.meshgrid(x_range, y_range)

# # COMPUTE PROBABILITIES ON GRID
# gridvalues_combined_tidy = np.vstack([x_gridvalues.flatten(), y_gridvalues.flatten()]).T
# # knn_class_probabilities = knn_classifier.predict_proba(gridvalues_combined_tidy)

# probability_postive_class = knn_class_probabilities[:,1]

# fig = go.Figure(data=[
#     go.Contour(
#         x = x_range
#         ,y = y_range
#         ,z = probability_postive_class.reshape(x_gridvalues.shape)
#         ,colorscale = 'RdBu'
#         #,alpha = .3
#     )
# ])
 
# x = np.linspace(-4, 4, 9)
# y = np.linspace(-5, 5, 11)
# random_data = np.random.random((11, 9))
# x_1, y_1 = np.meshgrid(df_g1["pca1"], df_g1["pca2"])
# data1,data2 = np.meshgrid(df_g2["umap1"],df_g2["umap2"])

# # data= np.concatenate(df_g2["umap1"].values,df_g2["umap2"].values)
# st.write(data1)
# fig6 = plt.plot_surface(x_1,y_1,data1, cmap = 'jet')

# st.pyplot(fig6.axes.figure)