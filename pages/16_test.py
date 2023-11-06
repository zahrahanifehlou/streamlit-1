import streamlit as st
from streamlit_plotly_events import plotly_events
import numpy as np
import plotly.express as px
from streamlib import sql_df
import psycopg2
import pandas as pd
from sklearn.cluster import KMeans
profile_conn = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksilink_cpds",
    password="12345",
)


conn = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksi_cpds",
    password="12345",
)


def toto():
    sql_umqpemd =  f"SELECT * FROM aggcombatprofile where metasource='Ksilink_625'"
    df_src_emd = sql_df(sql_umqpemd, profile_conn)
    return df_src_emd

df=toto()
all_cpds=["Staurosporine","Saccharin","Sorbitol",   "CA-074","BAYK8644",
    "Lys05","Cucurbitacin","FCCP","Rapamycin","Cladribine","Cytarabine",
    "Etoposide","Berberine","Fluphenazine","Docetaxel","Oxibendazole",
    "Ethoxyquin","Rotenone","GSK2879552","BAY K8644","NMDA","Tetrandrine",
    'Aloxistatin','Dexamethasone','Quinidine','LY2109761','AMG900','NVS-PAK1-1','Mirk-IN-1','Daporinad']
df_dcps=pd.DataFrame()
df_dcps['cpds']=all_cpds
df_dcps.to_csv('/mnt/shares/L/Temp/cpds_tox.csv',index=None)
cols = [c for c in df.columns if  not c.startswith("meta")]
meta_cols = [c for c in df.columns if  c.startswith("meta")]
X = df[cols]
# KNN
nb_cluster=st.slider('Number of clusters',min_value=2,max_value=30,value=15,step=1)
#kmeans = KMeans(n_clusters=nb_cluster, random_state=0).fit(X)


## PCA
# from sklearn.decomposition import PCA
from cuml import PCA
pca_2d = PCA(n_components=2)
projection_2d = pca_2d.fit_transform(X)
emd_pca = pd.DataFrame()
emd_pca["pca1"] = projection_2d[:, :1].flatten()
emd_pca["pca2"] = projection_2d[:, 1:2].flatten()
kmeans = KMeans(n_clusters=nb_cluster, random_state=0).fit(emd_pca)
# emd_pca[meta_cols] = df[meta_cols]

emd_pca['Cluster']=kmeans.labels_
emd_pca['Cluster'] = emd_pca['Cluster'].astype(str)
fig1 = px.scatter(
                emd_pca,
                x="pca1",
                y="pca2",
                color="Cluster",
            
               
                title=f" PCA ",
                # hover_data=["metabatchid"],
            )
            
st.plotly_chart(fig1, theme="streamlit",
                                use_container_width=True)

## UMAP
import umap

umap_2d = umap.UMAP(n_components=2, n_neighbors=30,
                min_dist=0, verbose=True, random_state=42)
mapper=umap_2d.fit(X)
projection_2d = umap_2d.fit_transform(X)
emd_umap = pd.DataFrame()
emd_umap["umap1"] = projection_2d[:, :1].flatten()
emd_umap["umap2"] = projection_2d[:, 1:2].flatten()
# emd_umap[meta_cols] = df[meta_cols]
emd_umap['Cluster']=kmeans.labels_
emd_umap['Cluster'] = emd_umap['Cluster'].astype(str)
fig = px.scatter(
                emd_umap,
                x="umap1",
                y="umap2",
                color="Cluster",
              
                
                title=f" UMAP ",
                # hover_data=["metabatchid"],
            )
            
st.plotly_chart(fig, theme="streamlit",
                                use_container_width=True)


st.write(df)
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
df_g1=emd_pca.groupby("Cluster").mean()
df_g2=emd_umap.groupby("Cluster").mean()
st.write(df_g1)
st.write(df_g2)
fig3=ff.create_quiver( df_g1["pca1"],df_g1["pca2"],(df_g2["umap1"]-df_g1["pca1"]),df_g2["umap2"]-df_g1["pca2"],scale=1)
st.plotly_chart(fig3)
import umap.plot
fig6 = umap.plot.diagnostic(mapper, diagnostic_type='pca')
# st.write(fig6)
st.pyplot(fig6.axes.figure)

fig7=umap.plot.diagnostic(mapper, diagnostic_type='vq')
st.pyplot(fig7.axes.figure)

fig8=umap.plot.diagnostic(mapper, diagnostic_type='local_dim')
st.pyplot(fig8.axes.figure)

fig9=umap.plot.diagnostic(mapper, diagnostic_type='neighborhood')
st.pyplot(fig9.axes.figure)

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