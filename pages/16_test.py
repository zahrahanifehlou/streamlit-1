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
cols = [c for c in df.columns if  not c.startswith("meta")]
meta_cols = [c for c in df.columns if  c.startswith("meta")]
X = df[cols]
# KNN
nb_cluster=st.slider('Number of clusters',min_value=2,max_value=30,value=15,step=1)
#kmeans = KMeans(n_clusters=nb_cluster, random_state=0).fit(X)


## PCA
from sklearn.decomposition import PCA
pca_2d = PCA(n_components=2)
projection_2d = pca_2d.fit_transform(X)
emd_pca = pd.DataFrame()
emd_pca["pca1"] = projection_2d[:, :1].flatten()
emd_pca["pca2"] = projection_2d[:, 1:2].flatten()
kmeans = KMeans(n_clusters=nb_cluster, random_state=0).fit(emd_pca)
emd_pca[meta_cols] = df[meta_cols]

emd_pca['Cluster']=kmeans.labels_
emd_pca['Cluster'] = emd_pca['Cluster'].astype(str)
fig1 = px.scatter(
                emd_pca,
                x="pca1",
                y="pca2",
                color="Cluster",
            
               
                title=f" PCA ",
                hover_data=["metabatchid"],
            )
            
st.plotly_chart(fig1, theme="streamlit",
                                use_container_width=True)

## UMAP
from umap import UMAP

umap_2d = UMAP(n_components=2, n_neighbors=30,
                min_dist=0, verbose=True, random_state=42)

projection_2d = umap_2d.fit_transform(X)
emd_umap = pd.DataFrame()
emd_umap["umap1"] = projection_2d[:, :1].flatten()
emd_umap["umap2"] = projection_2d[:, 1:2].flatten()
emd_umap[meta_cols] = df[meta_cols]
emd_umap['Cluster']=kmeans.labels_
emd_umap['Cluster'] = emd_umap['Cluster'].astype(str)
fig = px.scatter(
                emd_umap,
                x="umap1",
                y="umap2",
                color="Cluster",
              
                
                title=f" UMAP ",
                hover_data=["metabatchid"],
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