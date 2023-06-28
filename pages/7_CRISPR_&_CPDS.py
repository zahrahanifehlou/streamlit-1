import warnings
from functools import partial
from pydoc import describe

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import scipy.cluster.hierarchy as sch
import streamlit as st
from plotly.figure_factory import create_dendrogram
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics.pairwise import cosine_similarity

pd.set_option("display.max_rows", 100)
pd.set_option("display.max_columns", 4000)
pio.templates.default = "plotly_dark"
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.filterwarnings("ignore")
import psycopg2
import umap

st.write("## Based on Kegg Annotations....")

if "df_profiles" not in st.session_state:
    st.write("Connect DB First")
else:
    all_df = st.session_state["df_profiles"]
    crisper_pro = st.session_state["df_crisper"]
    list_sources = all_df["metasource"].unique().tolist()
    st.write("Data from DB", all_df)
    choix_source = st.selectbox("Select the Source", list_sources)
    df_source = all_df[all_df["metasource"] == choix_source]
    tab1, tab2 = st.tabs(["compounds profile", "Crisper profile"])
    tab1.write(df_source)
    tab2.write(crisper_pro)


    ################################# COMPUTATION ######################################
    st.write('## Computation part')

    numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
    df_num_crispr = df_crispr_merge.select_dtypes(include=numerics)
    df_num_cpd = df_cpds_merge_meta.select_dtypes(include=numerics)
    list_common_cols = list(set(df_num_crispr.columns).intersection(df_num_cpd.columns))
    df_cpds_merge_meta.dropna(subset='efficacy',inplace=True)

    df_num_cpd["meta_efficacy"] = df_cpds_merge_meta["efficacy"].apply(
        lambda x: x.split(",")[0]
    )
    df_num_cpd["meta_target"] = df_cpds_merge_meta["geneid"]+'_'+df_cpds_merge_meta['keggid']
    df_num_crispr["meta_target"] = df_crispr_merge["geneid"]
    df_num_cpd["meta_Type"] = "CPDS"
    df_num_crispr["meta_Type"] = "CRISPR"
    df_num_crispr["meta_efficacy"] = "Unknown"


    # st.write('join',df_join_kegg_geneid_prof_filt)
    list_common_cols2 = list(set(df_num_crispr.columns).intersection(df_num_cpd.columns))
    # st.write(list_common_cols2)
    df_cpd_gene = pd.concat(
        [df_num_crispr[list_common_cols2], df_num_cpd[list_common_cols2]]
    ).reset_index(drop=True)
    # st.write('CPD_GENES',df_cpd_gene)
    list_df = ["cpds only", "crispr only", "cpds vs crispr"]

    sel_data = st.selectbox("Select the data", list_df)

    if sel_data == "cpds only":
        df = df_num_cpd
    elif sel_data == "crispr only":
        df = df_num_crispr
    elif sel_data == "cpds vs crispr":
        df = df_cpd_gene

    st.write('df',len(df),len(df.columns))

    # df=df.dropna(subset='Dr_Gene')
    # st.write('df',df['target'])
    cols = [c for c in list_common_cols2 if not c.startswith("meta")]
    # df_num_cpd.to_csv('625.csv')
    # df.drop(columns=["Unnamed: 0"],inplace=True)
    # df = df[cols]
    st.write("## UMAP" )

    reducer = umap.UMAP(densmap=True, random_state=42, verbose=True)
    embedding = reducer.fit_transform(df.select_dtypes(include=numerics))
    df_emb = pd.DataFrame()
    df_emb["x"] = embedding[:, 0]
    df_emb["y"] = embedding[:, 1]
    df_emb["target"] = df["meta_target"]
    df_emb["Type"] = df["meta_Type"]
    df_emb["efficacy"] = df["meta_efficacy"]
    fig1 = px.scatter(df_emb, x="x", y="y", hover_data=["target"], color="efficacy")
    st.plotly_chart(fig1, theme="streamlit", use_container_width=True)



    st.write("## Similarity")
    simi = cosine_similarity(df[cols])

    sim_df = pd.DataFrame(simi)
    sim_df['target']=df['meta_target']
    st.write('sim_df',sim_df)
    fig = px.imshow(sim_df, title="Profiles Similarities")
    st.plotly_chart(fig)

    # pw_func = partial(pdist, metric="euclidean")
    # fig = create_dendrogram(
    #     df_emb.select_dtypes(include=numerics),
    #     orientation="left",
    #     distfun=pw_func,
    #     linkagefun=lambda x: sch.linkage(x, "ward"),
    #     labels=df_emb.target.values,
    # )
    # # fig.update_layout(width=2000, height=1600,template='plotly_dark')
    # st.plotly_chart(fig, theme="streamlit", use_container_width=True)
