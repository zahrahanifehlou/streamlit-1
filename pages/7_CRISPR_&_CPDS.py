import warnings
import pandas as pd
import streamlit as st
import umap
import plotly.express as px
from sklearn.metrics.pairwise import cosine_similarity

pd.set_option("display.max_rows", 100)
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.filterwarnings("ignore")
import psycopg2

st.write("## Based on Kegg Annotations....")

if "df_profiles" and  "df_crisper" not in st.session_state:
    st.write("Connect DB First")
else:
    df_source = st.session_state["df_profiles"]
    df_crisper = st.session_state["df_crisper"]
    df_efficacy = st.session_state["df_efficacy"]
    df_cpdGene = st.session_state["df_cpdGene"]
    df_cpdGene = df_cpdGene[df_cpdGene["server"] == "KEGG"]
    

    list_sources = df_source["metasource"].unique().tolist()
    st.write("Data from DB", df_source)
    choix_source = st.selectbox("Select the Source", list_sources)
    df_source = df_source[df_source["metasource"] == choix_source]
    df_source=df_source[df_source["metageneid"].notna()].reset_index(drop=True)
    

    tab1, tab2 = st.tabs(["compounds profile", "Crisper profile"])
    tab1.write(df_source)
    tab2.write(df_crisper)

    def find_umap(df, title):
        filter_cols = [col for col in df.columns if not col.startswith("meta")]
        reducer = umap.UMAP(densmap=True, random_state=42, verbose=True)
        embedding = reducer.fit_transform(df[filter_cols])
        df_emb = pd.DataFrame({"x": embedding[:, 0], "y": embedding[:, 1]})
        df_emb["target"] = df["metageneid"]
        df_emb["Type"] = df["metatype"]
        df_emb["efficacy"] = df["metaefficacy"]
        fig = px.scatter(df_emb, x="x", y="y", hover_data=["target"], color="efficacy", title=title)
        st.plotly_chart(fig, theme="streamlit", use_container_width=True)

    # cols = st.columns(3)
    # with cols[0]:
    #     find_umap(df_source, "UMAP in CPD profile")
    # with cols[1]:
    #     find_umap(df_crisper, "UMAP in Crisper profile")
    # with cols[2]:
    #     df = pd.concat([df_source, df_crisper]).reset_index(drop=True)
    #     find_umap(df, "UMAP in CPD and Crisper profile")

    # st.write("## UMAP")

    
    # st.write("## Similarity")
    # filter_cols = [col for col in df.columns if not col.startswith("meta")]
    # simi = cosine_similarity(df[filter_cols])

    # sim_df = pd.DataFrame(simi)
    # sim_df['target']=df['metageneid'] +' _ ' + df['metatype']
    # st.write('sim_df',sim_df)
    # fig = px.imshow(sim_df, title="Profiles Similarities")
    # st.plotly_chart(fig)

