import gseapy as gp
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from gseapy import barplot, dotplot

st.header("Get biological infos from multiple databases", divider="rainbow")


def convert_df(df):
    return df.to_csv(index=False).encode("utf-8")


def get_infos(dftiti):
    list_cols = df.columns.tolist()

    # dftiti.dropna(subset="symbol", axis=0, inplace=True)

    human = gp.get_library_name(organism="Human")
    col1, col2 = st.columns(2)
    list_cols = dftiti.columns.tolist()
    with col1:
        sel_col = st.selectbox(":red[Choose column] :worried:", list_cols)

    with col2:
        sel = st.selectbox(":green[Choose DB]", human)

    dftiti.dropna(subset=sel_col, axis=0, inplace=True)
    if len(dftiti) > 1:
        # import matplotlib.pyplot as plt

        gene_list = dftiti[sel_col].squeeze().str.strip().to_list()

        enr = gp.enrichr(
            gene_list=gene_list,  # or "./tests/data/gene_list.txt",
            gene_sets=[sel],
            organism="human",  # don't forget to set organism to the one you desired! e.g. Yeast
            outdir=None,  # don't write to disk
        )
        # st.write(gene_list)
        df_enr = enr.results

        df_enr["Log_10_Pv"] = -np.log10(df_enr["Adjusted P-value"])
        # df_enr.rename(columns={"Term": "Pathways"}, inplace=True)
        st.write(df_enr)

        fig = px.bar(df_enr, x="Term", y="Log_10_Pv")
        st.plotly_chart(fig, theme="streamlit", use_container_width=True)
        # fig_c, ax1 = plt.subplots()
        top_t = st.slider("Number of Pathways", 1, 30, 10)
        circle = st.slider("Size of circles", 1, 30, 10)
        ax1 = dotplot(
            enr.res2d,
            title=sel,
            cmap="viridis_r",
            size=circle,
            figsize=(3, 5),
            top_term=top_t,
        )
        # ax1 = dotplot(
        #     enr.results,
        #     column="P-value",
        #     x="Gene_set",  # set x axis, so you could do a multi-sample/library comparsion
        #     size=10,
        #     top_term=5,
        #     figsize=(3, 5),
        #     # title="KEGG",
        #     xticklabels_rot=45,  # rotate xtick labels
        #     show_ring=True,  # set to False to revmove outer ring
        #     marker="o",
        # )
        st.pyplot(ax1.figure)

        ax = barplot(
            enr.res2d, title=sel, figsize=(4, 5), color="darkred", top_term=top_t
        )
        st.pyplot(ax.figure)
    else:
        st.warning("Not enough Data")


def upload_files():
    uploaded_files = st.file_uploader(
        "## Choose files with symbol or geneid in cols, if df null",
        accept_multiple_files=True,
    )
    # st.write(uploaded_file)
    list_df = []
    for uploaded_file in uploaded_files:
        if ".csv" in uploaded_file.name:
            list_df.append(pd.read_csv(uploaded_file))

        if ".fth" in uploaded_file.name:
            list_df.append(pd.read_feather(uploaded_file))
    return list_df


disp_data = st.sidebar.toggle("Display Data")

list_data = []
for key in st.session_state.keys():
    list_data.append(key)
    if disp_data:
        st.write(f"{key}", st.session_state[key].head(2))


if len(list_data) == 0:
    list_df = upload_files()
    if len(list_df) > 0:
        df = pd.concat(list_df)
        st.write("Data from file:", df)
        get_infos(df)

else:
    sel_data = st.selectbox("select your data", list_data)
    df_all = st.session_state[sel_data]
    server_name = st.radio(
        "select server",
        ("all", "pubchem", "KEGG"),
        horizontal=True,
    )

    df = df_all
    if server_name != "all":
        df = df_all[df_all["server"] == server_name]
    st.write("Data :", df)

    if len(df) > 0:
        get_infos(df)
    else:
        st.write("# No data in df: check output of Get Structure, load file(s)")
        list_df = upload_files()

        if len(list_df) > 0:
            df = pd.concat(list_df)
            st.write("Data from file:", df)
            get_infos(df)
