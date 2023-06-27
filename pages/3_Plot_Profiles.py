import pandas as pd
import plotly.express as px
import pqdm
import psycopg2
import streamlit as st
from pqdm.processes import pqdm
from sklearn.metrics.pairwise import cosine_similarity

def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])
conn = init_connection()
profile_conn = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksilink_cpds",
    password="12345",
)

if "df_profiles" not in st.session_state:
    st.write("Connect DB First")
else:
    all_df = st.session_state["df_profiles"]
    list_sources = all_df["metasource"].unique().tolist()

    st.write("Data from DB", all_df)

    choix_source = st.selectbox("Select the Source", list_sources)
    batchs = all_df[all_df["metasource"] == choix_source]["metabatchid"].values

    sql_profile = (
        "select * from aggcombatprofile where metasource=" + "'" + choix_source + "'"
    )
    df_source = pd.read_sql(sql_profile, profile_conn)

    sql_umqpemd = "select * from umapemd where metasource=" + "'" + choix_source + "'"
    df_src_emd = pd.read_sql(sql_umqpemd, profile_conn)
    df_src_emd["color"] = "others"
    df_src_emd.loc[
        df_src_emd["metabatchid"].isin(batchs), "color"
    ] = "selected compounds"

    fig = px.scatter(
        df_src_emd,
        x="umap1",
        y="umap2",
        color="color",
        title=f"{choix_source} profiles ",
        hover_data=["metabatchid", "metaplate", "metawell"],
    )
    st.plotly_chart(fig, theme="streamlit", use_container_width=True)

    col3, col4 = st.columns(2)
    with col3:
        choix = st.selectbox("Select the Profile", batchs)

    df_sel = all_df[
        (all_df["metasource"] == choix_source) & (all_df["metabatchid"] == choix)
    ]
    with col4:
        st.write("Selected Profile", df_sel)

    list_cpd = []
    progress_text = "Computing similarities. Please wait"
    list_keep = []
    list_cpd_keep = []

    # @st.cache_data
    def find_sim_cpds(df1, df2):
        filter_col1 = [col for col in df1.columns if not col.startswith("meta")]
        filter_col2 = [col for col in df2.columns if not col.startswith("meta")]
        filter_col = list(set(filter_col1) & set(filter_col2))
        simi = cosine_similarity(df1[filter_col], df2[filter_col])
        return simi

    sim_cpds = find_sim_cpds(df_source, df_sel)
    df_hist = pd.DataFrame()
    df_hist["sim"] = sim_cpds.flatten().tolist()
    df_hist["metabatchid"] = df_source["metabatchid"]
    col3, col4 = st.columns(2)
    with col3:
        thres = st.slider("Threshold", -1.0, 1.0, 0.85)
    with col4:
        thresq = st.slider("Cardinal Threshold", 0, 200, 10)
    fig = px.histogram(df_hist, x="sim")
    fig.add_vline(x=thres)

    df_keep = df_hist[df_hist["sim"] > thres]
    df_keep = df_keep.head(thresq).reset_index()
    df_keep = df_keep.sort_values(by="sim", ascending=False)
    df_keep.reset_index(drop=True, inplace=True)

    batch_list = df_keep["metabatchid"].tolist()
    b_list = []
    df_cpd = pd.DataFrame()
    for b in batch_list:
        if "jcp2022_800" not in b:
            b_list.append("'" + b + "'")

    if len(b_list) > 0:
        sql_meta = (
            "select cpd.pubchemid,cpd.keggid, cpd.name, cpd.smile from cpd inner join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid where cpdbatchs.batchid in "
            + "("
            + ",".join(b_list)
            + ")"
        )

        df_cpd = pd.read_sql(sql_meta, conn)
        

        if len(df_cpd) > 0:
            df_cpd.fillna("No result", inplace=True)
            df_cpd = df_cpd[df_cpd["smile"] != "No result"].reset_index()
            df_cpd = df_cpd.drop(columns=["index"])

            tab_list = st.tabs(["Histogram", "Data", "UMAP", "MetaData", "Summary"])
            with tab_list[0]:
                st.plotly_chart(fig, theme="streamlit", use_container_width=True)
            with tab_list[1]:
                st.write(df_keep)
            with tab_list[2]:
                df_src_emd["color"] = "others"

                df_src_emd.loc[
                    df_src_emd["metabatchid"].isin(batch_list), "color"
                ] = "similar compounds"
                df_src_emd.loc[
                    df_src_emd["metabatchid"] == choix, "color"
                ] = "selectd compounds"
                fig = px.scatter(
                    df_src_emd,
                    x="umap1",
                    y="umap2",
                    color="color",
                    title=f"similar cpds to {choix} profiles  ",
                    hover_data=["metabatchid", "metaplate", "metawell"],
                )
                st.plotly_chart(fig, theme="streamlit", use_container_width=True)

            with tab_list[3]:
                st.write(df_cpd)

            with tab_list[4]:
                st.write(df_keep.describe())

            df_keep_prof = df_source[
                df_source["metabatchid"].isin(df_keep["metabatchid"].values)
            ]
            st.write(df_keep_prof)

            if len(df_keep_prof) < 11:
                cpd_names = df_keep_prof.metabatchid.values
                df_plt = df_keep_prof.set_index("metabatchid")
                filter_col = [
                    col for col in df_plt.columns if not col.startswith("meta")
                ]
                df_plt = df_plt[filter_col].T
                fig2 = px.line(
                    df_plt, x=filter_col, y=cpd_names, width=1400, height=1000
                )
                st.plotly_chart(fig2, theme="streamlit", use_container_width=True)
                st.session_state["df_cpds"] = df_cpd
                profile_conn.close()
                conn.close()
            # st.session_state["df_all_data"] = df_hist
