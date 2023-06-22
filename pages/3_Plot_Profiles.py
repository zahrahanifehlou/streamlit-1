import pandas as pd
import plotly.express as px
import pqdm
import psycopg2
import streamlit as st
from pqdm.processes import pqdm
from sklearn.metrics.pairwise import cosine_similarity

if "df_profiles" not in st.session_state:
    st.write("Connect DB First")
else:
    all_df = st.session_state["df_profiles"]
    st.write("Data from DB", all_df)
    list_cpds = all_df["metasource"] + "+" + all_df["metabatchid"]
    col1, col2 = st.columns(2)
    with col1:
        choix = st.selectbox("Select the Profile", list_cpds)
        list_sources = all_df["metasource"].unique()
        choix_source = st.selectbox("Select the Source", list_sources)
        
    source = choix.split("+")[0]
    bid = choix.split("+")[1]
    df_sel = all_df[(all_df["metasource"] == source) & (all_df["metabatchid"] == bid)]
    with col2:
        st.write("Selected Profile", df_sel)

    conn2 = psycopg2.connect(
        host="192.168.2.131",
        port="5432",
        user="arno",
        database="ksilink_cpds",
        password="12345",
    )
    sql_profile = (
        "select * from aggcombatprofile where metasource=" + "'" + choix_source + "'"
    )

    df_source = pd.read_sql(sql_profile, conn2)

    list_cpd = []
    progress_text = "Computing similarities. Please wait"
    list_keep = []
    list_cpd_keep = []

    # @st.cache_data
    def crisper_find_sim_from_cpd2(df1, df2):
        simi = cosine_similarity(df1[filter_col], df2[filter_col])
        return simi

 
    filter_col1 = [col for col in df_source.columns if not col.startswith("meta")]
    filter_col2 = [col for col in df_sel.columns if not col.startswith("meta")]
    filter_col = list(set(filter_col1) & set(filter_col2))
    sim_test = crisper_find_sim_from_cpd2(df_source, df_sel)
  
    df_hist = pd.DataFrame()
    df_hist["sim"] = sim_test.flatten().tolist()
    df_hist["metabatchid"] = df_source["metabatchid"]
 
    thres = st.slider("Threshold", -1.0, 1.0, 0.75)
    thresq = st.slider("Cardinal Threshold", 0, 200, 50)
    fig = px.histogram(df_hist, x="sim")
    fig.add_vline(x=thres)
    df_keep = df_hist[df_hist["sim"] > thres]
    df_keep = df_keep.head(thresq).reset_index()

    tab1, tab2, tab3 = st.tabs(["Histogram", "Data", "Summary"])
    with tab1:
        st.plotly_chart(fig, theme="streamlit", use_container_width=True)
    with tab2:
        st.write(df_keep)
    with tab3:
        st.write(df_keep.describe())
    # l700_KEGG[l700_KEGG["gene_symbol"].isin(u2os_gene_list)]
    df_keep_prof = df_source[df_source["metabatchid"].isin(df_keep["metabatchid"].values)]
    st.write(df_keep_prof)

    if len(df_keep_prof) < 11:
        cpd_names = df_keep_prof.metabatchid.values
        df_plt = df_keep_prof.set_index("metabatchid")
        df_plt = df_plt[filter_col].T
        fig2 = px.line(df_plt, x=filter_col, y=cpd_names, width=1400, height=1000)
        st.plotly_chart(fig2, theme="streamlit", use_container_width=True)
        # chart = st.line_chart(df_keep_prof,)

    if len(df_keep_prof) > 0:
        st.session_state["df_keep"] = df_keep_prof
    st.session_state["df_all_data"] = df_hist
    # df_hist.to_csv("sim_selected_source2.csv", index=None)
