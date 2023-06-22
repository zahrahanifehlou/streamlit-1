import pandas as pd
import plotly.express as px
import pqdm
import psycopg2
import streamlit as st
from pqdm.processes import pqdm
from sklearn.metrics.pairwise import cosine_similarity

if "df" not in st.session_state:
    st.write("Connect DB First")
else:
    df = st.session_state["df"]
    st.write("Data from DB", df)
    list_cpds = df["metasource"] + "+" + df["metabatchid"]
    choix = st.selectbox("Select the Profile", list_cpds)
    source = choix.split("+")[0]
    bid = choix.split("+")[1]
    df_sel = df[(df["metasource"] == source) & (df["metabatchid"] == bid)]
    st.write("Selected Profile", df_sel)
    list_sources = [
        "source_1",
        "source_2",
        "source_3",
        "source_4",
        "source_5",
        "source_6",
        "source_7_25",
        "source_7_625",
        "source_8",
        "source_9",
        "source_10",
        "source_11",
        "source_13",
    ]
    choix_source = st.selectbox("Select the Source", list_sources)
   
    x_list=[col for col in df.columns if not col.startswith('meta')]
    

    conn2 = psycopg2.connect(
        host="192.168.2.131",
        port="5432",
        user="arno",
        database="ksilink_cpds",
        password="12345",
    )
    sql_profile = "select * from aggcombatprofile where metasource=" + "'" + choix_source + "'"
    df_prof = pd.read_sql(sql_profile, conn2)

    # df_pro = pd.read_feather(f"/mnt/shares/L/DB/updateDB/clean_profiles/{source}_feature_select.fth")
    filter_col = [col for col in df_prof.columns if not col.startswith('meta')]
   

    list_cpd = []
    progress_text = "Computing similarities. Please wait"
    list_keep = []
    list_cpd_keep = []

    # @st.cache_data
    def crisper_find_sim_from_cpd2(df, df_se2l):
        simi = cosine_similarity(df[filter_col], df_se2l[filter_col])
        return simi
    df_prof=df_prof[df_prof.notna()]
    df_sel=df_sel[df_sel.notna()]

    sim_test = crisper_find_sim_from_cpd2(df_prof, df_sel)

    # st.write(sim_test.shape)
    # sim_all = crisper_find_sim_from_cpd(df_prof, "batchid")
    df_hist = pd.DataFrame()
    df_hist["sim"] = sim_test.flatten().tolist()
    df_hist["metabatchid"] = df_prof["metabatchid"]
    # st.write(df_hist.head())
    # thres = st.text_input("Enter your search")
    thres = st.slider("Threshold", -1.0, 1.0, 0.75)
    thresq = st.slider("Cardinal Threshold", 0, 200, 50)
    fig = px.histogram(df_hist, x="sim")
    fig.add_vline(x=thres)
    df_keep = df_hist[df_hist["sim"] > thres]
    df_keep = df_keep.head(thresq).reset_index()
    # df_keep = pd.DataFrame()
    # df_keep["cpds"] = list_cpd_keep
    # df_keep["sim"] = list_keep
    tab1, tab2, tab3 = st.tabs(["Histogram", "Data", "Summary"])
    with tab1:
        st.plotly_chart(fig, theme="streamlit", use_container_width=True)
    with tab2:
        st.write(df_keep)
    with tab3:
        st.write(df_keep.describe())
    # l700_KEGG[l700_KEGG["gene_symbol"].isin(u2os_gene_list)]
    df_keep_prof = df_prof[df_prof["metabatchid"].isin(df_keep["metabatchid"].values)]
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
    df_hist.to_csv("sim_selected_source2.csv", index=None)
