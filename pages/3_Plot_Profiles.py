import pandas as pd
import plotly.express as px
import pqdm
import psycopg2
import streamlit as st
from pqdm.processes import pqdm
from sklearn.metrics.pairwise import cosine_similarity

  # @st.cache_data
def find_sim_cpds(df1, df2):
        filter_col1 = [col for col in df1.columns if not col.startswith("meta")]
        filter_col2 = [col for col in df2.columns if not col.startswith("meta")]
        filter_col = list(set(filter_col1) & set(filter_col2))
        simi = cosine_similarity(df1[filter_col], df2[filter_col])
        return simi
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
    cpd_pro = st.session_state["df_profiles"]
    crisper_pro = st.session_state["df_crisper"]
    list_sources = cpd_pro["metasource"].unique().tolist()
    st.write("Data from DB", cpd_pro)
    choix_source = st.selectbox("Select the Source", list_sources)

    # -------------------umap
    sql_umqpemd = f"select * from umapemd where metasource='{choix_source}'"
    batchs = cpd_pro[cpd_pro["metasource"] == choix_source]["metabatchid"].values
    df_src_emd = pd.read_sql(sql_umqpemd, profile_conn)
    df_src_emd["color"] = "others"
    df_src_emd.loc[df_src_emd["metabatchid"].isin(batchs), "color"] = "selected compounds"

    fig = px.scatter(
        df_src_emd,
        x="umap1",
        y="umap2",
        color="color",
        title=f"{choix_source} profiles ",
        hover_data=["metabatchid"],
    )
    st.plotly_chart(fig, theme="streamlit", use_container_width=True)
    
    sql_umqpemd_crips = f"select * from umapemd where metasource='source_13'"
    df_crisper_emd = pd.read_sql(sql_umqpemd_crips, profile_conn)

    # ---- source and crisper profile ----------------------------------------------------------------------------
    sql_profile = f"select * from aggcombatprofile where metasource='{choix_source}'"
    df_source = pd.read_sql(sql_profile, profile_conn)

    sql_crisper_profile = f"SELECT * FROM aggcombatprofile WHERE metasource='source_13'"
    df_prof_crisper = pd.read_sql(sql_crisper_profile, profile_conn)

    # ---------------------------------------------------------------------------

    col3, col4 = st.columns(2)
    with col3:
        choix = st.selectbox("Select the Profile", batchs)

    df_sel = cpd_pro[
        (cpd_pro["metasource"] == choix_source) & (cpd_pro["metabatchid"] == choix)
    ]
    with col4:
        st.write("Selected Profile", df_sel)
    col3, col4 = st.columns(2)
    with col3:
        thres = st.slider("Threshold", -1.0, 1.0, 0.85)
    with col4:
        thresq = st.slider("Cardinal Threshold", 0, 200, 10)
    # ---------------------------------------------------------------------------

    sim_crisper = find_sim_cpds(df_prof_crisper, df_sel)
    df_hist_crisper = pd.DataFrame({"sim": sim_crisper.flatten().tolist(), "metabatchid": df_prof_crisper["metabatchid"]})
    df_keep_crisper = df_hist_crisper[df_hist_crisper ["sim"] > thres].head(thresq).sort_values(by="sim", ascending=False).reset_index(drop=True)
    batch_list_crisper = df_keep_crisper["metabatchid"].tolist()
    b_list_crisper = [f"'{b}'" for b in batch_list_crisper if "jcp2022_800" not in b]
    if len(b_list_crisper) > 0:
        sql_crisper = f"select gene.geneid,gene.symbol,crisperbatchs.batchid  from gene inner join crisperbatchs on crisperbatchs.geneid=gene.geneid where crisperbatchs.batchid in ({','.join(b_list_crisper)})"
        df_results_cripser = pd.read_sql(sql_crisper, conn)
        if len(df_results_cripser) > 0:
            st.session_state["df_crisper"] = df_results_cripser
            df_keep_prof_crisper = df_prof_crisper[df_prof_crisper["metabatchid"].isin(df_keep_crisper["metabatchid"].values)]
            st.session_state["df_crisper_profile"] = df_keep_prof_crisper
            
    sim_cpds = find_sim_cpds(df_source, df_sel)
    df_hist_cpd = pd.DataFrame({"sim": sim_cpds.flatten().tolist(), "metabatchid": df_source["metabatchid"]})
    df_keep_cpd = df_hist_cpd[df_hist_cpd["sim"] > thres].head(thresq).sort_values(by="sim", ascending=False).reset_index(drop=True)
    batch_list_cpd = df_keep_cpd["metabatchid"].tolist()
    b_list_cpd = [f"'{b}'" for b in batch_list_cpd if "jcp2022_800" not in b]
    if len(b_list_cpd) > 0:
        sql_cpds = f"select cpd.pubchemid,cpd.keggid, cpd.name, cpd.smile from cpd inner join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid where cpdbatchs.batchid in ({','.join(b_list_cpd)})"
        df_results_cpd = pd.read_sql(sql_cpds, conn)
        if len(df_results_cpd) > 0:
                st.session_state["df_cpds"] = df_results_cpd
                df_keep_prof = df_source[df_source["metabatchid"].isin(df_keep_cpd["metabatchid"].values)]
                st.session_state["df_cpds_profile"] = df_keep_prof
        
    
                

#--------------------------------------------------------------------
    tab1, tab2 = st.tabs(["Similar profiles in compounds", "Similar profiles in Cripser"])
    with tab1:
        fig1 = px.histogram(df_hist_cpd, x="sim")
        fig1.add_vline(x=thres)
        if len(df_results_cpd) > 0:
                tab_list = st.tabs(["Histogram", "Data", "UMAP", "MetaData", "Summary"])
                with tab_list[0]:
                    st.plotly_chart(fig1, theme="streamlit", use_container_width=True)
                with tab_list[1]:
                    st.write(df_keep_cpd)
                with tab_list[2]:
                    df_src_emd["color"] = "others"
                    df_src_emd.loc[df_src_emd["metabatchid"].isin(batch_list_cpd), "color"] = "similar compounds"
                    df_src_emd.loc[df_src_emd["metabatchid"] == choix, "color"] = "selected compounds"
                    fig = px.scatter(
                        df_src_emd,
                        x="umap1",
                        y="umap2",
                        color="color",
                        title=f"similar cpds to {choix} profiles  ",
                        hover_data=["metabatchid"],
                    )
                    st.plotly_chart(fig, theme="streamlit", use_container_width=True)

                with tab_list[3]:
                    st.write(df_results_cpd)

                with tab_list[4]:
                    st.write(df_keep_cpd.describe())

                
                st.write(df_keep_prof)

                if len(df_keep_prof) < 11:
                    cpd_names = df_keep_prof.metabatchid.values
                    df_plt = df_keep_prof.set_index("metabatchid")
                    filter_col = [col for col in df_plt.columns if not col.startswith("meta")]
                    df_plt = df_plt[filter_col].T
                    fig1 = px.line(df_plt, x=filter_col, y=cpd_names, width=1400, height=1000)
                    st.plotly_chart(fig1, theme="streamlit", use_container_width=True)
    with tab2:
        fig2 = px.histogram(df_hist_crisper, x="sim")
        fig2.add_vline(x=thres)
        if len(df_results_cripser) > 0:
                tab_list = st.tabs(["Histogram", "Data", "UMAP", "MetaData", "Summary"])
                with tab_list[0]:
                    st.plotly_chart(fig2, theme="streamlit", use_container_width=True)
                with tab_list[1]:
                    st.write(df_keep_crisper)
                with tab_list[2]:
                    df_crisper_emd["color"] = "others"
                    df_crisper_emd.loc[df_crisper_emd["metabatchid"].isin(batch_list_crisper), "color"] = "similar profile"
                    fig = px.scatter(
                        df_crisper_emd,
                        x="umap1",
                        y="umap2",
                        color="color",
                        title=f"similar cpds to {choix} CRISPER profiles  ",
                        hover_data=["metabatchid"],
                    )
                    st.plotly_chart(fig, theme="streamlit", use_container_width=True)

                with tab_list[3]:
                    st.write(df_results_cripser)

                with tab_list[4]:
                    st.write(df_keep_crisper.describe())

                
                st.write(df_keep_crisper)

                if len(df_keep_crisper) < 11:
                    cpd_names = df_keep_crisper.metabatchid.values
                    df_plt = df_keep_crisper.set_index("metabatchid")
                    filter_col = [col for col in df_plt.columns if not col.startswith("meta")]
                    df_plt = df_plt[filter_col].T
                    fig1 = px.line(df_plt, x=filter_col, y=cpd_names, width=1400, height=1000)
                    st.plotly_chart(fig1, theme="streamlit", use_container_width=True)


