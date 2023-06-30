# NMDA GSK2879552 Bay K
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import pqdm
import psycopg2
import seaborn as sns
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
    df_cpds = st.session_state["df_cpds"]
    cpd_pro = st.session_state["df_profiles"]
    crisper_pro = st.session_state["df_crisper"]
    list_sources = cpd_pro["metasource"].unique().tolist()
    st.write("Data from DB", cpd_pro)
    choix_source = st.selectbox("Select the Source", list_sources)

    # umap--------------------------------------------------------
    sql_umqpemd_crips = f"select * from umapemd where metasource='source_13'"
    df_crisper_emd = pd.read_sql(sql_umqpemd_crips, profile_conn)
    sql_umqpemd = f"select * from umapemd where metasource='{choix_source}'"
    batchs = cpd_pro[cpd_pro["metasource"] == choix_source]["metabatchid"].values
    df_src_emd = pd.read_sql(sql_umqpemd, profile_conn)
    df_src_emd["color"] = "others"
    df_src_emd.loc[df_src_emd["metabatchid"].isin(batchs), "color"] = "selected compounds"

    # clustermap-------------------------------------------------------------
    list_col = [col for col in cpd_pro.columns if not col.startswith("meta")]
    df_plt = pd.DataFrame()
    df_plt = cpd_pro[list_col]
    df_plt["name"] = cpd_pro["metabatchid"] + "_" + cpd_pro["metasource"]
    df_plt.set_index("name", inplace=True)
    fig, ax = plt.subplots()
    fig = sns.clustermap(
        df_plt,
         metric="cosine",
       # method="ward",
        xticklabels=False,
        yticklabels=True,
        col_cluster=False,
        cmap="vlag",
        center=0,
        vmin=-5,
        vmax=5,
    )
    st.pyplot(fig)

    # ---- source and crisper profile ----------------------------------------------------------------------------
    sql_profile = f"select * from aggcombatprofile where metasource='{choix_source}'"
    df_source = pd.read_sql(sql_profile, profile_conn)

    sql_crisper_profile = f"SELECT * FROM aggcombatprofile WHERE metasource='source_13'"
    df_prof_crisper = pd.read_sql(sql_crisper_profile, profile_conn)

    # ---------------------------------------------------------------------------

    col3, col4 = st.columns(2)
    with col3:
        choix = st.selectbox("Select the Profile", batchs)

    df_sel = cpd_pro[(cpd_pro["metasource"] == choix_source) & (cpd_pro["metabatchid"] == choix)]
    with col4:
        st.write("Selected Profile", df_sel)
    col3, col4 = st.columns(2)
    with col3:
        thres = st.slider("Threshold", -1.0, 1.0, 0.85)
    with col4:
        thresq = st.slider("Cardinal Threshold", 0, 200, 10)
        
    # ---------------------------------------------------------------------------
    # crisper
    sim_crisper = find_sim_cpds(df_prof_crisper, df_sel)
    df_hist_crisper = pd.DataFrame(
        {"sim": sim_crisper.flatten().tolist(), "metabatchid": df_prof_crisper["metabatchid"]}
    )
    df_keep_crisper = (
        df_hist_crisper[df_hist_crisper["sim"] > thres]
        .head(thresq)
        .sort_values(by="sim", ascending=False)
        .reset_index(drop=True)
    )
    batch_list_crisper = df_keep_crisper["metabatchid"].tolist()
    b_list_crisper = [f"'{b}'" for b in batch_list_crisper if "jcp2022_800" not in b]
    df_results_cripser = pd.DataFrame()
    if len(b_list_crisper) > 0:
        sql_crisper = f"select gene.geneid,gene.symbol,crisperbatchs.batchid \
        from gene inner join crisperbatchs on crisperbatchs.geneid=gene.geneid \
        where crisperbatchs.batchid in ({','.join(b_list_crisper)})"
        df_results_cripser = pd.read_sql(sql_crisper, conn)
        if len(df_results_cripser) > 0:
            st.session_state["df_crisper"] = df_results_cripser
            df_keep_prof_crisper = df_prof_crisper[
                df_prof_crisper["metabatchid"].isin(df_keep_crisper["metabatchid"].values)
            ]
            df_keep_prof_crisper.reset_index(inplace=True,drop=True)
            dic_gene = df_results_cripser.set_index("batchid").to_dict()["geneid"]
            df_keep_prof_crisper["metageneid"] = df_keep_prof_crisper["metabatchid"].map(dic_gene)
            df_keep_prof_crisper["metatype"] = "CRISPER"   
            df_keep_prof_crisper["metaefficacy"]=None
            st.session_state["df_crisper_profile"] = df_keep_prof_crisper
    
    # cpds
    sim_cpds = find_sim_cpds(df_source, df_sel)
    df_hist_cpd = pd.DataFrame({"sim": sim_cpds.flatten().tolist(), "metabatchid": df_source["metabatchid"]})
    df_keep_cpd = (
        df_hist_cpd[df_hist_cpd["sim"] > thres]
        .head(thresq)
        .sort_values(by="sim", ascending=False)
        .reset_index(drop=True)
    )
    df_keep_cpd.loc[len(df_keep_cpd.index)] = [1, choix]
    batch_list_cpd = df_keep_cpd["metabatchid"].tolist()
    b_list_cpd = [f"'{b}'" for b in batch_list_cpd if "jcp2022_800" not in b]
    df_results_cpd = pd.DataFrame()
    if len(b_list_cpd) > 0:
        sql_cpds = f"select cpd.pubchemid,cpd.keggid, cpd.name, cpd.smile,cpdgene.geneid,cpdbatchs.batchid,keggcpd.efficacy from cpd \
        inner join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid \
        left join cpdgene on cpdbatchs.pubchemid=cpdgene.pubchemid \
        left join keggcpd on cpd.keggid=keggcpd.keggid \
        where cpdbatchs.batchid in ({','.join(b_list_cpd)}) group by cpd.pubchemid,cpd.keggid, cpd.name, cpd.smile,cpdgene.geneid,cpdbatchs.batchid,keggcpd.efficacy"
        df_results_cpd = pd.read_sql(sql_cpds, conn)
        df_results_cpd.drop_duplicates(subset=["pubchemid"], inplace=True)
        if len(df_results_cpd) > 0:
            st.session_state["df_cpds"] = pd.concat([df_results_cpd, df_cpds])
            df_keep_prof_cpd = df_source[df_source["metabatchid"].isin(df_keep_cpd["metabatchid"].values)]
            df_keep_prof_cpd.reset_index(inplace=True,drop=True)
            dic_gene = df_results_cpd.set_index("batchid").to_dict()["geneid"]
            dic_efficacy= df_results_cpd.set_index("batchid").to_dict()["efficacy"]
            df_keep_prof_cpd["metageneid"] = df_keep_prof_cpd["metabatchid"].map(dic_gene)
            df_keep_prof_cpd["metaefficacy"] = df_keep_prof_cpd["metabatchid"].map(dic_efficacy)
            df_keep_prof_cpd["metatype"] = "CPD"  
            st.session_state["df_cpds_profile"] = df_keep_prof_cpd

    # --------------------------------------------------------------------
    tab1, tab2 = st.tabs(["Similar profiles in compounds", "Similar profiles in Cripser"])
    with tab1:
        fig1 = px.histogram(df_hist_cpd, x="sim")
        fig1.add_vline(x=thres)
        if len(df_results_cpd) > 0:
            tab_list = st.tabs(["Histogram", "Data", "UMAP", "MetaData", "Summary","similar cpds profile"])
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
            with tab_list[5]:
                st.write(df_keep_prof_cpd)

            if len(df_keep_prof_cpd) < 11:
                cpd_names = df_keep_prof_cpd.metabatchid.values
                df_plt = df_keep_prof_cpd.set_index("metabatchid")
                filter_col = [col for col in df_plt.columns if not col.startswith("meta")]
                df_plt = df_plt[filter_col].T
                fig1 = px.line(df_plt, x=filter_col, y=cpd_names, width=1400, height=1000)
                st.plotly_chart(fig1, theme="streamlit", use_container_width=True)
    with tab2:
        fig2 = px.histogram(df_hist_crisper, x="sim")
        fig2.add_vline(x=thres)
        if len(df_results_cripser) > 0:
            tab_list = st.tabs(["Histogram", "Data", "UMAP", "MetaData", "Summary","similar crisper profile"])
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
            with tab_list[5]:
                st.write(df_keep_prof_crisper)

#compare CPD and CRISPER---------------------------------------------------------------------------------------------------------------
    st.write("\n")
    st.write("## compare CPD and CRISPER")
    
    import umap
    compare_cols = st.columns(2)
    meta_cols = [col for col in df_keep_prof_cpd.columns if  col.startswith("meta")]
    with compare_cols[0]:
        st.write("compounds profile")
        st.write(df_keep_prof_cpd[meta_cols].describe().T)
    with compare_cols[1]:
        st.write("Crisper profile")
        st.write(df_keep_prof_crisper[meta_cols].describe().T)
   

    def find_umap(df, title):
        filter_cols = [col for col in df.columns if not col.startswith("meta")]
        reducer = umap.UMAP(densmap=True, random_state=42, verbose=True)
        embedding = reducer.fit_transform(df[filter_cols])
        df_emb = pd.DataFrame({"x": embedding[:, 0], "y": embedding[:, 1]})
        df_emb["target"] = df["metageneid"]
        df_emb["Type"] = df["metatype"]
        df_emb["efficacy"] = df["metaefficacy"]
       
        fig = px.scatter(df_emb, x="x", y="y", hover_data=["target","efficacy","Type"], color="Type", title=title)
        st.plotly_chart(fig, theme="streamlit", use_container_width=True)

    cols = st.columns(3)
    with cols[0]:
        find_umap(df_keep_prof_cpd, "UMAP in CPD profile")
    with cols[1]:
        find_umap(df_keep_prof_crisper, "UMAP in Crisper profile")
    with cols[2]:
        df = pd.concat([df_keep_prof_cpd, df_keep_prof_crisper]).reset_index(drop=True)
        find_umap(df, "UMAP in CPD and Crisper profile")

    list_col = [col for col in df.columns if not col.startswith("meta")]
    df_plt = pd.DataFrame()
    df_plt = df[list_col]
    df_plt["name"] = df["metabatchid"] + "_" + df["metatype"]
    df_plt.set_index("name", inplace=True)
    fig, ax = plt.subplots()
    fig = sns.clustermap(
        df_plt,
         metric="cosine",
       # method="ward",
        xticklabels=False,
        yticklabels=True,
        col_cluster=False,
        cmap="vlag",
        center=0,
        vmin=-5,
        vmax=5,
    )
    st.pyplot(fig)
    st.write("## Similarity")
    filter_cols = [col for col in df.columns if not col.startswith("meta")]
    simi = cosine_similarity(df[filter_cols])
    sim_df = pd.DataFrame(simi)
    sim_df['target']=df['metageneid'] +' _ ' + df['metatype']
    sim_cols = st.columns(2)
    with sim_cols[0]:     
        st.write('sim_df',sim_df)
    with sim_cols[1]: 
        fig = px.imshow(sim_df, title="Profiles Similarities")
        st.plotly_chart(fig)
