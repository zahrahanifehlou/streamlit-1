# NMDA GSK2879552 Bay K
from streamlib import sql_df,  find_sim_cpds, convert_df, find_umap, get_col_colors
import sys
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import psycopg2
import seaborn as sns
import streamlit as st
import umap
from sklearn.metrics.pairwise import cosine_similarity

sys.path.append('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/')


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


# -------------------------------------------------------------------------------------------------
st.set_page_config(
    layout="wide",
)
st.header("Find Similar profiles", divider='rainbow')
if "df_profiles" not in st.session_state:
    st.write("Connect DB First")
else:
    df_cpds = st.session_state["df_cpds"]
    cpd_pro = st.session_state["df_profiles"]
    list_sources = cpd_pro["metasource"].unique().tolist()
    mainCols = st.columns(2)
    with mainCols[0]:
        choix_source = st.selectbox("Select the Source", list_sources)
        cpdnames = cpd_pro[cpd_pro["metasource"]
                           == choix_source]["metacpdname"].values
        choix = st.selectbox("Select the Profile", cpdnames)
    df_sel = cpd_pro[(cpd_pro["metasource"] == choix_source)
                     & (cpd_pro["metacpdname"] == choix)].head(1)
    with mainCols[1]:
        st.write("Selected Profile", df_sel)

    # ---- source and crisper profile ----------------------------------------------------------------------------
    sql_profile = f"select * from aggcombatprofile where metasource='{choix_source}'"
    df_source = sql_df(sql_profile, profile_conn)

    sql_crisper_profile = f"SELECT * FROM aggcombatprofile WHERE metasource='CRISPER'"
    df_prof_crisper = sql_df(sql_crisper_profile, profile_conn)

    # umap--------------------------------------------------------
    sql_umqpemd_crips = f"select * from umapemd where source='CRISPER'"
    df_crisper_emd = sql_df(sql_umqpemd_crips, profile_conn)
    sql_umqpemd = f"select * from umapemd where source='{choix_source}'"
    df_src_emd = sql_df(sql_umqpemd, profile_conn)

    # --------------------------------------------------------------------
    tab1, tab2 = st.tabs(["Similar profiles in compounds",
                         "Similar profiles in Cripser"])

    # cpds---------------------------------------------------------------------------------------------------------------------
    with tab1:  # cpds -------------------
        cols_cpds = st.columns(2)
        with cols_cpds[0]:
            thres_cpd = st.slider("Threshold cpd", -1.0, 1.0, 0.75)
        with cols_cpds[1]:
            thresq_cpd = st.slider("Cardinal Threshold cpd", 0, 1000, 10)
        sim_cpds = find_sim_cpds(df_source, df_sel)

        df_hist_cpd = pd.DataFrame(
            {"sim": sim_cpds.flatten().tolist(
            ), "metabatchid": df_source["metabatchid"]}
        )
        df_keep_cpd = (
            df_hist_cpd[df_hist_cpd["sim"] > thres_cpd]

            .sort_values(by="sim", ascending=False)
            .head(thresq_cpd)
            .reset_index(drop=True)
        )

        batch_list_cpd = df_keep_cpd["metabatchid"].tolist()
        b_list_cpd = [
            f"'{b}'" for b in batch_list_cpd if "jcp2022_800" not in b]
        df_results_cpd = pd.DataFrame()
        df_keep_prof_cpd = pd.DataFrame()
        if len(b_list_cpd) > 0:
            sql_cpds = f"select cpd.pubchemid,cpd.keggid, cpd.cpdname, cpd.smile,cpdgene.geneid,cpdbatchs.batchid,keggcpd.efficacy from cpd \
            inner join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid \
            left join cpdgene on cpdbatchs.pubchemid=cpdgene.pubchemid \
            left join keggcpd on cpd.keggid=keggcpd.keggid \
            where cpdbatchs.batchid in ({','.join(b_list_cpd)}) group by cpd.pubchemid,cpd.keggid, cpd.cpdname, cpd.smile,cpdgene.geneid,cpdbatchs.batchid,keggcpd.efficacy"
            df_results_cpd = sql_df(sql_cpds, conn)
            df_results_cpd.drop_duplicates(subset=["pubchemid"], inplace=True)
            if len(df_results_cpd) > 0:
                st.session_state["df_cpds"] = pd.concat(
                    [df_results_cpd, df_cpds])
                df_keep_prof_cpd = df_source[df_source["metabatchid"].isin(
                    df_keep_cpd["metabatchid"].values)]
                df_keep_prof_cpd.reset_index(inplace=True, drop=True)
                df_keep_prof_cpd = df_keep_prof_cpd.merge(df_results_cpd.add_prefix(
                    'meta'), left_on='metabatchid', right_on='metabatchid').reset_index(drop=True)
                df_keep_prof_cpd.loc[df_keep_prof_cpd.metacpdname ==
                                     "No result", 'metacpdname'] = None
                df_keep_prof_cpd['metacpdname'] = df_keep_prof_cpd['metacpdname'].str[:30]
                df_keep_prof_cpd['metacpdname'] = df_keep_prof_cpd['metacpdname'].fillna(
                    df_keep_prof_cpd['metabatchid'])

        fig_clusmap_cpd = px.histogram(df_hist_cpd, x="sim")
        fig_clusmap_cpd.add_vline(x=thres_cpd)

        tab_list = st.tabs(["Histogram", "Similar CPDs with cosine", "Similar CPDs with UMAP",
                           "Summary", "similar cpds profile"])
        with tab_list[0]:  # Histogram
            st.plotly_chart(fig_clusmap_cpd, theme="streamlit",
                            use_container_width=True)
        if len(df_results_cpd) > 0:
            with tab_list[1]:  # Similar CPDs with cosine

                df_keep_cpd = df_keep_cpd.merge(
                    df_results_cpd, left_on='metabatchid', right_on='batchid').reset_index(drop=True)
                df_keep_cpd = df_keep_cpd.drop(["metabatchid"], axis=1)

                fig_cols1 = st.columns(3)
                with fig_cols1[0]:
                    st.write(df_keep_cpd)
                    
                with fig_cols1[1]:
                    fig = px.pie(df_keep_cpd,  names='geneid',
                                 title=' geneid',
                                 )
                    fig.update_traces(textposition='inside',
                                      textinfo='percent+label')
                    st.plotly_chart(fig, theme="streamlit",
                                    use_container_width=True)
                with fig_cols1[2]:
                    fig = px.pie(df_keep_cpd,  names='efficacy',
                                 title=' efficacy',
                                 )
                    fig.update_traces(textposition='inside',
                                      textinfo='percent+label')
                    st.plotly_chart(fig, theme="streamlit",
                                    use_container_width=True)
                st.download_button(
                        label="Save", data=convert_df(df_keep_cpd), file_name=f"{df_keep_cpd.cpdname[0]}.csv", mime='csv',)
                # ----------plot sim cpds in UMAP
                df_src_emd["color"] = "others"
                df_src_emd.loc[df_src_emd["batchid"].isin(
                    batch_list_cpd), "color"] = "similar compounds"
                df_src_emd.loc[df_src_emd["name"] ==
                               choix, "color"] = "selected compounds"
                fig = px.scatter(
                    df_src_emd,
                    x="umap1",
                    y="umap2",
                    color="color",
                    color_discrete_sequence=["blue", "red", "green"],
                    title=f"similar cpds to {choix} profiles  ",
                    hover_data=["batchid"],
                )
                st.plotly_chart(fig, theme="streamlit",
                                use_container_width=True)

            with tab_list[2]:  # Similar CPDs with UMAP
         
                from sklearn.neighbors import NearestNeighbors
                nb_cluster = st.slider(
                    'Number of neighbors', min_value=2, max_value=30, value=2, step=1)
                X=df_src_emd[["umap1","umap2"]].to_numpy()
             
                neigh = NearestNeighbors(n_neighbors=nb_cluster, n_jobs=-1)
                neigh.fit(X)
                points=df_src_emd[df_src_emd["name"]==choix][["umap1","umap2"]]
                distances, indices = neigh.kneighbors(points)
                nearest_neighbor_name = df_src_emd.loc[indices[ 0,1:], 'name']
                similar_df=df_src_emd[df_src_emd["name"].isin(nearest_neighbor_name)]
                fig_cols2 = st.columns(3)
                with fig_cols2[0]:
                    st.write(similar_df)
                    
                with fig_cols2[1]:
                    fig = px.pie(similar_df,  names='symbol',
                                 title=' geneid',
                                 )
                    fig.update_traces(textposition='inside',
                                      textinfo='percent+label')
                    st.plotly_chart(fig, theme="streamlit",
                                    use_container_width=True)
                with fig_cols2[2]:
                    fig = px.pie(similar_df,  names='efficacy',
                                 title=' efficacy',
                                 )
                    fig.update_traces(textposition='inside',
                                      textinfo='percent+label')
                    st.plotly_chart(fig, theme="streamlit",
                                    use_container_width=True)
            
                
           

                df_src_emd["color"] = "others"
                df_src_emd.loc[df_src_emd["name"].isin(
                    nearest_neighbor_name), "color"] = "similar compounds"
                df_src_emd.loc[df_src_emd["name"] ==
                               choix, "color"] = "selected compounds"
                figUMAP = px.scatter(
                    df_src_emd,
                    x="umap1",
                    y="umap2",
                    color="color",
                   
                    title=f"similar cpds to {choix} profiles  ",
                    hover_data=["batchid"],
                )
                for trace in figUMAP.data:
                    if trace.name=='similar profile':
                        trace.marker.opacity=0.9
                        trace.marker.size=15

                st.plotly_chart(figUMAP, theme="streamlit",
                                use_container_width=True)
             
        


                


            with tab_list[3]:  # Summary
                st.write(df_keep_cpd.describe())
            with tab_list[4]:  # similar cpds profil
                st.write(df_keep_prof_cpd)

            if len(df_keep_prof_cpd) < 11:
                cpd_names = df_keep_prof_cpd.metacpdname.values
                df_plt = df_keep_prof_cpd.set_index("metacpdname")
                filter_col = [
                    col for col in df_plt.columns if not col.startswith("meta")]
                df_plt = df_plt[filter_col].T
                fig_clusmap = px.line(
                    df_plt, x=filter_col, y=cpd_names, width=1400, height=1000)
                st.plotly_chart(fig_clusmap, theme="streamlit",
                                use_container_width=True)

    # crisper---------------------------------------------------------------------------------------------------------------------
    with tab2:  # crisper ----------------------
        cols_crisper = st.columns(2)
        with cols_crisper[0]:
            thres_crisper = st.slider("Threshold crisper ", -1.0, 1.0, 0.85)
        with cols_crisper[1]:
            thresq_crisper = st.slider(
                "Cardinal Threshold crisper", 0, 200, 10)
        sim_crisper = find_sim_cpds(df_prof_crisper, df_sel)
        df_hist_crisper = pd.DataFrame(
            {"sim": sim_crisper.flatten().tolist(
            ), "metabatchid": df_prof_crisper["metabatchid"]}
        )
        df_keep_crisper = (
            df_hist_crisper[df_hist_crisper["sim"] > thres_crisper]
            .sort_values(by="sim", ascending=False)
            .head(thresq_crisper)
            .reset_index(drop=True)
        )
        batch_list_crisper = df_keep_crisper["metabatchid"].tolist()
        b_list_crisper = [f"'{b}'" for b in batch_list_crisper]
    # b_list_crisper = [f"'{b}'" for b in batch_list_crisper if "jcp2022_800" not in b]
        df_results_cripser = pd.DataFrame()
        df_keep_prof_crisper = pd.DataFrame()
        if len(b_list_crisper) > 0:
            sql_crisper = f"select gene.geneid,gene.symbol,crisperbatchs.batchid \
            from gene inner join crisperbatchs on crisperbatchs.geneid=gene.geneid \
            where crisperbatchs.batchid in ({','.join(b_list_crisper)})"

            df_results_cripser = sql_df(sql_crisper, conn)
            if len(df_results_cripser) > 0:
                st.session_state["df_crisper"] = df_results_cripser
                df_keep_prof_crisper = df_prof_crisper[
                    df_prof_crisper["metabatchid"].isin(
                        df_keep_crisper["metabatchid"].values)
                ]
                df_keep_prof_crisper.reset_index(inplace=True, drop=True)
                dic_gene = df_results_cripser.set_index(
                    "batchid").to_dict()["geneid"]
                df_keep_prof_crisper["metageneid"] = df_keep_prof_crisper["metabatchid"].map(
                    dic_gene)
                df_keep_prof_crisper["metacpdname"] = df_keep_prof_crisper["metabatchid"]
                df_keep_prof_crisper["metaefficacy"] = None
                st.session_state["df_crisper_profile"] = df_keep_prof_crisper

        fig_clusmap_crisper = px.histogram(df_hist_crisper, x="sim")
        fig_clusmap_crisper.add_vline(x=thres_crisper)

        tab_list = st.tabs(["Histogram", "Data", "UMAP",
                           "Summary", "similar crisper profile"])
        with tab_list[0]:  # Histogram
            st.plotly_chart(fig_clusmap_crisper,
                            theme="streamlit", use_container_width=True)
        if len(df_results_cripser) > 0:
            with tab_list[1]:  # Data
                df_keep_crisper = df_keep_crisper.merge(
                    df_results_cripser, left_on='metabatchid', right_on='batchid').reset_index(drop=True)
                df_keep_crisper = df_keep_crisper.drop(["metabatchid"], axis=1)
                st.write(df_keep_crisper)
                st.download_button(
                    label="Save", data=convert_df(df_keep_crisper), file_name=f"sim_crisper.csv", mime='csv',)

            with tab_list[2]:  # UMAP
                df_crisper_emd["color"] = "others"
                df_crisper_emd.loc[df_crisper_emd["batchid"].isin(
                    batch_list_crisper), "color"] = "similar profile"
                fig = px.scatter(
                    df_crisper_emd,
                    x="umap1",
                    y="umap2",
                    color="color",
                    color_discrete_sequence=["blue", "red", "green"],
                    title=f"similar cpds to {choix} CRISPER profiles  ",
                    hover_data=["batchid"],
                )
                st.plotly_chart(fig, theme="streamlit",
                                use_container_width=True)

            with tab_list[3]:  # Summary
                st.write(df_keep_crisper.describe())
            with tab_list[4]:
                st.write(df_keep_prof_crisper)
            if len(df_keep_prof_crisper) < 11:
                cpd_names = df_keep_prof_crisper.metabatchid.values
                df_plt = df_keep_prof_crisper.set_index("metabatchid")
                filter_col = [
                    col for col in df_plt.columns if not col.startswith("meta")]
                df_plt = df_plt[filter_col].T
                fig_clusmap = px.line(
                    df_plt, x=filter_col, y=cpd_names, width=1400, height=1000)
                st.plotly_chart(fig_clusmap, theme="streamlit",
                                use_container_width=True)

# compare CPD and CRISPER---------------------------------------------------------------------------------------------------------------
    st.write("\n")

    tmp = pd.DataFrame()
    if len(df_keep_prof_crisper) > 0 and len(df_keep_prof_cpd) > 0:
        tmp = pd.concat([df_keep_prof_cpd, df_keep_prof_crisper]
                        ).reset_index(drop=True)
        cols = st.columns(3)
        with cols[0]:
            find_umap(df_keep_prof_cpd, "UMAP in CPD profile")
        with cols[1]:
            find_umap(df_keep_prof_crisper, "UMAP in Crisper profile")
        with cols[2]:

            find_umap(tmp, "UMAP in CPD and Crisper profile")

    elif len(df_keep_prof_cpd) > 0:
        st.session_state["df_cpds_profile"] = df_keep_prof_cpd
        tmp = df_keep_prof_cpd.copy()
    elif len(df_keep_prof_crisper) > 0:
        tmp = df_keep_prof_crisper.copy()

    if len(tmp) > 1:
        import matplotlib.pyplot as plt
        import seaborn as sns
        plt_src, col_colors = get_col_colors(tmp)

        fig_clusmap, ax1 = plt.subplots()
        fig_clusmap = sns.clustermap(
            plt_src,
            metric="cosine",
            col_colors=col_colors,
            # method="ward",
            xticklabels=False,
            yticklabels=True,
            col_cluster=False,
            cmap="vlag",
            center=0,
            vmin=-5,
            vmax=5,
            figsize=(16, len(plt_src)/2),
        )

        st.pyplot(fig_clusmap)


conn.close()
# profile_conn.close()
