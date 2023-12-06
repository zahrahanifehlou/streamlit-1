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
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics.pairwise import cosine_similarity
from streamlit_plotly_events import plotly_events
from PIL import Image
sys.path.append('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/')
import os


def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])


conn = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"
profile_conn = "postgres://arno:12345@192.168.2.131:5432/ksilink_cpds"


# -------------------------------------------------------------------------------------------------
st.set_page_config(
    layout="wide",
)
st.header("Find Similar profiles", divider='rainbow')
if "df_profiles" not in st.session_state:
    st.write("Connect DB First")
else:
   
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
    title = f' similar compounds to {choix} in  {choix_source}'
    with mainCols[1]:
        st.write("Selected Profile", df_sel)

    # ---- source and crisper profile ----------------------------------------------------------------------------
    sql_profile = f"select * from aggcombatprofile where metasource='{choix_source}'"
    df_source = sql_df(sql_profile, profile_conn)

    sql_crisper_profile = f"SELECT * FROM aggcombatprofile WHERE metasource='CRISPER'"
    df_prof_crisper = sql_df(sql_crisper_profile, profile_conn)

    # umap--------------------------------------------------------
    @st.cache_data
    def get_umap(choix_source):
        sql_umqpemd = f"select * from umapemd where metasource='{choix_source}'"
        df_src_emd = sql_df(sql_umqpemd, profile_conn)
        return df_src_emd
    df_src_emd = get_umap(choix_source=choix_source)



    # cpds---------------------------------------------------------------------------------------------------------------------

    cosine_sim_tab, knn_sim_tab = st.tabs(
        ["Similar CPDs with cosine", "Similar CPDs with KNN in UMAP  "])
    with cosine_sim_tab:  # Similar CPDs with cosine
        thr_cols = st.columns(2)
        with thr_cols[0]:
            thres_cpd = st.slider("Threshold cpd", -1.0, 1.0, 0.85)
        with thr_cols[1]:
            thresq_cpd = st.slider("Cardinal Threshold cpd", 0, 1000, 15)
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
            sql_cpds = f"select cpd.pubchemid,cpd.keggid, cpd.cpdname, gene.symbol, cpd.smile,cpdgene.geneid,cpdbatchs.batchid,keggcpd.efficacy from cpd \
            inner join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid \
            left join cpdgene on cpdbatchs.pubchemid=cpdgene.pubchemid \
            left join gene on gene.geneid=cpdgene.geneid \
            left join keggcpd on cpd.keggid=keggcpd.keggid \
            where cpdbatchs.batchid in ({','.join(b_list_cpd)}) group by cpd.pubchemid,cpd.keggid, cpd.cpdname, gene.symbol,cpd.smile,cpdgene.geneid,cpdbatchs.batchid,keggcpd.efficacy"
            df_results_cpd = sql_df(sql_cpds, conn)
            df_results_cpd.drop_duplicates(subset=["pubchemid"], inplace=True)
            df_keep_prof_cpd = df_source[df_source["metabatchid"].isin(
                df_keep_cpd["metabatchid"].values)]
            df_keep_prof_cpd.reset_index(inplace=True, drop=True)
            
           
            if len(df_results_cpd) > 0:
                df_keep_prof_cpd = df_keep_prof_cpd.merge(df_results_cpd.add_prefix(
                    'meta'), left_on='metabatchid', right_on='metabatchid').reset_index(drop=True)
                df_keep_prof_cpd.loc[df_keep_prof_cpd.metacpdname ==
                                     "No result", 'metacpdname'] = None
                df_keep_prof_cpd['metacpdname'] = df_keep_prof_cpd['metacpdname'].str[:30]
                df_keep_prof_cpd['metacpdname'] = df_keep_prof_cpd['metacpdname'].fillna(
                    df_keep_prof_cpd['metabatchid'])

        fig_clusmap_cpd = px.histogram(df_hist_cpd, x="sim")
        fig_clusmap_cpd.add_vline(x=thres_cpd)

        st.plotly_chart(fig_clusmap_cpd, theme="streamlit",
                        use_container_width=True)
        st.session_state["df_sim_cpd"] = df_results_cpd
        if choix_source == "CRISPER":
            sql_crisper2 = f"SELECT gene.symbol, gene.geneid,crisperbatchs.batchid  FROM crisperbatchs  inner join gene \
            on gene.geneid=crisperbatchs.geneid  where crisperbatchs.batchid in ({','.join(b_list_cpd)})  group by gene.symbol, gene.geneid,crisperbatchs.batchid "
            df_keep_cpd = sql_df(
                sql_crisper2, conn).drop_duplicates(subset=["batchid"])
        
            df_prof_crisper = df_prof_crisper.merge(df_keep_cpd.add_prefix(
                'meta'), left_on='metabatchid', right_on='metabatchid').reset_index(drop=True)
      
            meta_cols = [
                col for col in df_keep_prof_cpd.columns if col.startswith("meta")]
            df_keep_cpd[meta_cols] = df_keep_prof_cpd[meta_cols]
    
            df_keep_cpd["metacpdname"] = df_keep_cpd["batchid"]
            df_keep_prof_cpd["metacpdname"] = df_keep_prof_cpd["metabatchid"]
            df_keep_cpd["efficacy"] = None
            st.session_state["df_sim_crisper"] = df_keep_cpd

        if len(df_keep_cpd) > 0:
            st.write("\n")
            if choix_source != "CRISPER":
                df_keep_cpd = df_keep_cpd.merge(
                    df_results_cpd, left_on='metabatchid', right_on='batchid').reset_index(drop=True)
                df_keep_cpd = df_keep_cpd.drop(["metabatchid"], axis=1)
                df_keep_cpd["source"] = choix_source

            fig_cols1 = st.columns(2)
            name = choix_source+choix
            with fig_cols1[0]:
                st.write(f" {title} :  metadata")
                st.write(df_keep_cpd)
            
            with fig_cols1[1]:  # Profile
                st.write(f"  Profile")
                st.write(df_keep_prof_cpd.head(10))
               
        
            
         
            st.write("\n")
            fig_cols2 = st.columns(2)
            with fig_cols2[0]:
                fig = px.pie(df_keep_cpd,  names='symbol',
                             title=f'   {title} : symbol',
                             )
                fig.update_traces(textposition='inside',
                                  textinfo='percent+label')
                st.plotly_chart(fig, theme="streamlit",
                                use_container_width=True)
            with fig_cols2[1]:
                fig = px.pie(df_keep_cpd,  names='efficacy',
                             title=f'{title} : efficacy',
                             )
                fig.update_traces(textposition='inside',
                                  textinfo='percent+label')
                st.plotly_chart(fig, theme="streamlit",
                                use_container_width=True)

            st.write("\n")  # ----------plot sim cpds in UMAP
            df_src_emd["color"] = "others"
            df_src_emd.loc[df_src_emd["metabatchid"].isin(
                batch_list_cpd), "color"] = "similar compounds"
            df_src_emd.loc[df_src_emd["metacpdname"] ==
                           choix, "color"] = "selected compounds"
           
            fig = px.scatter(
                df_src_emd,
                x="umap1",
                y="umap2",
                color="color",
                color_discrete_sequence=["blue", "red", "green"],
                title=f"{title} : UMAP ",
                hover_data=["metabatchid", "metaefficacy",
                            "metageneid", "metakeggid"],
            )
            if choix_source not in(["Ksilink_625","Ksilink_25","CRISPER"]):
                st.plotly_chart(fig, theme="streamlit",
                                use_container_width=True)
            
            else:
                selected_points = plotly_events(fig,click_event=True)
               
                if selected_points:
                    x=selected_points[0]['x']
                    y=selected_points[0]['y']
                
                    tmp = df_src_emd.loc[(df_src_emd['umap1']==x) & (df_src_emd['umap2'] == y) ]
                    batch=tmp.metabatchid.values[0]
                    name=tmp.metacpdname.values[0]
             
                    sql_point = f"select * from platemap where batchid='{batch}' and assay='{choix_source}'"
                    df_plates= sql_df(sql_point, conn)
                  
                    plt_len=len(df_plates)
                    if plt_len>0:
                        br_cols = st.columns(plt_len)
                    
                    for i in range(len(df_plates)):
                    
                        plate=df_plates.plate[i]
                        well=df_plates.well[i]
                        fpath=f"/mnt/shares/L/PROJECTS/JUMP-CP/Checkout_Results/BirdView/{plate}/{plate}_{well}.jpg"
                        if choix_source=="CRISPER":
                            fpath=f"/mnt/shares/L/PROJECTS/JUMP-CRISPR/Checkout_Results/BirdView/{plate}/{plate}_{well}.jpg"
                        if os.path.isfile(fpath):
                            image = Image.open(fpath)
                            with br_cols[i]:
                                st.image(image, caption=f"{name} : {plate} {well}", width =256)

                
    
    
           
      

            st.write("\n")  # ----------plot PROFILE
            tmp = df_keep_prof_cpd.head(15)
            cpd_names = tmp.metacpdname.values
            df_plt = tmp.set_index("metacpdname")
            filter_col = [
                col for col in df_plt.columns if not col.startswith("meta")]
            df_plt = df_plt[filter_col].T
            fig_clusmap = px.line(
                df_plt, x=filter_col, y=cpd_names, width=1400, height=1000, title=f"{title} : 15 similar compounds")
            st.plotly_chart(fig_clusmap, theme="streamlit",
                            use_container_width=True)

            st.write("\n")  # ----------plot PROFILE heatmap
            tmp = df_keep_prof_cpd.copy()
            if len(tmp) > 1:
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

                    figsize=(16, len(plt_src)/2),
                )

                st.pyplot(fig_clusmap)

    with knn_sim_tab:  # Similar CPDs with UMAP
        from sklearn.neighbors import NearestNeighbors
        nb_cluster = st.slider(
            'Number of neighbors', min_value=2, max_value=30, value=10, step=1)

        X = df_src_emd[["umap1", "umap2"]].to_numpy()
        neigh = NearestNeighbors(n_neighbors=nb_cluster, n_jobs=-1)
        neigh.fit(X)
        if choix in df_src_emd["metacpdname"].values:

            points = df_src_emd[df_src_emd["metacpdname"]
                                == choix][["umap1", "umap2"]]

            distances, indices = neigh.kneighbors(points)
            knn_sim_df=pd.concat([df_src_emd[df_src_emd["metacpdname"]==choix],df_src_emd.loc[indices[0,1:]]])
            knn_sim_df.reset_index(inplace=True, drop=True)

            # knn_sim_df.drop(columns=["metasource"],inplace=True)
            df_keep_prof_cpd_knn =df_source[df_source['metabatchid'].isin(knn_sim_df['metabatchid'].values)]  
        
            df_keep_prof_cpd_knn.reset_index(inplace=True, drop=True)
            
            

            df_keep_prof_cpd_knn = df_keep_prof_cpd_knn.merge(knn_sim_df[["metabatchid", "metaefficacy", "metageneid", "metagenesymbol",
                                                              "metasmile", "metakeggid", "metacpdname"]], left_on='metabatchid', right_on='metabatchid').reset_index(drop=True)
            df_keep_prof_cpd_knn.loc[df_keep_prof_cpd_knn.metacpdname ==
                                     "No result", 'metacpdname'] = None
            df_keep_prof_cpd_knn['metacpdname'] = df_keep_prof_cpd_knn['metacpdname'].str[:30]
            df_keep_prof_cpd_knn['metacpdname'] = df_keep_prof_cpd_knn['metacpdname'].fillna(
                df_keep_prof_cpd_knn['metabatchid'])

            if len(knn_sim_df) > 0:
                st.write("\n")
                fig_cols3 = st.columns(2)
                name = choix_source+choix
                with fig_cols3[0]:
                    st.write(f" {title} MetaData NN")
                    st.write(knn_sim_df)
                
                with fig_cols3[1]:  # Profile
                    st.write(f" {title} Profile NN")
                    st.write(df_keep_prof_cpd_knn.head(10))
       
                    

                st.write("\n")
                fig_cols4 = st.columns(2)
                with fig_cols4[0]:
                    fig = px.pie(knn_sim_df,  names='metagenesymbol',
                                 title=f' {title}  : gene symbol',
                                 )
                    fig.update_traces(textposition='inside',
                                      textinfo='percent+label')
                    st.plotly_chart(fig, theme="streamlit",
                                    use_container_width=True)
                with fig_cols4[1]:
                    fig = px.pie(knn_sim_df,  names='metaefficacy',
                                 title=f' {title} : efficacy',
                                 )
                    fig.update_traces(textposition='inside',
                                      textinfo='percent+label')
                    st.plotly_chart(fig, theme="streamlit",
                                    use_container_width=True)

                # -------------------------------------------
                st.write("\n")
                df_src_emd["color"] = "others"
                df_src_emd.loc[indices[0,1:], "color"] = "similar compounds"
                # df_src_emd.loc[df_src_emd["metacpdname"].isin(
                #     nearest_neighbor_name), "color"] = "similar compounds"
                df_src_emd.loc[df_src_emd["metacpdname"] ==
                               choix, "color"] = "selected compounds"
                figUMAP_knn = px.scatter(
                    df_src_emd,
                    x="umap1",
                    y="umap2",
                    color="color",
                     color_discrete_sequence=["blue", "red", "green"],

                    title=f" {title}  :UMAP ",
                    hover_data=["metabatchid", "metaefficacy",
                                "metageneid", "metakeggid"]
                )
                if choix_source not in(["Ksilink_625","Ksilink_25","CRISPER"]):
                    st.plotly_chart(figUMAP_knn, theme="streamlit",
                                    use_container_width=True)

                else:
                    selected_point_knn = plotly_events(figUMAP_knn,click_event=True)
                    if selected_point_knn:
                        x=selected_point_knn[0]['x']
                        y=selected_point_knn[0]['y']
                    
                        tmp = df_src_emd.loc[(df_src_emd['umap1']==x) & (df_src_emd['umap2'] == y) ]
                        batch_knn=tmp.metabatchid.values[0]
                        name_knn=tmp.metacpdname.values[0]
                    
                       
                        sql_point = f"select * from platemap where batchid='{batch_knn}' and assay='{choix_source}'"
                        df_plates= sql_df(sql_point, conn)                        
                        plt_len=len(df_plates)
                        if plt_len>0:
                            br_cols = st.columns(plt_len)
                        for i in range(len(df_plates)):
                            plate=df_plates.plate[i]
                            well=df_plates.well[i]
                            fpath=f"/mnt/shares/L/PROJECTS/JUMP-CP/Checkout_Results/BirdView/{plate}/{plate}_{well}.jpg"
                            if choix_source=="CRISPER":
                                 fpath=f"/mnt/shares/L/PROJECTS/JUMP-CRISPR/Checkout_Results/BirdView/{plate}/{plate}_{well}.jpg"
                            if os.path.isfile(fpath):
                                image = Image.open(fpath)
                                with br_cols[i]:
                                    st.image(image, caption=f"{name_knn} : {plate} {well}", width =256)

    
                    

                #
                st.write("\n")  # ----------plot PROFILE
                tmp = df_keep_prof_cpd_knn.head(15)
                cpd_names = tmp.metacpdname.values
                df_plt = tmp.set_index("metacpdname")
                filter_col = [
                    col for col in df_plt.columns if not col.startswith("meta")]
                df_plt = df_plt[filter_col].T
                fig_clusmap = px.line(
                    df_plt, x=filter_col, y=cpd_names, width=1400, height=1000, title=f"{title} : Profile")
                st.plotly_chart(fig_clusmap, theme="streamlit",
                                use_container_width=True)

                st.write("\n")  # ----------plot PROFILE heatmap
                tmp = df_keep_prof_cpd_knn.copy()
                if len(tmp) > 1:
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


