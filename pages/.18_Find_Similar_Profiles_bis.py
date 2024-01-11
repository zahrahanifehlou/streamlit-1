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
    title = f' compounds to {choix} in  {choix_source}'
    with mainCols[1]:
        st.write("Selected Profile", df_sel)

    # ---- source and crisper profile ----------------------------------------------------------------------------
    sql_profile = f"select * from aggprofile where metasource='{choix_source}'"
    df_source = sql_df(sql_profile, profile_conn)

    # sql_crisper_profile = f"SELECT * FROM aggprofile WHERE metasource='CRISPER' or metasource='ORF-Broad'"
    # df_prof_crisper = sql_df(sql_crisper_profile, profile_conn)

    # umap--------------------------------------------------------
    @st.cache_data
    def get_umap(choix_source):
        sql_umqpemd = f"select * from umapemd where metasource='{choix_source}'"
        df_src_emd = sql_df(sql_umqpemd, profile_conn)
        return df_src_emd
    
    
    # 


    cosine_sim_tab, knn_sim_tab = st.tabs(
    ["Similar CPDs with cosine", "Similar CPDs with KNN in UMAP  "])
    with cosine_sim_tab: 
        

            
        rad = st.radio('Pos or Neg',['Correl', 'AntiCorrel'])# Similar CPDs with cosine
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
        if rad=='Correl':
            df_keep_cpd = (
                df_hist_cpd[df_hist_cpd["sim"] > thres_cpd]
                .sort_values(by="sim", ascending=False)
                .head(thresq_cpd)
                .reset_index(drop=True)
            )
        else:
            df_keep_cpd = (
                df_hist_cpd[df_hist_cpd["sim"] <= thres_cpd]
                .sort_values(by="sim", ascending=True)
                .head(thresq_cpd)
                .reset_index(drop=True)
            )
    fig_clusmap_cpd = px.histogram(df_hist_cpd, x="sim")
    fig_clusmap_cpd.add_vline(x=thres_cpd)
    st.plotly_chart(fig_clusmap_cpd, theme="streamlit",
                        use_container_width=True) 
    
    st.write('df_keep',df_keep_cpd)
    
    batch_list_cpd = df_keep_cpd["metabatchid"].tolist()
    b_list_cpd = [f"'{b}'" for b in batch_list_cpd]
    df_results_cpd = pd.DataFrame()
    df_keep_prof_cpd = pd.DataFrame()
    
    if len(b_list_cpd) > 0 and ('ORF' in choix_source or 'CRIS' in choix_source) :
        sql_crisper2 = f"SELECT * FROM genebatchs  WHERE genebatchs.batchid IN ({','.join(b_list_cpd)})"
        df_results_cpd = sql_df(sql_crisper2, conn)
        st.write('df_results_cpd',df_results_cpd)
        b_list_gene = [f"'{b}'" for b in df_results_cpd.geneid]
        sql_kegg=f"SELECT * FROM gene WHERE gene.geneid IN ({','.join(b_list_gene)})"
        df_genes=sql_df(sql_kegg,conn)
        
        df_merge = df_results_cpd.merge(
                df_genes, left_on='geneid', right_on='geneid').reset_index(drop=True)
        # df_merge=df_merge.add_prefix('meta_')
        # list_meta_cols=df_merge.columns
        
        st.write("\n")  # ----------plot sim cpds in UMAP
        df_src_emd = get_umap(choix_source=choix_source)
        st.write('df_src_emd',df_src_emd)
        df_src_emd["color"] = "others"
        df_src_emd.loc[df_src_emd["metageneid"].isin(
            df_results_cpd.geneid), "color"] = "similar compounds"
        choix_batch=str(df_sel['metageneid'].values[0])
        df_src_emd.loc[df_src_emd["metageneid"] ==
                        choix_batch, "color"] = "selected compounds"
        
        # st.write(choix_batch)
        # st.write('choix',df_src_emd[df_src_emd["metabatchid"]==
        #                 df_sel['metabatchid'].values[0]])
        fig = px.scatter(
        df_src_emd,
        x="umap1",
        y="umap2",
        color="color",
        opacity=0.5,
        color_discrete_sequence=["blue", "red", "green"],
        title=f"{rad} {title}:UMAP ",
        hover_data=['metabatchid'],
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

        