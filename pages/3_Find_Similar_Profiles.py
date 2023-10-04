# NMDA GSK2879552 Bay K
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
from streamlib import sql_df

st.set_page_config(
    layout="wide",
)
st.header("Find Similar profiles",divider='rainbow')
def convert_df(df):
       return df.to_csv(index=False).encode('utf-8')
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

def find_sim_cpds(df1, df2):
    filter_col1 = [col for col in df1.columns if not col.startswith("meta")]
    filter_col2 = [col for col in df2.columns if not col.startswith("meta")]
    filter_col = list(set(filter_col1) & set(filter_col2))
    simi = cosine_similarity(df1[filter_col], df2[filter_col])
    return simi

def find_umap(df, title):
    filter_cols = [col for col in df.columns if not col.startswith("meta")]
    meta_cols = [col for col in df.columns if  col.startswith("meta")]
    reducer = umap.UMAP(densmap=True, random_state=42, verbose=True)
    embedding = reducer.fit_transform(df[filter_cols])
    df_emb = pd.DataFrame({"x": embedding[:, 0], "y": embedding[:, 1]})
    df_emb[meta_cols] = df[meta_cols]
    fig = px.scatter(df_emb, x="x", y="y", hover_data=meta_cols, color="metasource", title=title)
    st.plotly_chart(fig, theme="streamlit", use_container_width=True)

def get_col_colors(df):
    list_col = [col for col in df.columns if not col.startswith("Meta")]
    ER = [
        x
        for x in list_col
        if "ER" in x and all(y not in x for y in ["RNA", "Mito", "AGP", "DNA"])
    ]
    RNA = [
        x
        for x in list_col
        if "RNA" in x and all(y not in x for y in ["ER", "Mito", "AGP", "DNA"])
    ]
    Mito = [
        x
        for x in list_col
        if "Mito" in x and all(y not in x for y in ["ER", "RNA", "AGP", "DNA"])
    ]
    mito = [
        x
        for x in list_col
        if "mito" in x and all(y not in x for y in ["ER", "RNA", "AGP", "DNA"])
    ]
    AGP = [
        x
        for x in list_col
        if "AGP" in x and all(y not in x for y in ["ER", "RNA", "Mito", "DNA"])
    ]
    DNA = [
        x
        for x in list_col
        if "DNA" in x and all(y not in x for y in ["ER", "RNA", "Mito", "AGP"])
    ]
    list_fin = []
    list_fin.extend(DNA)
    list_fin.extend(RNA)
    list_fin.extend(ER)
    list_fin.extend(AGP)
    list_fin.extend(Mito)
    list_fin.extend(mito)

    list_fin = list(dict.fromkeys(list_fin))


    list_fin.append("name")
    df["name"] = df["metacpdname"] + "_" + df["metasource"]

    df_plt = df[list_fin]
    df_plt.set_index("name", inplace=True)

    col_colors = []


    for col in  df_plt.columns:
        if col in ER:
            col_colors.append("red")
        elif col in DNA:
            col_colors.append("blue")
        elif col in RNA:
            col_colors.append("green")
        elif col in AGP:
            col_colors.append("orange")
        elif col in Mito or col in mito:
            col_colors.append("pink")
        else:
            col_colors.append("white")
    return df_plt,col_colors


#-----------------------------------------------------------------------------------------------------------------------------------------
if "df_profiles" not in st.session_state:
    st.write("Connect DB First")
else:
    df_cpds = st.session_state["df_cpds"]
    cpd_pro = st.session_state["df_profiles"]
    list_sources = cpd_pro["metasource"].unique().tolist()


    cols1= st.columns(2)
    with cols1[0]:
        choix_source = st.selectbox ("Select the Source", list_sources)
        batchs = cpd_pro[cpd_pro["metasource"] == choix_source]["metabatchid"].values
        choix = st.selectbox("Select the Profile", batchs)
    df_sel = cpd_pro[(cpd_pro["metasource"] == choix_source) & (cpd_pro["metabatchid"] == choix)]
    with cols1[1]:
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


    ## --------------------------------------------------------------------
    tab1, tab2 = st.tabs(["Similar profiles in compounds", "Similar profiles in Cripser"])

    ## cpds---------------------------------------------------------------------------------------------------------------------
    with tab1: # cpds -------------------
        cols_cpds = st.columns(2)
        with cols_cpds[0]:
            thres_cpd = st.slider("Threshold cpd", -1.0, 1.0, 0.75)
        with cols_cpds[1]:
            thresq_cpd = st.slider("Cardinal Threshold cpd", 0, 1000, 10)
        sim_cpds = find_sim_cpds(df_source, df_sel)


        df_hist_cpd = pd.DataFrame(
            {"sim": sim_cpds.flatten().tolist(), "metabatchid": df_source["metabatchid"]}
        )
        df_keep_cpd = (
            df_hist_cpd[df_hist_cpd["sim"] > thres_cpd]
            .head(thresq_cpd)
            .sort_values(by="sim", ascending=False)
            .reset_index(drop=True)
        )

        batch_list_cpd = df_keep_cpd["metabatchid"].tolist()
        b_list_cpd = [f"'{b}'" for b in batch_list_cpd if "jcp2022_800" not in b]
        df_results_cpd = pd.DataFrame()
        df_keep_prof_cpd=pd.DataFrame()
        if len(b_list_cpd) > 0:

            sql_cpds = f"select cpd.pubchemid,cpd.keggid, cpd.cpdname, cpd.smile,cpdgene.geneid,cpdbatchs.batchid,keggcpd.efficacy from cpd \
            inner join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid \
            left join cpdgene on cpdbatchs.pubchemid=cpdgene.pubchemid \
            left join keggcpd on cpd.keggid=keggcpd.keggid \
            where cpdbatchs.batchid in ({','.join(b_list_cpd)}) group by cpd.pubchemid,cpd.keggid, cpd.cpdname, cpd.smile,cpdgene.geneid,cpdbatchs.batchid,keggcpd.efficacy"
            df_results_cpd = sql_df(sql_cpds, conn)
            df_results_cpd.drop_duplicates(subset=["pubchemid"], inplace=True)
            if len(df_results_cpd) > 0:
                st.session_state["df_cpds"] = pd.concat([df_results_cpd, df_cpds])


                df_keep_prof_cpd = df_source[df_source["metabatchid"].isin(df_keep_cpd["metabatchid"].values)]
                df_keep_prof_cpd.reset_index(inplace=True,drop=True)
                df_keep_prof_cpd = df_keep_prof_cpd.merge(df_results_cpd.add_prefix('meta'),left_on='metabatchid',right_on='metabatchid').reset_index(drop=True)
                df_keep_prof_cpd.loc[df_keep_prof_cpd.metacpdname =="No result", 'metacpdname'] = None
                df_keep_prof_cpd['metacpdname'] = df_keep_prof_cpd['metacpdname'].str[:30]
                df_keep_prof_cpd['metacpdname'] = df_keep_prof_cpd['metacpdname'].fillna(df_keep_prof_cpd['metabatchid'])




        fig_clusmap_cpd = px.histogram(df_hist_cpd, x="sim")
        fig_clusmap_cpd.add_vline(x=thres_cpd)

        tab_list = st.tabs(["Histogram", "Data", "UMAP",  "Summary","similar cpds profile"])
        with tab_list[0]:#Histogram
            st.plotly_chart(fig_clusmap_cpd, theme="streamlit", use_container_width=True)
        if len(df_results_cpd) > 0:
            with tab_list[1]:#Data

                df_keep_cpd = df_keep_cpd.merge(df_results_cpd,left_on='metabatchid',right_on='batchid').reset_index(drop=True)
                df_keep_cpd=df_keep_cpd.drop(["metabatchid"],axis=1)
                st.write(df_keep_cpd)
                st.download_button(
                        label="Save",data=convert_df(df_keep_cpd),file_name=f"{df_keep_cpd.cpdname[0]}.csv",mime='csv',)
            with tab_list[2]:#UMAP
                df_src_emd["color"] = "others"
                df_src_emd.loc[df_src_emd["batchid"].isin(batch_list_cpd), "color"] = "similar compounds"
                df_src_emd.loc[df_src_emd["batchid"] == choix, "color"] = "selected compounds"
                fig = px.scatter(
                    df_src_emd,
                    x="umap1",
                    y="umap2",
                    color="color",
                    color_discrete_sequence=["blue", "red","green" ],
                    title=f"similar cpds to {choix} profiles  ",
                    hover_data=["batchid"],
                )
                st.plotly_chart(fig, theme="streamlit", use_container_width=True)

            with tab_list[3]:#Summary
                st.write(df_keep_cpd.describe())
            with tab_list[4]: #similar cpds profil
                st.write(df_keep_prof_cpd)

            if len(df_keep_prof_cpd) < 11:
                cpd_names = df_keep_prof_cpd.metacpdname.values
                df_plt = df_keep_prof_cpd.set_index("metacpdname")
                filter_col = [col for col in df_plt.columns if not col.startswith("meta")]
                df_plt = df_plt[filter_col].T
                fig_clusmap = px.line(df_plt, x=filter_col, y=cpd_names, width=1400, height=1000)
                st.plotly_chart(fig_clusmap, theme="streamlit", use_container_width=True)

    ## crisper---------------------------------------------------------------------------------------------------------------------
    with tab2: # crisper ----------------------
        cols_crisper = st.columns(2)
        with cols_crisper[0]:
            thres_crisper = st.slider("Threshold crisper ", -1.0, 1.0, 0.85)
        with cols_crisper[1]:
            thresq_crisper = st.slider("Cardinal Threshold crisper", 0, 200, 10)
        sim_crisper = find_sim_cpds(df_prof_crisper, df_sel)
        df_hist_crisper = pd.DataFrame(
            {"sim": sim_crisper.flatten().tolist(), "metabatchid": df_prof_crisper["metabatchid"]}
        )
        df_keep_crisper = (
            df_hist_crisper[df_hist_crisper["sim"] > thres_crisper]
            .head(thresq_crisper)
            .sort_values(by="sim", ascending=False)
            .reset_index(drop=True)
        )
        batch_list_crisper = df_keep_crisper["metabatchid"].tolist()
        b_list_crisper = [f"'{b}'" for b in batch_list_crisper ]
    # b_list_crisper = [f"'{b}'" for b in batch_list_crisper if "jcp2022_800" not in b]
        df_results_cripser = pd.DataFrame()
        df_keep_prof_crisper=pd.DataFrame()
        if len(b_list_crisper) > 0:
            sql_crisper = f"select gene.geneid,gene.symbol,crisperbatchs.batchid \
            from gene inner join crisperbatchs on crisperbatchs.geneid=gene.geneid \
            where crisperbatchs.batchid in ({','.join(b_list_crisper)})"

            df_results_cripser = sql_df(sql_crisper, conn)
            if len(df_results_cripser) > 0:
                st.session_state["df_crisper"] = df_results_cripser
                df_keep_prof_crisper = df_prof_crisper[
                    df_prof_crisper["metabatchid"].isin(df_keep_crisper["metabatchid"].values)
                ]
                df_keep_prof_crisper.reset_index(inplace=True,drop=True)
                dic_gene = df_results_cripser.set_index("batchid").to_dict()["geneid"]
                df_keep_prof_crisper["metageneid"] = df_keep_prof_crisper["metabatchid"].map(dic_gene)
                df_keep_prof_crisper["metacpdname"] = df_keep_prof_crisper["metabatchid"]
                df_keep_prof_crisper["metaefficacy"]=None
                st.session_state["df_crisper_profile"] = df_keep_prof_crisper

        fig_clusmap_crisper = px.histogram(df_hist_crisper, x="sim")
        fig_clusmap_crisper.add_vline(x=thres_crisper)

        tab_list = st.tabs(["Histogram", "Data", "UMAP",  "Summary","similar crisper profile"])
        with tab_list[0]: #Histogram
            st.plotly_chart(fig_clusmap_crisper, theme="streamlit", use_container_width=True)
        if len(df_results_cripser) > 0:
            with tab_list[1]: #Data
                df_keep_crisper = df_keep_crisper.merge(df_results_cripser,left_on='metabatchid',right_on='batchid').reset_index(drop=True)
                df_keep_crisper=df_keep_crisper.drop(["metabatchid"],axis=1)
                st.write(df_keep_crisper)
                st.download_button(
                        label="Save",data=convert_df(df_keep_crisper),file_name=f"sim_crisper.csv",mime='csv',)

            with tab_list[2]: #UMAP
                df_crisper_emd["color"] = "others"
                df_crisper_emd.loc[df_crisper_emd["batchid"].isin(batch_list_crisper), "color"] = "similar profile"
                fig = px.scatter(
                        df_crisper_emd,
                        x="umap1",
                        y="umap2",
                        color="color",
                        color_discrete_sequence=["blue", "red","green" ],
                        title=f"similar cpds to {choix} CRISPER profiles  ",
                        hover_data=["batchid"],
                    )
                st.plotly_chart(fig, theme="streamlit", use_container_width=True)


            with tab_list[3]: #Summary
                    st.write(df_keep_crisper.describe())
            with tab_list[4]:
                    st.write(df_keep_prof_crisper)
        if len(df_keep_prof_crisper) < 11:
                cpd_names = df_keep_prof_crisper.metabatchid.values
                df_plt = df_keep_prof_crisper.set_index("metabatchid")
                filter_col = [col for col in df_plt.columns if not col.startswith("meta")]
                df_plt = df_plt[filter_col].T
                fig_clusmap = px.line(df_plt, x=filter_col, y=cpd_names, width=1400, height=1000)
                st.plotly_chart(fig_clusmap, theme="streamlit", use_container_width=True)

#compare CPD and CRISPER---------------------------------------------------------------------------------------------------------------
    st.write("\n")

    # if len(df_results_cripser) > 0:
    #     st.write("## compare CPD and CRISPER")

    #     compare_cols = st.columns(2)

    #     with compare_cols[0]:
    #         st.write("compounds profile")
    #         meta_cols = [col for col in df_keep_prof_cpd.columns if  col.startswith("meta")]
    #         st.write(df_keep_prof_cpd[meta_cols].describe().T)
    #     with compare_cols[1]:
    #         st.write("Crisper profile")
    #         meta_cols_crs = [col for col in df_keep_prof_crisper.columns if  col.startswith("meta")]
    #         st.write(df_keep_prof_crisper[meta_cols_crs].describe().T)






    tmp=pd.DataFrame()
    if len(df_keep_prof_crisper) > 0 and len(df_keep_prof_cpd) > 0:
        tmp = pd.concat([df_keep_prof_cpd, df_keep_prof_crisper]).reset_index(drop=True)
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

    if len(tmp)>1:
        import matplotlib.pyplot as plt
        import seaborn as sns
        plt_src,col_colors=get_col_colors(tmp)

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
#profile_conn.close()