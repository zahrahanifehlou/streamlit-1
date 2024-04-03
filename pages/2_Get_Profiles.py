import warnings
import numpy as np
import pandas as pd
import psycopg2
import streamlit as st
import plotly.express as px
import sys

sys.path.append("/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit_ZH2/lib/")
# import umap
from streamlib import sql_df, get_sql_jump, get_col_colors

warnings.filterwarnings("ignore")


def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])


conn = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"


table_mapping = {
    "KEGG": {
        "gene": ("keggcpdgene", "geneid"),
        "pathway": ("keggcpdpath", "pathid"),
        "disease": ("keggcpddis", "disid"),
        "cpd": ("keggcpd", "keggid"),
    },
    "SelectChem and Jump DataSet": {
        "gene": ("cpdgene", "geneid"),
        "pathway": ("cpdpath", "pathid"),
        "disease": ("keggcpddis", "disid"),
        "cpd": ("cpd", "pubchemid"),
        "cpdbatchs": ("cpdbatchs", "batchid"),
    },
}

st.set_page_config(
    layout="wide",
)
# -------------------------------------------------------------------------
st.header("Search for Jump profiles", divider="rainbow")
list_tables = ["gene", "cpd", "pathway", "disease", "cpdbatchs"]
col1, col2 = st.columns(2)
with col1:
    option = st.selectbox("Pick one table", list_tables)
    x = st.slider("Number of rows to display", 1, 20)
    st.write("You selected the table: ", option)

data = sql_df("select * from " + str(option), conn)
data = data.apply(lambda x: x.astype(str).str.upper())

with col2:
    tab1, tab2 = st.tabs(["Dataframe", "Summary"])
    tab1.write(data.sample(x))
    tab2.write(data.describe().T)
st.write(f"\n")

# -------------------------------------------------------------------------------------------------------------------
col3, col4 = st.columns(2)
with col3:
    rad_match = st.radio(
        "Exact Match?",
        ["Yes", "No"],
        horizontal=True,
    )
    var_text = st.text_area("Enter your search", help="Name or ID separated by enter")

var_t = var_text.split("\n")
var_t = [t.strip().upper() for t in var_t if t != ""]
var_t = list(set(var_t))
# st.write(var_t)
df_res = pd.DataFrame()
if len(var_t) > 0:
    if rad_match == "Yes":
        mask = np.column_stack(
            [
                data[col].str.fullmatch(v.upper(), na=False)
                for col in data.columns
                for v in var_t
                if v != ""
            ]
        )
    else:
        mask = np.column_stack(
            [
                data[col].str.match(v.upper(), na=False)
                for col in data.columns
                for v in var_t
                if v != ""
            ]
        )

    df_res = data.loc[mask.any(axis=1)]
    df_res.reset_index(inplace=True, drop=True)
    with col4:
        st.write(f"Result in  {str(option)} table:", df_res)
    st.write(f"\n")

# -------------------------------------------------------------------------------------------------------------------
db_name = st.radio(
    "Find compounds in which DataSet",
    ("SelectChem and Jump DataSet", "KEGG"),
    horizontal=True,
)

df_cpds = pd.DataFrame()
sql_query = []
if len(df_res) > 0:
    if str(option) in table_mapping["SelectChem and Jump DataSet"]:
        table_name, col_name = table_mapping["SelectChem and Jump DataSet"][str(option)]
        list_geneid = [f"'{t}'" for t in df_res[col_name]]
        sql_query = get_sql_jump(
            table_name=table_name, col_name=col_name, list_geneid=list_geneid
        )

        df_cpds = (
            sql_df(sql_query, conn)
            .drop_duplicates(subset=["batchid"])
            .reset_index(drop=True)
        )

        if db_name == "KEGG":
            df_cpds = df_cpds.dropna(subset=["keggid"]).reset_index(drop=True)

# get genes are in crispr but not in cpd
df_prof_crisper = pd.DataFrame()
df_prof = pd.DataFrame()
if str(option) == "gene" and len(df_res) > 0:
    geneid_lis = [f"'{geneid}'" for geneid in df_res["geneid"]]

    sql_crisper = f"SELECT gene.symbol, gene.geneid,genebatchs.batchid  FROM genebatchs  inner join gene \
            on gene.geneid=genebatchs.geneid  WHERE genebatchs.geneid IN ({','.join(geneid_lis)})  group by gene.symbol, gene.geneid,genebatchs.batchid "
    df_crisperBatchs = sql_df(sql_crisper, conn).drop_duplicates(subset=["batchid"])
    batch_list = [f"'{batchid}'" for batchid in df_crisperBatchs["batchid"]]
    if len(batch_list) > 0:
        sql_crisper_profile = (
            f"SELECT * FROM aggprofile WHERE batchid IN ({','.join(batch_list)})"
        )

        df_prof_crisper = sql_df(sql_crisper_profile, conn)
        df_prof_crisper = df_prof_crisper.merge(
            df_crisperBatchs,
            on="batchid",
        ).reset_index(drop=True)
        df_prof_crisper["type"] = "CRISPR"

        df_prof_crisper["efficacy"] = "Unknown"
        df_prof_crisper["name"] = df_prof_crisper["symbol"]
        df_prof = pd.concat([df_prof, df_prof_crisper])

        df_prof["name"] = df_prof["name"].fillna(df_prof["batchid"])


# get cpd gene/target info-------------------------------------------------------------------------------------------------------------------

if len(df_cpds) > 0:
    list_pubchemid = [f"'{t}'" for t in df_cpds["pubchemid"]]
    list_batchid = [f"'{batch}'" for batch in df_cpds["batchid"]]

    # compounds GENE info
    sql = f"SELECT cpdgene.pubchemid, gene.geneid, gene.symbol, cpdgene.server FROM gene\
            INNER JOIN cpdgene ON cpdgene.geneid=gene.geneid\
            WHERE cpdgene.pubchemid IN ({','.join(list_pubchemid)})"

    df_cpdGene = (
        sql_df(sql, conn)
        .drop_duplicates(subset=["pubchemid", "geneid", "server"])
        .reset_index(drop=True)
    )
    df_cpdGene = df_cpdGene.merge(
        df_cpds[["cpdname", "batchid", "pubchemid"]],
        left_on="pubchemid",
        right_on="pubchemid",
    ).reset_index(drop=True)

    # compounds efficacy info
    sql = f"SELECT cpd.pubchemid, keggcpd.efficacy, keggcpd.keggid FROM keggcpd\
            INNER JOIN cpd ON keggcpd.keggid=cpd.keggid\
            WHERE cpd.pubchemid IN ({','.join(list_pubchemid)})"

    df_efficacy = (
        sql_df(sql, conn)
        .drop_duplicates(subset=["pubchemid", "efficacy", "keggid"])
        .reset_index(drop=True)
    )
    df_efficacy = df_efficacy.merge(
        df_cpds[["cpdname", "batchid", "pubchemid"]],
        left_on="pubchemid",
        right_on="pubchemid",
    ).reset_index(drop=True)
    df_efficacy = df_efficacy.drop_duplicates(
        subset=["pubchemid", "efficacy", "keggid"]
    ).reset_index(drop=True)

    # compounds PATHWAY info
    sql = f"SELECT cpdpath.pubchemid, pathway.pathid, pathway.pathname, cpdpath.server FROM pathway\
            INNER JOIN cpdpath ON cpdpath.pathid=pathway.pathid\
            WHERE cpdpath.pubchemid IN ({','.join(list_pubchemid)})"

    df_cpdPath = (
        sql_df(sql, conn)
        .drop_duplicates(subset=["pubchemid", "pathid", "server"])
        .reset_index(drop=True)
    )
    df_cpdPath = df_cpdPath.merge(
        df_cpds[["cpdname", "batchid", "pubchemid"]],
        left_on="pubchemid",
        right_on="pubchemid",
    ).reset_index(drop=True)
    df_cpdPath = df_cpdPath.drop_duplicates(
        subset=["pubchemid", "pathid", "server"]
    ).reset_index(drop=True)

    tabs = st.tabs(
        [
            "compounds info",
            "compounds GENE info",
            "compounds PATHWAY info",
            "compounds efficacy info",
            "compounds info in Pubchem",
        ]
    )
    tabs[0].write(df_cpds)

    with tabs[1]:
        server_name = st.radio(
            "select server",
            ("all", "pubchem", "KEGG", "broad", "SelectChem"),
            horizontal=True,
        )
        df_cpdGene_tmp = df_cpdGene.copy()
        if server_name != "all":
            df_cpdGene_tmp = df_cpdGene[
                df_cpdGene["server"] == server_name
            ].drop_duplicates(subset=["pubchemid", "geneid", "batchid"])
            df_cpdGene_tmp.reset_index(inplace=True, drop=True)
        fig_cols = st.columns(2)
        with fig_cols[0]:
            st.write(df_cpdGene_tmp)
        with fig_cols[1]:
            fig = px.pie(
                df_cpdGene_tmp,
                names="symbol",
                title=" symbol",
            )
            fig.update_traces(textposition="inside", textinfo="percent+label")
            st.plotly_chart(fig, theme="streamlit", use_container_width=True)
    with tabs[2]:
        fig_cols2 = st.columns(2)
        with fig_cols2[0]:
            st.write(df_cpdPath)
        with fig_cols2[1]:
            fig = px.pie(
                df_cpdPath,
                names="pathname",
                title=" pathname",
            )
            fig.update_traces(textposition="inside", textinfo="percent+label")
            st.plotly_chart(fig, theme="streamlit", use_container_width=True)
    with tabs[3]:
        fig_cols3 = st.columns(2)
        with fig_cols3[0]:
            st.write(df_efficacy)
        with fig_cols3[1]:
            fig = px.pie(
                df_efficacy,
                names="efficacy",
                title=" efficacy",
            )
            fig.update_traces(textposition="inside", textinfo="percent+label")
            st.plotly_chart(fig, theme="streamlit", use_container_width=True)
    with tabs[4]:
        rad_match = st.radio(
            "check in All compounds from compounds info ?can be long...",
            ["No", "Yes"],
            horizontal=True,
        )
        if rad_match == "No":
            var_text_cpdid = st.text_area(
                "Or just for few compounds", help="cpd id separated by enter"
            )
            var_cpd_id = var_text_cpdid.split("\n")
            var_cpd_id = [t.strip().upper() for t in var_cpd_id if t != ""]
            list_pub = list(set(var_cpd_id))

            if len(list_pub) > 0:
                st.session_state["PubChem"] = list_pub
                st.switch_page("pages/19_PubChem.py")
        else:
            list_pub = list(df_cpdGene.pubchemid.unique())

            st.switch_page("pages/19_PubChem.py")

    # get profile------------------------------------------------------------------------------------------------------------------
    st.write("Profiles in JUMP CP DATA SET")

    # get CPD profiles

    list_batchid = [f"'{batch}'" for batch in df_cpds["batchid"]]
    if server_name == "KEGG":  ## JUST KEGG
        list_batchid = [f"'{batch}'" for batch in df_cpdGene_tmp["batchid"]]
    if len(list_batchid) > 0:
        sql_profile = (
            "select * from aggprofile where batchid in (" + ",".join(list_batchid) + ")"
        )

        df_prof = sql_df(sql_profile, conn)
        df_prof.reset_index(inplace=True, drop=True)

        df_prof = df_prof.merge(df_cpds, on="batchid").reset_index(drop=True)

        df_prof.loc[df_prof.cpdname == "No result", "cpdname"] = None
        df_prof["name"] = df_prof["cpdname"].str[:30]
        df_prof["name"] = df_prof["cpdname"].fillna(df_prof["batchid"])
        if len(df_prof_crisper) > 0:
            df_prof = pd.concat([df_prof, df_prof_crisper])
            df_prof["name"] = df_prof["name"].fillna(df_prof["batchid"])

tab1, tab2, tab3, tab4 = st.tabs(
    [
        "compounds Profiles",
        "compounds Summary",
        "Crispr/ORF Profiles",
        "Crispr/ORF Summary",
    ]
)

if len(df_prof_crisper) > 0:
    tab3.write(df_prof_crisper)
    tab4.write(df_prof_crisper.describe().T)
    st.session_state["df_crisper"] = df_crisperBatchs


if len(df_prof) > 0:
    tab1.write(df_prof)
    tab2.write(df_prof.describe().T)

    list_sources = df_prof.source.unique().tolist()
    options2 = st.text_area("Enter sources")
    if options2:
        var_t2 = options2.split("\n")
        var_t2 = [t.strip() for t in var_t2]

    if not options2:
        tmp = df_prof.copy()
    else:
        tmp = df_prof.loc[df_prof["source"].isin(var_t2)]
    tmp.reset_index(inplace=True, drop=True)
    meta_cols = ["name", "batchid", "source", "symbol"]

    tmp["nameAndsource"] = tmp["name"] + "_" + tmp["source"]
    cpd_names = tmp.nameAndsource.values
    df_plt = tmp.drop_duplicates(subset="nameAndsource")
    df_plt = df_plt.set_index("nameAndsource")

    # df_plt = df_plt.drop('name',axis=1)

    filter_col = df_plt.select_dtypes(include=[int, float]).columns.tolist()
    df_plt2 = df_plt[filter_col].T

    sty = st.radio("Line Style", ["linear", "spline"])
    fig_line = px.line(
        df_plt2,
        x=filter_col,
        y=cpd_names,
        width=1400,
        height=1000,
        line_shape=sty,
        title="Profiles",
    )
    st.plotly_chart(fig_line, theme="streamlit", use_container_width=True)

    if len(tmp) > 1:
        import matplotlib.pyplot as plt
        import seaborn as sns

        plt_src, col_colors = get_col_colors(tmp, inex_col_name="nameAndsource")

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
            figsize=(16, len(plt_src) / 2),
        )

        st.pyplot(fig_clusmap)

    st.session_state["df_profiles"] = tmp

    @st.cache_data
    def get_umap(choix_source):
        sql_umqpemd = f"select umap.* from umap where umap.source='{choix_source}' "
        # sql_umqpemd="SELECT umap.*, batchs.name  FROM umap  INNER JOIN batchs ON batchs.batchid = umap.batchid  WHERE umap.source = '{choix_source}'"

        df_src_emd = sql_df(sql_umqpemd, conn)
        return df_src_emd

    # st.write(var_t2)
    if options2:
        df_src_emd = get_umap(choix_source=var_t2[0])
        df_src_emd["color"] = "others"
        df_src_emd.loc[
            df_src_emd["batchid"].isin(tmp["batchid"].to_list()), "color"
        ] = "selected compounds"

        fig = px.scatter(
            df_src_emd,
            x="umap1",
            y="umap2",
            color="color",
            opacity=0.5,
            color_discrete_sequence=["blue", "red", "green"],
            title="UMAP ",
            hover_data=["batchid"],
        )

        st.plotly_chart(fig, theme="streamlit", use_container_width=True)
