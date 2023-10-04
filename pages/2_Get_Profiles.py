import sys
import warnings

import numpy as np
import pandas as pd
import psycopg2
import streamlit as st

sys.path.append("/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/")
from streamlib import sql_df

st.set_page_config(
    layout="wide",
)
st.header("Search for Jump profiles", divider="rainbow")


warnings.filterwarnings("ignore")
# -------------------------------------------------------------------------
def convert_df(df):
    return df.to_csv(index=False).encode("utf-8")


def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])


conn = init_connection()

conn_profileDB = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksilink_cpds",
    password="12345",
)


def get_sql_kegg(table_name="keggcpdgene", col_name="geneid"):
    st.write(table_name)
    list_geneid = [f"'{t}'" for t in df_res[col_name]]
    sql_last_line = (
        f" WHERE UPPER({table_name}.{col_name} ) IN ({','.join(list_geneid)})"
    )

    if table_name == "keggcpd":
        sql = f"SELECT * FROM keggcpd {sql_last_line}"
    else:
        sql = f"SELECT  keggcpd.keggid,keggcpd.keggcpdname , {table_name}.{col_name} FROM keggcpd\
                INNER JOIN {table_name} ON {table_name}.keggid=keggcpd.keggid {sql_last_line}\
                GROUP BY keggcpd.keggcpdname, keggcpd.keggid, {table_name}.{col_name}"
    return sql


def get_sql_jump(table_name="cpdgene", col_name="geneid"):
    list_geneid = [f"'{t}'" for t in df_res[col_name]]
    sql_last_line = (
        f" WHERE UPPER({table_name}.{col_name}) IN ({','.join(list_geneid)})"
    )
    sql_first_line = "SELECT cpd.pubchemid, cpdbatchs.batchid, cpd.synonyms, cpd.keggid, cpd.cpdname, cpd.smile"

    if table_name == "keggcpddis":
        # st.write(table_name)
        sql = f"{sql_first_line}, keggcpddis.disid FROM cpdbatchs\
                INNER JOIN cpd ON cpdbatchs.pubchemid=cpd.pubchemid\
                INNER JOIN keggcpddis ON keggcpddis.keggid=cpd.keggid {sql_last_line}"
    elif table_name in ["cpd", "cpdbatchs"]:
        # st.write(table_name)
        sql = f"{sql_first_line} FROM cpdbatchs right JOIN cpd ON cpd.pubchemid=cpdbatchs.pubchemid {sql_last_line}"
    else:
        # st.write(table_name)
        sql = f"{sql_first_line}, {table_name}.{col_name} FROM cpdbatchs\
                INNER JOIN {table_name} ON {table_name}.pubchemid=cpdbatchs.pubchemid\
                INNER JOIN cpd ON cpdbatchs.pubchemid=cpd.pubchemid {sql_last_line}\
                GROUP BY cpd.pubchemid, cpdbatchs.batchid, cpd.synonyms, cpd.keggid, cpd.cpdname, cpd.smile, {table_name}.{col_name}"
    return sql


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

    for col in df_plt.columns:
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
    return df_plt, col_colors


# -------------------------------------------------------------------------
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
    var_text = st.text_area("Enter your search", help="Name or ID separated by enter")
var_t = var_text.split("\n")
var_t = [t.strip().upper() for t in var_t]
df_res = pd.DataFrame()
if len(var_t) > 0 and var_t[0] != "":
    mask = np.column_stack(
        [
            data[col].str.match(v.upper(), na=False)
            for col in data.columns
            for v in var_t
            if v != ""
        ]
    )
    df_res = data.loc[mask.any(axis=1)]
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
        sql_query = get_sql_jump(table_name=table_name, col_name=col_name)
        df_cpds = (
            sql_df(sql_query, conn)
            .drop_duplicates(subset=["batchid"])
            .reset_index(drop=True)
        )
        if db_name == "KEGG":
            # sql_query = get_sql_kegg(table_name=table_name, col_name=col_name)

            df_cpds = df_cpds.dropna(subset=["keggid"]).reset_index(drop=True)
            # st.write(df_cpds)

        # else:
        #     # sql_query = get_sql_jump(table_name=table_name, col_name=col_name)

        #     df_cpds = sql_df(sql_query, conn).drop_duplicates(subset=["batchid"]).reset_index(drop=True)

# get crisper profiles when search gene
df_prof_crisper = None
if str(option) == "gene" and len(df_res) > 0:
    geneid_lis = [f"'{geneid}'" for geneid in df_res["geneid"]]
    sql_crisper = f"SELECT * FROM crisperbatchs WHERE crisperbatchs.geneid IN ({','.join(geneid_lis)})"
    df_crisperBatchs = sql_df(sql_crisper, conn).drop_duplicates(subset=["batchid"])
    batch_list = [f"'{batchid}'" for batchid in df_crisperBatchs["batchid"]]
    if len(batch_list) >0:
        sql_crisper_profile = (
            f"SELECT * FROM aggcombatprofile WHERE metabatchid IN ({','.join(batch_list)})"
        )
       
        df_prof_crisper = sql_df(sql_crisper_profile, conn_profileDB)
        crisper_dic = df_crisperBatchs.set_index("batchid")["geneid"].to_dict()
        df_prof_crisper["metageneid"] = df_prof_crisper["metabatchid"].map(crisper_dic)


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
        ]
    )
    tabs[0].write(df_cpds)
    st.download_button(
        label="Save compounds",
        data=convert_df(df_cpds.drop_duplicates(subset=["pubchemid"])),
        file_name="cpds.csv",
        mime="csv",
    )
    with tabs[1]:
        server_name = st.radio(
            "select server",
            ("all", "pubchem", "KEGG"),
            horizontal=True,
        )
        df_cpdGene_tmp = df_cpdGene.copy()
        if server_name != "all":
            df_cpdGene_tmp = df_cpdGene[
                df_cpdGene["server"] == server_name
            ].drop_duplicates(subset=["pubchemid", "geneid", "batchid"])
        # if server_name == "KEGG":
        #     df_cpdGene_tmp = df_cpds[df_cpdGene["server"] == server_name].drop_duplicates(subset=["pubchemid","geneid"])
    tabs[1].write(df_cpdGene_tmp)
    tabs[2].write(df_cpdPath)
    tabs[3].write(df_efficacy)

    # get profile------------------------------------------------------------------------------------------------------------------

    st.write(f"Profiles in JUMP CP DATA SET")
    df_prof = pd.DataFrame()
    if server_name == "KEGG":
        list_batchid = [f"'{batch}'" for batch in df_cpdGene_tmp["batchid"]]
    else:
        list_batchid = [f"'{batch}'" for batch in df_cpds["batchid"]]
    sql_profile = (
        "select * from aggcombatprofile where metabatchid in ("
        + ",".join(list_batchid)
        + ")"
    )

    df_prof = sql_df(sql_profile, conn_profileDB)
    df_prof.reset_index(inplace=True, drop=True)

    df_prof = df_prof.merge(
        df_cpds.add_prefix("meta"), left_on="metabatchid", right_on="metabatchid"
    ).reset_index(drop=True)

    df_prof.loc[df_prof.metacpdname == "No result", "metacpdname"] = None
    df_prof["metacpdname"] = df_prof["metacpdname"].str[:30]
    df_prof["metacpdname"] = df_prof["metacpdname"].fillna(df_prof["metabatchid"])

    tab1, tab2, tab3, tab4 = st.tabs(
        [
            "compounds Profiles",
            "compounds Summary",
            "Crisper Profiles",
            "Crisper Summary",
        ]
    )
    tab1.write(df_prof)
    st.download_button(
        label="Save Profile",
        data=convert_df(df_prof),
        file_name="df_prof.csv",
        mime="csv",
    )
    tab2.write(df_prof.describe().T)
    if str(option) == "gene":
        df_prof_crisper["metatype"] = "CRISPR"
        df_prof_crisper["metaefficacy"] = "Unknown"
        df_prof = pd.concat([df_prof, df_prof_crisper])
        df_prof["metacpdname"] = df_prof["metacpdname"].fillna(df_prof["metabatchid"])
        tab3.write(df_prof_crisper)
        tab4.write(df_prof_crisper.describe().T)

    # plot------------------------------------------------------------------------------------------------------------------
    list_sources = df_prof.metasource.unique().tolist()
    options = st.text_area("Enter sources")
    var_t = options.split("\n")
    var_t = [t.strip() for t in var_t]

    if not options:
        tmp = df_prof.copy()
    else:
        tmp = df_prof.loc[df_prof["metasource"].isin(var_t)]

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
            vmin=-10,
            vmax=10,
            figsize=(16, len(plt_src) / 2),
        )

        st.pyplot(fig_clusmap)
        # st.write(len(plt_src)/2)

    st.session_state["df_profiles"] = df_prof
    st.session_state["df_cpds"] = df_cpds
    st.session_state["df_cpdGene"] = df_cpdGene
    st.session_state["df_cpdPath"] = df_cpdPath
    st.session_state["df_efficacy"] = df_efficacy
    st.session_state["df_crisper"] = df_prof_crisper


else:
    st.warning(" No Luck!! ")

# cleas cache
# data=None
# df_res=None
# df_cpds=None
# df_prof_crisper=None
# df_cpdGene=None
# df_efficacy=None
# df_cpdPath=None
# df_prof=None
conn_profileDB.close()
conn.close()
