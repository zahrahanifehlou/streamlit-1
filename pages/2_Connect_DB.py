import warnings
import numpy as np
import pandas as pd
import psycopg2
import streamlit as st

warnings.filterwarnings("ignore")


def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])


conn = init_connection()


def run_query(query):
    with conn.cursor() as cur:
        cur.execute(query)
        return cur.fetchall()


list_tables = ["gene", "cpd", "pathway", "disease", "keggcpd"]
# -------------------------------------------------------------------------------------------------------------------

col1, col2 = st.columns(2)
with col1:
    option = st.selectbox("Pick one table", list_tables)
    x = st.slider("Number of rows to display", 1, 20)
    st.write("You selected the table: ", option)
data = pd.read_sql("select * from " + str(option), conn)
data = data.apply(lambda x: x.astype(str).str.upper())

with col2:
    tab1, tab2 = st.tabs(["Dataframe", "Summary"])
    tab1.write(data.sample(x))
    tab2.write(data.describe().T)
st.write(f"\n")
# -------------------------------------------------------------------------------------------------------------------

col3, col4 = st.columns(2)
with col3:
    var_text = st.text_input("Enter your search")

mask = np.column_stack(
    [data[col].str.match(var_text.upper(), na=False) for col in data]
)

df_res = data.loc[mask.any(axis=1)].reset_index(drop=True)

with col4:
    st.write(f"Result in  {str(option)} table:", df_res)
st.write(f"\n")
# -------------------------------------------------------------------------------------------------------------------

db_name = st.radio(
    "Find compounds in which DataSet",
    ("SelectChem and Jump DataSet", "KEGG"),
    horizontal=True,
)
df_cpds = []


def get_sql_kegg(table_name="keggcpdgene", col_name="geneid"):
    list_tmp = df_res[col_name]
    list_geneid = []
    for t in list_tmp:
        list_geneid.append("'" + t + "'")
    sql_last_line = (
        " where  " + table_name + "." + col_name + " in (" + ",".join(list_geneid) + ")"
    )
    sql = (
        'select keggcpd.name as "KEGG Compound Name", keggcpd.keggid, '
        + table_name
        + "."
        + col_name
        + " from keggcpd  inner join "
        + table_name
        + " on "
        + table_name
        + ".keggid=keggcpd.keggid "
        + sql_last_line
        + " GROUP BY keggcpd.name, keggcpd.keggid, "
        + table_name
        + "."
        + col_name
    )
    return sql


def get_sql_jump(table_name="cpdgene", col_name="geneid"):
    list_tmp = df_res[col_name]
    list_geneid = []
    for t in list_tmp:
        list_geneid.append("'" + t + "'")
    sql_last_line = (
        " where  " + table_name + "." + col_name + " in (" + ",".join(list_geneid) + ")"
    )
    sql_first_line = "select cpd.pubchemid,cpd.synonyms,cpd.keggid, cpd.name, cpd.smile ,cpdbatchs.batchid "

    if table_name == "disease":
        st.write(table_name)
        sql = (
            sql_first_line
            + ", keggcpddis.disid "
            + " from cpdbatchs inner join cpd on cpdbatchs.pubchemid=cpd.pubchemid "
            + " inner  join keggcpd on keggcpd.keggid=cpd.keggid"
            + " inner  join keggcpddis on keggcpd.keggid=keggcpddis.keggid"
            + sql_last_line
        )

    if table_name == "cpd":
        st.write(table_name)

        sql = (
            sql_first_line
            + table_name
            + "."
            + col_name
            + " from cpdbatchs "
            + " inner join "
            + table_name
            + " on "
            + table_name
            + ".pubchemid=cpdbatchs.pubchemid   "
            + sql_last_line
            + " GROUP BY cpd.pubchemid,cpdbatchs.batchid,"
            + table_name
            + "."
            + col_name
            + ", cpd.name"
        )

    if table_name not in ["cpd", "disease"]:
        sql = (
            sql_first_line
            + ","
            + table_name
            + "."
            + col_name
            + " from cpdbatchs "
            + " inner join "
            + table_name
            + " on "
            + table_name
            + ".pubchemid=cpdbatchs.pubchemid   inner join cpd on cpdbatchs.pubchemid=cpd.pubchemid"
            + sql_last_line
            + " GROUP BY cpd.pubchemid,cpdbatchs.batchid,"
            + table_name
            + "."
            + col_name
            + ", cpd.name"
        )

    return sql


table_mapping = {
    "KEGG": {
        "gene": ("keggcpdgene", "geneid"),
        "pathway": ("keggcpdpath", "pathid"),
        "disease": ("keggcpddis", "disid"),
        "cpdkegg": ("cpdkegg", "keggid"),
    },
    "SelectChem and Jump DataSet": {
        "gene": ("cpdgene", "geneid"),
        "pathway": ("cpdpath", "pathid"),
        "disease": ("disease", "disid"),
        "cpd": ("cpd", "pubchemid"),
    },
}

list_df = []
df_cpds = []
sql_query = []
if len(df_res) > 0:
    if str(option) in table_mapping[db_name]:
        table_name, col_name = table_mapping[db_name][str(option)]
        if db_name == "KEGG":
            sql_query = get_sql_kegg(table_name=table_name, col_name=col_name)
        else:
            sql_query = get_sql_jump(table_name=table_name, col_name=col_name)


if sql_query != []:
    df_cpds = pd.read_sql(sql_query, conn)
    df_cpds = df_cpds.drop_duplicates(subset=["batchid"])
    df_cpds = df_cpds.reset_index(drop=True)

# get KEGG info
df_kegg = df_cpds[df_cpds["keggid"].notna()]
# get KEGG GENE and PATHWAY info
if len(df_cpds) > 0:
    list_batchid = []
    for t in df_cpds["pubchemid"]:
        list_batchid.append("'" + t + "'")
    sql_kegg_gene = (
        "select gene.geneid, gene.symbol,cpdgene.pubchemid,cpdgene.server from gene inner join  cpdgene on cpdgene.geneid=gene.geneid where  cpdgene.pubchemid"
        + " in ("
        + ",".join(list_batchid)
        + ")"
    )
    df_keggGene = pd.read_sql(sql_kegg_gene, conn)
    df_keggGene = df_keggGene.drop_duplicates(subset=["pubchemid","geneid","server"])
    df_keggGene = df_keggGene.reset_index(drop=True)

    sql_kegg_pathway = (
        "select pathway.pathid, pathway.name,cpdpath.pubchemid, cpdpath.server from pathway inner join  cpdpath on cpdpath.pathid=pathway.pathid  where  cpdpath.pubchemid"
        + " in ("
        + ",".join(list_batchid)
        + ")"
    )
    df_keggPath = pd.read_sql(sql_kegg_pathway, conn)
    df_keggPath = df_keggPath.drop_duplicates(subset=["pubchemid","pathid","server"])
    df_keggPath = df_keggPath.reset_index(drop=True)


if len(df_cpds) > 0:
    tab1, tab2, tab3= st.tabs(["compounds info", "compounds GENE info",  "compounds PATHWAY info"])
    tab1.write(df_cpds)
    tab2.write(df_keggGene)
    tab3.write(df_keggPath)
    
    st.write(f"Profiles in JUMP CP DATA SET")
    if db_name == "SelectChem and Jump DataSet":
        df_prof = pd.DataFrame()
        list_batch = df_cpds["batchid"]
        batch_list = []
        for batch in list_batch:
            batch_list.append("'" + batch + "'")

        conn2 = psycopg2.connect(
            host="192.168.2.131",
            port="5432",
            user="arno",
            database="ksilink_cpds",
            password="12345",
        )
        sql_profile = (
            "select * from aggcombatprofile where metabatchid in ("
            + ",".join(batch_list)
            + ")"
        )

        df_prof = pd.read_sql(sql_profile, conn2)
        tab1, tab2 = st.tabs(["Profiles", "Summary"])
        tab1.write(df_prof)
        tab2.write(df_prof.describe().T)
     

        st.session_state["df_profiles"] = df_prof
        st.session_state["df_cpds"] = df_cpds
        st.session_state["df_kegg"] = df_kegg
        st.session_state["df_keggGene"] = df_keggGene
        st.session_state["df_keggPath"] = df_keggPath
else:
    st.write(" try  something else :) ")
