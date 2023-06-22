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


# rows = run_query(
#     "select table_name from information_schema.tables where table_schema='public'"
# )
# list_tables = []

# for row in rows:
#     list_tables.append(str(row).split(",")[0].split("(")[1].replace("'", ""))
list_tables = ["gene", "cpd", "pathway", "disease", "keggcpd"]

# -------------------------------------------------------------------------------------------------------------------
st.image("DB.png")
st.write(f"\n")
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
df_cp = []

def get_last_line(col_name="geneid",table_name="keggcpdgene"):
    list_tmp = df_res[col_name]
    list_geneid = []
    for t in list_tmp:
        list_geneid.append("'" + t + "'")
    sql_last_line = f" where {table_name}.{col_name} in (" + ",".join(list_geneid) + ")"
    return sql_last_line

if db_name == "KEGG":
    st.write(f"Result in KEGG drugs and compunds  table:")
    list_df = []
    sql_cpd = None
    df_cp = []
    if len(df_res) > 0:
        if str(option) == "gene":
            sql_last_line = get_last_line(col_name="geneid",table_name="keggcpdgene")
            
            sql_cpd = (
                'select keggcpd.name as "KEGG Compound Name", keggcpd.keggid, keggcpdgene.geneid   '
                + " from keggcpd "
                + " inner join keggcpdgene on keggcpdgene.keggid=keggcpd.keggid "
                + sql_last_line
                + " GROUP BY keggcpd.name, keggcpd.keggid, keggcpdgene.geneid "
            )

        if str(option) == "pathway":
            sql_last_line = get_last_line(col_name="pathid",table_name="keggcpdpath")
            sql_cpd = (
                'select keggcpd.name as "KEGG Compound Name", keggcpd.keggid, keggcpdpath.pathid '
                + " from keggcpd "
                + " inner join keggcpdpath on keggcpdpath.keggid=keggcpd.keggid "
                + sql_last_line
                + " GROUP BY keggcpd.name, keggcpd.keggid, keggcpdpath.pathid"
            )

        if str(option) == "disease":
            sql_last_line = get_last_line(col_name="disid",table_name="keggcpddis")
            sql_cpd = (
                'select keggcpd.name as "KEGG Compound Name", keggcpd.keggid, keggcpddis.disid '
                + " from keggcpd "
                + " inner join keggcpddis on keggcpddis.keggid=keggcpd.keggid "
                + sql_last_line
                + " GROUP BY keggcpd.name, keggcpd.keggid, keggcpddis.disid"
            )


if db_name == "SelectChem and Jump DataSet":
    st.write(f"Result in SelectChem compunds  table:")
    list_df = []
    sql_cpd = None
    df_cp = []
    if len(df_res) > 0:
        if str(option) == "gene":
            sql_last_line = get_last_line(col_name="geneid",table_name="cpdgene")
            sql_cpd = (
                'select cpd.synonyms as "Compound Name"  ,cpd.pubchemid as "Compound PubChemID" ,cpdbatchs.batchid as "Compound BatchID" ,cpd.name as "Compound FullName",cpdgene.geneid '
                + " from cpdbatchs "
                + " inner join cpdgene on cpdgene.pubchemid=cpdbatchs.pubchemid "
                + " inner join cpd on cpdbatchs.pubchemid=cpd.pubchemid"
                + sql_last_line
                + " GROUP BY cpd.pubchemid,cpdbatchs.batchid,cpdgene.geneid, cpd.name"
            )

        if str(option) == "pathway":
            sql_last_line = get_last_line(col_name="pathid",table_name="cpdpath")
            
            sql_cpd = (
                'select cpd.synonyms as "Compound synonyms"  ,cpd.pubchemid as "Compound PubChemID" ,cpdbatchs.batchid as "Compound BatchID" ,cpd.name as "Compound FullName",cpdpath.pathid '
                + " from cpdbatchs "
                + " inner join cpdpath on cpdpath.pubchemid=cpdbatchs.pubchemid "
                + " inner join cpd on cpdbatchs.pubchemid=cpd.pubchemid"
                + sql_last_line
                + " GROUP BY cpd.pubchemid,cpdbatchs.batchid,cpdpath.pathid, cpd.name"
            )

        if str(option) == "disease":
            sql_last_line = get_last_line(col_name="disid",table_name="keggcpddis")
            sql_cpd = (
                'select cpd.synonyms as "Compound synonyms"  ,cpd.pubchemid as "Compound PubChemID" ,cpdbatchs.batchid as "Compound BatchID" ,cpd.name as "Compound FullName",keggcpddis.disid '
                + " from cpdbatchs "
                + " INNER JOIN cpd ON cpdbatchs.pubchemid=cpd.pubchemid "
                + " INNER JOIN keggcpddis ON keggcpddis.keggid=cpd.keggid"
                + sql_last_line
                + " GROUP BY cpd.pubchemid,cpdbatchs.batchid,keggcpddis.disid, cpd.name"
            )

        if str(option) == "cpd":
            sql_last_line = get_last_line(col_name="pubchemid",table_name="cpd")
 

            sql_cpd = (
                'select cpd.synonyms as "Compound synonyms"  ,cpd.pubchemid as "Compound PubChemID" ,cpdbatchs.batchid as "Compound BatchID" ,cpd.name as "Compound FullName" '
                + " from cpdbatchs "
                + " INNER JOIN cpd ON cpdbatchs.pubchemid=cpd.pubchemid "
                + sql_last_line
                + " GROUP BY cpd.pubchemid,cpdbatchs.batchid,cpd.pubchemid, cpd.name"
            )


df_cp = pd.read_sql(sql_cpd, conn)


if len(df_cp) > 0:
    st.write(df_cp)
    st.write(f"Profiles ")
    if db_name == "SelectChem and Jump DataSet":
        df_prof = pd.DataFrame()
        list_batch = df_cp["Compound BatchID"]
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
