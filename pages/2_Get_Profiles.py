
import warnings

import numpy as np
import pandas as pd
import psycopg2
import streamlit as st

warnings.filterwarnings("ignore")


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


list_tables = ["gene", "cpd", "pathway", "disease","cpdbatchs"]
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
    var_text = st.text_area("Enter your search",help='Name or ID separated by enter')
    
var_t=var_text.split('\n')
var_t = [t.strip().upper() for t in var_t]

df_res=pd.DataFrame()

# for var_text in var_t:
#     mask = np.column_stack([data[col].str.match(var_text.upper(), na=False) for col in data])
# df_res= data.loc[mask.any(axis=1)]

if len(var_t)>0 and var_t[0]!="":
    mask = np.column_stack([data[col].str.match(v.upper(), na=False) for col in data.columns for v in var_t if v!=''])
    df_res= data.loc[mask.any(axis=1)]
    with col4:
        st.write(f"Result in  {str(option)} table:", df_res)
    st.write(f"\n")
# -------------------------------------------------------------------------------------------------------------------

db_name = st.radio(
    "Find compounds in which DataSet",
    ("SelectChem and Jump DataSet", "KEGG"),
    horizontal=True,
)



def get_sql_kegg(table_name="keggcpdgene", col_name="geneid"):
    list_geneid = [f"'{t}'" for t in df_res[col_name]]
    sql_last_line = f" WHERE UPPER({table_name}.{col_name} ) IN ({','.join(list_geneid)})"

    if table_name == "keggcpd":
        sql = f"SELECT * FROM keggcpd {sql_last_line}"
    else:
        sql = f"SELECT  keggcpd.keggid,keggcpd.name , {table_name}.{col_name} FROM keggcpd INNER JOIN {table_name} ON {table_name}.keggid=keggcpd.keggid {sql_last_line} GROUP BY keggcpd.name, keggcpd.keggid, {table_name}.{col_name}"

    return sql


def get_sql_jump(table_name="cpdgene", col_name="geneid"):
    list_geneid = [f"'{t}'" for t in df_res[col_name]]
    sql_last_line = f" WHERE UPPER({table_name}.{col_name}) IN ({','.join(list_geneid)})"
    sql_first_line = "SELECT cpd.pubchemid, cpdbatchs.batchid, cpd.synonyms, cpd.keggid, cpd.name, cpd.smile"

    if table_name == "keggcpddis":
        st.write(table_name)
        sql = f"{sql_first_line}, keggcpddis.disid FROM cpdbatchs INNER JOIN cpd ON cpdbatchs.pubchemid=cpd.pubchemid INNER JOIN keggcpddis ON keggcpddis.keggid=cpd.keggid {sql_last_line}"
    elif table_name in[ "cpd","cpdbatchs" ]:
        st.write(table_name)
        sql = f"{sql_first_line} FROM cpdbatchs right JOIN cpd ON cpd.pubchemid=cpdbatchs.pubchemid {sql_last_line}"
    else:
        st.write(table_name)
        sql = f"{sql_first_line}, {table_name}.{col_name} FROM cpdbatchs INNER JOIN {table_name} ON {table_name}.pubchemid=cpdbatchs.pubchemid INNER JOIN cpd ON cpdbatchs.pubchemid=cpd.pubchemid {sql_last_line} GROUP BY cpd.pubchemid, cpdbatchs.batchid, cpd.synonyms, cpd.keggid, cpd.name, cpd.smile, {table_name}.{col_name}"

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


df_cpds = []
sql_query = []
if len(df_res) > 0:
    if str(option) in table_mapping[db_name]:
        table_name, col_name = table_mapping[db_name][str(option)]
        if db_name == "KEGG":
            sql_query = get_sql_kegg(table_name=table_name, col_name=col_name)

            df_kegg = pd.read_sql(sql_query, conn).drop_duplicates(subset=["keggid"]).reset_index(drop=True)
            st.write(df_kegg)

        else:
            sql_query = get_sql_jump(table_name=table_name, col_name=col_name)
           
            df_cpds = pd.read_sql(sql_query, conn).drop_duplicates(subset=["batchid"]).reset_index(drop=True)

# get crisper profiles when search gene
df_prof_crisper = None
if str(option) == "gene" and len(df_res)>0:
    geneid_lis = [f"'{geneid}'" for geneid in df_res["geneid"]]
    sql_crisper = f"SELECT * FROM crisperbatchs WHERE crisperbatchs.geneid IN ({','.join(geneid_lis)})"
    df_crisperBatchs = pd.read_sql(sql_crisper, conn).drop_duplicates(subset=["batchid"])
    batch_list = [f"'{batchid}'" for batchid in df_crisperBatchs["batchid"]]
    sql_crisper_profile = f"SELECT * FROM aggcombatprofile WHERE metabatchid IN ({','.join(batch_list)})"
    df_prof_crisper = pd.read_sql(sql_crisper_profile, conn_profileDB)
    crisper_dic = df_crisperBatchs.set_index("batchid")["geneid"].to_dict()
    df_prof_crisper["metageneid"] = df_prof_crisper["metabatchid"].map(crisper_dic)

    
    
    

if len(df_cpds) > 0:
    list_pubchemid = [f"'{t}'" for t in df_cpds["pubchemid"]]
    list_batchid = [f"'{batch}'" for batch in df_cpds["batchid"]]
    

    # compounds GENE info
    sql = f"SELECT cpdgene.pubchemid, gene.geneid, gene.symbol, cpdgene.server FROM gene INNER JOIN cpdgene ON cpdgene.geneid=gene.geneid WHERE cpdgene.pubchemid IN ({','.join(list_pubchemid)})"
    df_cpdGene = pd.read_sql(sql, conn).drop_duplicates(subset=["pubchemid", "geneid", "server"]).reset_index(drop=True)
    df_cpdGene = df_cpdGene.merge(df_cpds[["batchid","pubchemid"]],left_on='pubchemid',right_on='pubchemid').reset_index(drop=True)


    # compounds efficacy info
    sql = f"SELECT cpd.pubchemid, keggcpd.efficacy, keggcpd.keggid FROM keggcpd INNER JOIN cpd ON keggcpd.keggid=cpd.keggid WHERE cpd.pubchemid IN ({','.join(list_pubchemid)})"
    df_efficacy = pd.read_sql(sql, conn).drop_duplicates(subset=["pubchemid", "efficacy", "keggid"]).reset_index(drop=True)
    df_efficacy = df_efficacy.merge(df_cpds[["batchid","pubchemid"]],left_on='pubchemid',right_on='pubchemid').reset_index(drop=True)

    # compounds PATHWAY info
    sql = f"SELECT cpdpath.pubchemid, pathway.pathid, pathway.name, cpdpath.server FROM pathway INNER JOIN cpdpath ON cpdpath.pathid=pathway.pathid WHERE cpdpath.pubchemid IN ({','.join(list_pubchemid)})"
    df_cpdPath = pd.read_sql(sql, conn).drop_duplicates(subset=["pubchemid", "pathid", "server"]).reset_index(drop=True)
    df_cpdPath = df_cpdPath.merge(df_cpds[["batchid","pubchemid"]],left_on='pubchemid',right_on='pubchemid').reset_index(drop=True)

    tabs = st.tabs(["compounds info", "compounds GENE info", "compounds PATHWAY info", "compounds efficacy info"])
    tabs[0].write(df_cpds)
    with tabs[1]:
        server_name = st.radio(
        "select server",
        ("all", "pubchem", "KEGG"),
        horizontal=True,
    )
        df = df_cpdGene
        if server_name != "all":
            df = df_cpdGene[df_cpdGene["server"] == server_name]

    tabs[1].write(df)
    tabs[2].write(df_cpdPath)
    tabs[3].write(df_efficacy)




    st.write(f"Profiles in JUMP CP DATA SET")

    df_prof = pd.DataFrame()


    
    sql_profile = (
        "select * from aggcombatprofile where metabatchid in ("
        + ",".join(list_batchid)
        + ")"
    )

    df_prof = pd.read_sql(sql_profile, conn_profileDB)

    dic_efficacy = df_efficacy.set_index("batchid").to_dict()["efficacy"]
    dic_Gene = df_cpdGene.set_index("batchid").to_dict()["geneid"]
    df_prof["metageneid"] = df_prof["metabatchid"].map(dic_Gene)
    df_prof["metaefficacy"] = df_prof["metabatchid"].map(dic_efficacy)
    df_prof["metaefficacy"].fillna("Unknown", inplace=True)
    df_prof["metatype"] = "CPDS"



    tab1, tab2,tab3,tab4 = st.tabs(["compounds Profiles", "compounds Summary","Crisper Profiles", "Crisper Summary"])
    tab1.write(df_prof)
    tab2.write(df_prof.describe().T)
    if str(option)=="gene":
        df_prof_crisper["metatype"] = "CRISPR"
        df_prof_crisper["metaefficacy"] = "Unknown"
        tab3.write(df_prof_crisper)
        tab4.write(df_prof_crisper.describe().T)

    #-------------------------------------------------------------
    

    list_sources = df_prof.metasource.unique().tolist()

    options = st.multiselect('Select Code', list_sources )
    st.write(options)
    if not options:
        tmp = df_prof.copy()
    else:
        tmp = df_prof.loc[data["metasource"].isin(options)]
  

    if len(tmp)>0:
        import matplotlib.pyplot as plt
        import seaborn as sns
        list_col = [col for col in tmp.columns if not col.startswith("meta")]
        plt_src = pd.DataFrame()
        plt_src = tmp[list_col]
        plt_src["name"] = tmp["metabatchid"] + "_" + tmp["metasource"]
        plt_src.set_index("name", inplace=True)
        fig_clusmap, ax1 = plt.subplots()
        fig_clusmap = sns.clustermap(
                    plt_src,
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
            
        st.pyplot(fig_clusmap)
            
    st.session_state["df_profiles"] = df_prof
    st.session_state["df_cpds"] = df_cpds
    st.session_state["df_cpdGene"] = df_cpdGene
    st.session_state["df_cpdPath"] = df_cpdPath
    st.session_state["df_efficacy"] = df_efficacy
    st.session_state["df_crisper"] = df_prof_crisper



else:
    st.warning(" No Luck!! ")

# tmap
sql = "SELECT * FROM tmap"
df_tmap = pd.read_sql(sql, conn_profileDB).reset_index(drop=True)
st.session_state["df_tmap"] = df_tmap
conn_profileDB.close()