import warnings

import numpy as np
import pandas as pd
import psycopg2
import streamlit as st

warnings.filterwarnings("ignore")


# Initialize connection.
# Uses st.cache_resource to only run once.
@st.cache_resource
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])


conn = init_connection()


# Perform query.
# Uses st.cache_data to only rerun when the query changes or after 10 min.
@st.cache_data(ttl=600)
def run_query(query):
    with conn.cursor() as cur:
        cur.execute(query)
        return cur.fetchall()


rows = run_query("select table_name from information_schema.tables where table_schema='public'")
list_tables = []
st.image("DB.png")
for row in rows:
    list_tables.append(str(row).split(",")[0].split("(")[1].replace("'", ""))

option = st.selectbox("Pick one table", list_tables)
st.write("You selected the table: ", option)
data = pd.read_sql("select * from " + str(option), conn)
data = data.apply(lambda x: x.astype(str).str.upper())

x = st.slider("Number of rows to display", 1, 20)
tab1, tab2 = st.tabs(["Dataframe", "Summary"])
tab1.write(data.sample(x))
tab2.write(data.describe())

var_text = st.text_input("Enter your search")

# df_res = data[data.eq(var_text.upper()).any(1)]
mask = np.column_stack([data[col].str.match(var_text.upper(), na=False) for col in data])
# mask = np.column_stack([data[col].to_string()==var_text.upper() for col in data])
# st.write(mask)
df_res= data.loc[mask.any(axis=1)]
# result = data[data.applymap(lambda x: var_text.upper() in str(x))]
#filt= result[result.any(axis=1)]
#df_res =data[data[filt.any(axis=1)]]
st.write(f"Result in selected {str(option)} table:", df_res)
list_df = []
g = st.radio("Select gene or compound", ["gene", "compound"])

for tab in list_tables:
    #st.write(tab)
    dfx = pd.read_sql("select * from " + str(tab), conn)
    dfx = dfx.apply(lambda x: x.astype(str).str.upper())
    mask = np.column_stack([dfx[col].str.match(var_text.upper(), na=False) for col in dfx])

    df_resx= dfx.loc[mask.any(axis=1)]
    # st.write(df_resx)
    # df_resx = dfx[dfx.eq(var_text.upper()).any(1)]
    if len(df_resx) > 0:
        st.write(f"Result in {str(tab)}  table:", df_resx)

        for col in df_resx:
            sql_last_line = f" where UPPER({tab}.{col}) like   UPPER('{var_text}%') "
            if g == "gene":
                sql_crisper = (
                    "select crispermeta.batchid,gene.geneid,gene.symbol,gene.synonyms"
                    + " from crispermeta "
                    + "inner join gene on crispermeta.geneid=gene.geneid "
                    + sql_last_line
                    + ""
                    + " GROUP BY crispermeta.batchid,gene.geneid,gene.symbol,gene.synonyms"
                )
                df_cris = pd.read_sql(sql_crisper, conn)
                if len(df_cris)>0:
                    list_df.append(df_cris)
                sql_cpd = (
                    "select  batchs.batchid, gene.geneid,gene.symbol,gene.synonyms"
                    + " from batchs"
                    + " inner join cpdgene on cpdgene.pubchemid=batchs.pubchemid"
                    + " inner join gene on cpdgene.geneid=gene.geneid"
                    + sql_last_line
                    + " GROUP BY batchs.batchid, gene.geneid,gene.symbol,gene.synonyms"
                )
                df_cp = pd.read_sql(sql_cpd, conn)
                if len(df_cp)>0:
                    list_df.append(df_cp)
            else:
                sql_cpd = (
                    "select  batchs.batchid, cpd.pubchemid,cpd.name "
                    + " from cpd"
                    + " inner join batchs on batchs.pubchemid=cpd.pubchemid "
                    + sql_last_line
                    # + " GROUP BY batchs.batchid, cpd.pubchemid,cpd.name"
                )
                df_cps = pd.read_sql(sql_cpd, conn)
                # st.write("Data for cpds",sql_cpd)
                if len(df_cps)>0:
                    list_df.append(df_cps)
# st.write(list_df)
if len(list_df) > 0:
    df_prof = pd.DataFrame()
    df_tot = pd.concat(list_df)
    # st.write(list_df)
    st.write(f"Compounds/crispr for {var_text.upper()}", df_tot)
    list_batch = df_tot["batchid"]
    batch_list = []
    for batch in list_batch:
        batch_list.append("'" + batch + "'")
    conn2 = psycopg2.connect(host="192.168.2.131", port="5432", user="arno", database="ksilink_cpds", password="12345")
    sql_profile = "select * from aggprofiles where batchid  in (" + ",".join(batch_list) + ")"
    df_prof = pd.read_sql(sql_profile, conn2)
    tab1, tab2 = st.tabs(["Profiles", "Summary"])
    tab1.write(df_prof)
    tab2.write(df_prof.describe())

    st.session_state["df"] = df_prof
# df_prof.to_csv(f"Data_for_{var_text}.csv", index=None)
# Print results.
