

import sys

import pandas as pd
import psycopg2
import streamlit as st
import polars as pl
import streamlit.components.v1 as components
from pygwalker.api.streamlit import init_streamlit_comm, get_streamlit_html

st.set_page_config(
    layout="wide",
)
sys.path.append('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/')
init_streamlit_comm()

@st.cache_resource
def get_pyg_html(df: pd.DataFrame) -> str:
    # When you need to publish your application, you need set `debug=False`,prevent other users to write your config file.
    # If you want to use feature of saving chart config, set `debug=True`
    html = get_streamlit_html(df,spec="./chart_meta_0.json", themeKey='vega',use_kernel_calc=True, debug=False)
    return html

@st.cache_resource
def sql_df(sql_str, _conn):
    df_d =pl.read_database_uri(query=sql_str, uri=_conn)

    return df_d.to_pandas()
# from streamlib import sql_df
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])

# conn_meta = init_connection()
conn_meta = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"
def convert_df(df):
       return df.to_csv(index=True).encode('utf-8')

profile_conn = "postgres://arno:12345@192.168.2.131:5432/ksilink_cpds"


st.header("In house Dataset",divider='rainbow')


sql_proj='select project from projectsprofile group by project'
res_proj=sql_df(sql_proj, profile_conn)


# st.write(res_proj)
project_name = st.selectbox("Select Project name",res_proj,help='DM1')
sql_profile = f"SELECT assay from projectsprofile WHERE projectsprofile.project='{project_name}' group by assay"
res_assay=sql_df(sql_profile, profile_conn)
# list_assay = res_assay.assay.unique()
sel_proj = st.selectbox("Select Assay name",res_assay)


@st.cache_data
def get_data_once(sel_projet):
    sql_assay=f"SELECT * from projectsprofile WHERE projectsprofile.assay='{sel_projet}'"
    df_pro=sql_df(sql_assay, profile_conn)
    original_columns = [ "name", "batchid", "tags", "plate", "well","project","assay","concentration"]  
 
    pivot_df = df_pro.pivot( index=original_columns, columns='feature', values='value').reset_index()
    return pivot_df


pivot_df=get_data_once(sel_proj)
pivot_df=pivot_df.apply(pd.to_numeric, errors='ignore')

components.html(get_pyg_html(pivot_df), height=1000, scrolling=True)
