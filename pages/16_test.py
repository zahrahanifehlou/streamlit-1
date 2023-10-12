import streamlit as st
from streamlit_plotly_events import plotly_events
import numpy as np
import plotly.express as px
from streamlib import sql_df
import psycopg2

profile_conn = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksilink_cpds",
    password="12345",
)

@st.cache_data
def toto():
    sql_umqpemd = f"select * from umapemd where metasource='CRISPER'"
    df_src_emd = sql_df(sql_umqpemd, profile_conn)
    return df_src_emd
# Writes a component similar to st.write()


df=toto()
fig = px.scatter(df,x="umap1",
                y="umap2",)
selected_points = plotly_events(fig,click_event=True)
st.write(df.iloc[selected_points[0]['pointIndex']])