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


conn = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksi_cpds",
    password="12345",
)

@st.cache_data
def toto():
    sql_umqpemd = f"select * from umapemd where metasource='Ksilink_625'"
    df_src_emd = sql_df(sql_umqpemd, profile_conn)
    return df_src_emd
# Writes a component similar to st.write()


df=toto()
# st.write(df)
fig = px.scatter(df,x="umap1",
                y="umap2",)
selected_points = plotly_events(fig,click_event=True,hover_event=True)
if selected_points:
    batch = df.iloc[selected_points[0]['pointIndex']].metabatchid
    sql_test = f"select * from platemap where batchid='{batch}' and source='Ksilink_625'"
    df_test_emd = sql_df(sql_test, conn)
    st.write(sql_test)
    st.write(df_test_emd)