import streamlit as st
from streamlit_plotly_events import plotly_events
import numpy as np
import plotly.express as px
from streamlib import sql_df
import psycopg2
import pandas as pd
from sklearn.cluster import KMeans
from streamlib import convert_df
st.set_page_config(
    layout="wide",
)

def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])
conn = init_connection()
df_results=pd.DataFrame()



st.image("DB.png")
st.write("Example:  SELECT cpd.pubchemid, cpd.keggid, cpd.cpdname, gene.* from cpd INNER JOIN keggcpdgene ON keggcpdgene.keggid=cpd.keggid INNER JOIN gene ON gene.geneid=keggcpdgene.geneid")


sql_line = st.text_area("Enter your search",help="select * from cpd")
if len(sql_line)>0:
    st.write("your SQL is :",sql_line)
    df_results= sql_df(sql_line, conn)


   
if len(df_results)>0:
            st.write(df_results)
            st.download_button(
            label="Save results",
            data=convert_df(df_results),
            file_name="df_results.csv",
            mime="csv",
    )
    