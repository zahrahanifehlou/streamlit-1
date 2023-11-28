import streamlit as st
from streamlib import sql_df

import pandas as pd

st.set_page_config(
    layout="wide",
)


conn = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"
df_results=pd.DataFrame()



st.image("DB.png")
st.write("Example:  SELECT cpd.pubchemid, cpd.keggid, cpd.cpdname, gene.* from cpd INNER JOIN keggcpdgene ON keggcpdgene.keggid=cpd.keggid INNER JOIN gene ON gene.geneid=keggcpdgene.geneid")


sql_line = st.text_area("Enter your search",help="select * from cpd")
if len(sql_line)>0:
    st.write("your SQL is :",sql_line)
    df_results= sql_df(sql_line, conn)


   
if len(df_results)>0:
            st.write(df_results)
      
    