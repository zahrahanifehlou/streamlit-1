import streamlit as st
from streamlib import sql_df

import pandas as pd

st.set_page_config(
    layout="wide",
)

profile_conn = "postgres://arno:12345@192.168.2.131:5432/ksilink_cpds"
conn = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"
df_results=pd.DataFrame()
df_profiles=pd.DataFrame()

tab1, tab2 = st.tabs(["cpds","profile"])
with tab1:
    st.image("DB.png")
    
    st.write("Example:  SELECT cpd.pubchemid, cpd.keggid, cpd.cpdname, gene.* from cpd INNER JOIN keggcpdgene ON keggcpdgene.keggid=cpd.keggid INNER JOIN gene ON gene.geneid=keggcpdgene.geneid")

    sql_line = st.text_area("Enter your search",help="select * from cpd")
    if len(sql_line)>0:
        st.write("your SQL is :",sql_line)
        df_results= sql_df(sql_line, conn)


    
    if len(df_results)>0:
                st.write(df_results)
        
with tab2:
    st.image("profile_db.png")
    st.write("for profiles")
    st.write("Example for projects:  SELECT * from projectsprofile where project='PD'")
    st.write("Example for jump and crispr:  SELECT * from aggcombatprofile where metasource='CRISPER'")
    sql_line = st.text_area("Enter your search",help="select * from aggcombatprofile")
    if len(sql_line)>0:
        st.write("your SQL is :",sql_line)
        df_profiles= sql_df(sql_line, profile_conn)

    if len(df_profiles)>0:
                st.write(df_profiles)
    

    