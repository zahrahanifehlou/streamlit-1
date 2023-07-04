import os
from urllib.parse import urlencode
from urllib.request import urlretrieve

import pandas as pd
import plotly.express as px
import psycopg2
import pubchempy as pcp
import streamlit as st
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw


def tanimoto_similarity(mol1, mol2):
    fp1 = AllChem.GetMorganFingerprint(mol1, 2)
    fp2 = AllChem.GetMorganFingerprint(mol2, 2)
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity


def get_struc(mol_list):
    similarity_matrix = []
    for i in range(len(mol_list)):
        row = []
        for j in range(len(mol_list)):
            similarity = tanimoto_similarity(mol_list[i], mol_list[j])
            row.append(similarity)
        similarity_matrix.append(row)
    sim_df = pd.DataFrame(similarity_matrix)
    return sim_df


def upload_files():
    uploaded_files = st.file_uploader(
        "Choose files with pubchemid if df null", accept_multiple_files=True
    )
    list_df = []
    for uploaded_file in uploaded_files:
        if ".csv" in uploaded_file.name:
            list_df.append(pd.read_csv(uploaded_file))

        if ".fth" in uploaded_file.name:
            list_df.append(pd.read_feather(uploaded_file))
    return list_df

conn_profileDB = psycopg2.connect(
        host="192.168.2.131",
        port="5432",
        user="arno",
        database="ksilink_cpds",
        password="12345",
    )
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])


        
#------------------------------------------------------------------------------------------
if "df_tmap" not in st.session_state:
    sql = "SELECT * FROM tmap"
    df_tmap = pd.read_sql(sql, conn_profileDB).reset_index(drop=True)
    conn_profileDB.close()
else:
    df_tmap = st.session_state["df_tmap"]
color_col = st.radio(
        "select color",
        ("source", "genetarget", "efficacy", "disname"),
        horizontal=True,
    )
st.write("TMAP plt of all compounds")
fig = px.scatter(df_tmap, x="x", y="y", hover_data=['pubchemid', 'name','genetarget', 'efficacy', 'disname', 'keggid'],opacity=0.3,color=color_col,width=1500,height=1000 ) # color_continuous_scale='RdBu_r',
st.plotly_chart(fig, theme="streamlit", use_container_width=True)
st.write("\n")
#-----------------------------------------------------------------------------------------------------------------------
st.write("Connect DB First or load a file")
df_cpds=pd.DataFrame()
if "df_cpds"  in st.session_state:
    df_cpds = st.session_state["df_cpds"]
    df_cpds = df_cpds[df_cpds["smile"] != "No result"].reset_index()
    df_cpds=df_cpds.drop_duplicates(subset=["pubchemid"])
    df_cpds=df_cpds[df_cpds["pubchemid"].notna()]
cols1= st.columns(3)
with cols1[0]:
    if len(df_cpds)>0:
        st.write("results from previous step")
        st.write(df_cpds)
    
with cols1[1]:
    list_pubchem = st.text_area("Enter your search",help='pubchemid separated by enter')
    if len(list_pubchem)>0:
        list_pubchem=df_cpds["pubchemid"].to_list()
        df_cpds=pd.DataFrame()
        df_cpds["pubchemid"]=list_pubchem
        smile_list=[]
        for pubchemid  in list_pubchem:
                try:
                    s=pcp.Compound.from_cid(pubchemid).canonical_smiles
                    smile_list.append(s)
                except:
                    smile_list.append(None)
        df_cpds["smile"]=smile_list
with cols1[2]:
    list_df2 = upload_files()
    if len(list_df2)>0:
        df_cpds=pd.concat(list_df2)
        df_cpds['pubchemid'] = df_cpds['pubchemid'].astype(str)
        df_cpds['pubchemid'] = df_cpds['pubchemid'].str.split('.').str[0]
        smile_list=[]
        for pubchemid  in df_cpds["pubchemid"]:
            try:
                s=pcp.Compound.from_cid(pubchemid).canonical_smiles
                smile_list.append(s)
            except:
                smile_list.append(None)
        df_cpds["smile"]=smile_list
        df_cpds=df_cpds[df_cpds["smile"].notna()]
        list_pubchem=df_cpds["pubchemid"]
#-get smiles and mol info -------------------------------------------------------- 
if len(df_cpds)  >0:
    mol_list=[]
    for smiles in df_cpds['smile']:
        try:
            mol_list.append(Chem.MolFromSmiles(smiles))
        except:
            st.write(smiles)

    mainTabs=st.tabs(["plot structure ","find similar in structure "])
    with mainTabs[0]:
        df_str=df_cpds.copy()
        df_tmap["selected compounds"]="others"
        df_tmap.loc[df_tmap['pubchemid'].isin(df_str.pubchemid), 'selected compounds'] = "selected compounds"
        st.write("TMAP plt of selected compounds")
        fig = px.scatter(df_tmap, x="x", y="y", hover_data=['pubchemid', 'name','genetarget', 'efficacy', 'disname', 'keggid'],color_discrete_sequence=["blue", "red" ],color="selected compounds",width=800,height=1000 )
        st.plotly_chart(fig, theme="streamlit", use_container_width=True)
        st.write("\n")

                # plot similarity-----------------------------------------------------------------------------------
        st.write("DATA frame",df_str)
        
        df_c = get_struc(mol_list)

        name_list = [f"`{t}`" for t in df_str["pubchemid"]]
        df_c = df_c.rename(index=dict(enumerate(name_list, 0)), columns=dict(enumerate(name_list, 0)))
        fig = px.imshow(df_c, color_continuous_scale='RdBu_r')
        st.plotly_chart(fig)

        # plot smile-----------------------------------------------------------------------------------

        ik = 0
        cpt = 0
        cols = st.columns(5)
        for mol in mol_list:
            if cpt == 5:
                cpt = 0
            try:
                cols[cpt].image(Draw.MolToImage(mol), df_str["pubchemid"][ik])
                ik = ik + 1
                cpt = cpt + 1
            except:
                st.write("")
  

    with mainTabs[1]:
        conn = init_connection()
        sql_cpd = f"select cpd.pubchemid, cpd.name, cpd.smile from cpd  "
        jump_df=pd.read_sql(sql_cpd,conn)
        conn.close()
        jump_df=jump_df[jump_df["smile"].notna()]
        mol_list_jump=[]
        find_list=[]
        bits=1024
        for smiles in jump_df['smile']:
            find=1
            try:
                mol= Chem.MolFromSmiles(smiles)
                morgan=AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=3, nBits=bits)
                mol_list_jump.append(morgan)
                # mol= mol=Chem.MolFromSmiles(smiles)
                # fp=Chem.RDKFingerprint(mol)
                # mol_list_jump.append(fp)
            except:
                find=0

            find_list.append(find)

        jump_df["find"]=find_list
        jump_df=jump_df[jump_df["find"]!=0]

        search_tbl=pd.DataFrame() 
        df_sim_result=pd.DataFrame()
        df_sim_result["pubchemid"]=jump_df.pubchemid
        
        for pubchemid  in list_pubchem:
            try:
                mol=Chem.MolFromSmiles(pcp.Compound.from_cid(pubchemid).canonical_smiles)
                fp=AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=3, nBits=bits)
                similarities = DataStructs.BulkTanimotoSimilarity(fp, mol_list_jump)  
                df_sim_result[pubchemid]=similarities 
            except:
                print("cant find this pubchemid", pubchemid)
        st.write(df_sim_result)       


 




