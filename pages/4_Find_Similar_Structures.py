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


#-----------------------------------------------------------------------------------------------------------------------
conn_profileDB = psycopg2.connect(
        host="192.168.2.131",
        port="5432",
        user="arno",
        database="ksilink_cpds",
        password="12345",
    )
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])
conn = init_connection()

if "df_tmap" not in st.session_state:
    st.write("Connect DB First or load a file")
    sql = "SELECT * FROM tmap"
    df_tmap = pd.read_sql(sql, conn_profileDB).reset_index(drop=True)


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
mainTabs=st.tabs(["plot structure ","find NOT similar in structure "])
if "df_cpds" not in st.session_state:
    list_df2 = upload_files()
    if len(list_df2)>0:
        df_file=pd.concat(list_df2)
        list_pubchem=df_file['pubchemid'].to_list()


        df_selleck_smile = df_tmap[df_tmap['source']=='SelectChem']
        df_selleck_smile=df_selleck_smile[df_selleck_smile["pubchemid"].notna()]
        st.write(df_selleck_smile)
        # list_sel_smile = df_selleck_smile['pubchemid'].to_list()
        name_pub = [f"'{t}'" for t in df_selleck_smile["pubchemid"]]
     
        sql_cpd = f"select * from cpd  WHERE cpd.pubchemid  IN ({','.join(name_pub)})"
        df_smile=pd.read_sql(sql_cpd,conn)
        mol_list1=[]
      
        for smiles in list_pubchem:
            try:
                 mol_list1.append(Chem.MolFromSmiles(pcp.Compound.from_cid(smiles).canonical_smiles))
            except:
                st.write(smiles)
        mol_list2=[]
        
        for smiles in df_smile['smile']:
            try:
                 mol_list2.append(Chem.MolFromSmiles(smiles))
            except:
                st.write(smiles)
            

        
        fps2 =[Chem.RDKFingerprint(mol) for mol in mol_list2]
        for mol in mol_list1 :
            fp=Chem.RDKFingerprint(mol)
            similarities = DataStructs.BulkTanimotoSimilarity(fp, fps2)
            for s in similarities:
                if s>0.95:
                    st.write(s)
            
        
        
        # fps = [Chem.RDKFingerprint(mol) for mol in mol_list1]
        # similarities = [DataStructs.BulkTanimotoSimilarity(fp, [Chem.RDKFingerprint(mol) for mol in mol_list2]) for fp in fps]
        # st.write(similarities)
        # for sim in similarities:
        #     if sim>0.95:
        #         st.write('sim')
    # exit(0)

    # if len(list_df) > 0:
    #     df_keep = pd.concat(list_df)
    #     df_c = get_struc(df_keep)
    #     if len(df_c) > 0:
    #         st.session_state["df_kegg"] = df_c
else:
        df_cpds = st.session_state["df_cpds"]
        df_str = df_cpds[df_cpds["smile"] != "No result"].reset_index()
        df_str=df_str.drop_duplicates(subset=["pubchemid"])
        df_str=df_str[df_str["pubchemid"].notna()]
        if len(df_str) > 0:
            with mainTabs[0] :
                df_tmap["selected compounds"]="others"
                df_tmap.loc[df_tmap['pubchemid'].isin(df_str.pubchemid), 'selected compounds'] = "selected compounds"
                st.write("TMAP plt of selected compounds")
                fig = px.scatter(df_tmap, x="x", y="y", hover_data=['pubchemid', 'name','genetarget', 'efficacy', 'disname', 'keggid'],color_discrete_sequence=["blue", "red" ],color="selected compounds",width=800,height=1000 )
                st.plotly_chart(fig, theme="streamlit", use_container_width=True)
                st.write("\n")

                # plot similarity-----------------------------------------------------------------------------------
                st.write("DATA frame",df_str)
                mol_list = [Chem.MolFromSmiles(smiles) for smiles in df_str["smile"].tolist()]
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
            with mainTabs[1] :
                cols1= st.columns(2)
                with cols1[0]:
                    pubchemid_t = st.text_input("Enter your search",help='PubChemID')
                with cols1[1]:
                    gens_text = st.text_area("search in special target",help='gene name')
                    gene_t=gens_text.split('\n')
                    gene_t = [t.strip().upper() for t in gene_t]
                if len(gene_t)>0 and gene_t[0]!="" and pubchemid_t!="":
                    gene_t = [t.strip().upper() for t in gene_t]






