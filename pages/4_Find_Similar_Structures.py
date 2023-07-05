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
st.set_page_config(
    layout="wide",
)

def plot_smile(mol_list,df_str):
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

def tanimoto_similarity(mol1, mol2):
    fp1 = AllChem.GetMorganFingerprint(mol1, 2)
    fp2 = AllChem.GetMorganFingerprint(mol2, 2)
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity

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


def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])


# Load or fetch data
def load_data():
    if "df_tmap" in st.session_state:
        df_tmap = st.session_state["df_tmap"]
    else:
        sql = "SELECT * FROM tmap"
        conn_profileDB = psycopg2.connect(
            host="192.168.2.131",
            port="5432",
            user="arno",
            database="ksilink_cpds",
            password="12345",
        )
        df_tmap = pd.read_sql(sql, conn_profileDB).reset_index(drop=True)
        conn_profileDB.close()
        st.session_state["df_tmap"] = df_tmap
    return df_tmap


# Plot TMAP
def plot_tmap(df_tmap, color_col,l ):
    
    fig = px.scatter(
        df_tmap,
        x="x",
        y="y",
        hover_data=[
            "pubchemid",
            "name",
            "genetarget",
            "efficacy",
            "disname",
            "keggid"
        ],
       color_discrete_sequence  =l ,
        opacity=0.3,
        color=color_col,
        width=1500,
        height=1000,
    )
    st.plotly_chart(fig, theme="streamlit", use_container_width=True)
    st.write("\n")


# Process compound data
def process_compound_data(df_cpds):
    df_cpds = df_cpds[df_cpds["pubchemid"].notna()]
    df_cpds = df_cpds[df_cpds["smile"] != "No result"].reset_index(drop=True)
    df_cpds = df_cpds[df_cpds["smile"].notna()]
    df_cpds = df_cpds.drop_duplicates(subset=["pubchemid"])
    
    return df_cpds

#------------------------------------------------------------------------------------------

df_tmap =load_data()
color_col = st.radio(
        "select color",
        ("source", "genetarget", "efficacy", "disname"),
        horizontal=True,
    )
st.write("TMAP plt of all compounds")
plot_tmap(df_tmap, color_col, px.colors.qualitative.D3)


#-----------------------------------------------------------------------------------------------------------------------
st.write("Connect DB First or load a file")
if "df_cpds" in st.session_state:
    df_cpds = st.session_state["df_cpds"]
else:
    df_cpds = pd.DataFrame()

cols1 = st.columns(3)

with cols1[0]:
    if not df_cpds.empty:
        st.write("Results from previous step")
        df_cpds.drop_duplicates(subset=["pubchemid"],inplace=True)
        st.write(df_cpds)

with cols1[1]:
    list_pubchem = st.text_area("Enter your search", help='pubchemid separated by enter')
    if len(list_pubchem) > 0:
        list_pubchem = list_pubchem.split('\n')
        df_cpds = pd.DataFrame({'pubchemid': list_pubchem})
  
with cols1[2]:
    list_df2 = upload_files()
    if len(list_df2) > 0:
        df_cpds=pd.concat(list_df2)
        
if len(df_cpds)>0:
    df_cpds['pubchemid'] = df_cpds['pubchemid'].astype(str).str.split('.').str[0]
    smile_list = []
    mol_list=[]
    for pubchemid in df_cpds['pubchemid'].values:
        try:
            smiles = pcp.Compound.from_cid(pubchemid).canonical_smiles
            mol_list.append(Chem.MolFromSmiles(smiles))
            smile_list.append(smiles)
        except:
            smile_list.append(None)
            st.write(f"not found {pubchemid}")
    df_cpds["smile"] = smile_list
    df_cpds = process_compound_data(df_cpds)
    list_pubchem = df_cpds["pubchemid"]
    


    mainTabs=st.tabs(["plot structure ","find similar in structure "])
    with mainTabs[0]:
                df_str = df_cpds.copy()
                df_tmap["selected compounds"]="others"
                df_tmap.loc[df_tmap['pubchemid'].isin(df_str.pubchemid), 'selected compounds'] = "selected compounds"
                st.write("TMAP plt of selected compounds")
                plot_tmap(df_tmap, "selected compounds",["red","blue"])

             # plot similarity-----------------------------------------------------------------------------------
                st.write("DATA frame",df_str)

                df_c = get_struc(mol_list)
                name_list = [f"`{t}`" for t in df_str["pubchemid"]]
                df_c = df_c.rename(index=dict(enumerate(name_list, 0)), columns=dict(enumerate(name_list, 0)))
                fig = px.imshow(df_c, color_continuous_scale='RdBu_r')
                st.plotly_chart(fig)

                # plot smile-----------------------------------------------------------------------------------
                plot_smile(mol_list,df_str)
      

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
        #17397405
       
        for pubchemid  in list_pubchem:
          
            try:
                mol=Chem.MolFromSmiles(pcp.Compound.from_cid(pubchemid).canonical_smiles)
                fp=AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=3, nBits=bits)
                similarities = DataStructs.BulkTanimotoSimilarity(fp, mol_list_jump) 
               
                sim_df=pd.DataFrame()
                name_list=[]
                value_list=[]
                mol_list2=[]
                for name, value in zip(jump_df.pubchemid, similarities):
                    
                    if value > 0.45:
                        name_list.append(name)
                        value_list.append(value)
                        smiles_2 = pcp.Compound.from_cid(name).canonical_smiles
                        mol_list2.append(Chem.MolFromSmiles(smiles_2))
               
                        #st.write(f"pubchemid: {name} , similarity:{value}")
                       
                sim_df["pubchemid"]=name_list
                sim_df["sim"]=value_list
                st.write("\n")
                
                if len(mol_list2) >0:
                    
                        s_mail = pcp.Compound.from_cid(pubchemid).canonical_smiles
                        mol_main=Chem.MolFromSmiles(s_mail)
                        st.write(f" similar compounds of {pubchemid} with  structure :")
                        st.image(Draw.MolToImage(mol_main), pubchemid)
                        plot_smile(mol_list2,sim_df) 
                        st.write(f"************************************************\n")
                
       
                        
                df_sim_result[pubchemid]=similarities 
              
            except:
                print("cant find this pubchemid", pubchemid)
        #st.write(df_sim_result)
        # for index, row in df_sim_result.iterrows():
        #     pubchemid = row['pubchemid']
        #     print(f"pubchemid: {pubchemid}")
            
        #     for column, value in row.items():
        #         if column != 'pubchemid' and value > 0.1:
        #             print(f"Column '{column}' has a value greater than 0.5: {value}")

            
    # with mainTabs[0]:
    #     df_str=df_cpds.copy()
    #     df_tmap["selected compounds"]="others"
    #     df_tmap.loc[df_tmap['pubchemid'].isin(df_str.pubchemid), 'selected compounds'] = "selected compounds"
    #     st.write("TMAP plt of selected compounds")
    #     fig = px.scatter(df_tmap, x="x", y="y", hover_data=['pubchemid', 'name','genetarget', 'efficacy', 'disname', 'keggid'],color_discrete_sequence=["blue", "red" ],color="selected compounds",width=800,height=1000 )
    #     st.plotly_chart(fig, theme="streamlit", use_container_width=True)
    #     st.write("\n")

    #     # plot similarity-----------------------------------------------------------------------------------
    #     st.write("DATA frame",df_str)
        
    #     df_c = get_struc(mol_list)
    #     name_list = [f"`{t}`" for t in df_str["pubchemid"]]
    #     df_c = df_c.rename(index=dict(enumerate(name_list, 0)), columns=dict(enumerate(name_list, 0)))
    #     fig = px.imshow(df_c, color_continuous_scale='RdBu_r')
    #     st.plotly_chart(fig)

    #     # plot smile-----------------------------------------------------------------------------------
    #     ik = 0
    #     cpt = 0
    #     cols = st.columns(5)
    #     for mol in mol_list:
    #         if cpt == 5:
    #             cpt = 0
    #         try:
    #             cols[cpt].image(Draw.MolToImage(mol), df_str["pubchemid"][ik])
    #             ik = ik + 1
    #             cpt = cpt + 1
    #         except:
    #             st.write("")
  

    # with mainTabs[1]:
    #     conn = init_connection()
    #     sql_cpd = f"select cpd.pubchemid, cpd.name, cpd.smile from cpd  "
    #     jump_df=pd.read_sql(sql_cpd,conn)
    #     conn.close()
    #     jump_df=jump_df[jump_df["smile"].notna()]
    #     mol_list_jump=[]
    #     find_list=[]
    #     bits=1024
    #     for smiles in jump_df['smile']:
    #         find=1
    #         try:
    #             mol= Chem.MolFromSmiles(smiles)
    #             morgan=AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=3, nBits=bits)
    #             mol_list_jump.append(morgan)
    #             # mol= mol=Chem.MolFromSmiles(smiles)
    #             # fp=Chem.RDKFingerprint(mol)
    #             # mol_list_jump.append(fp)
    #         except:
    #             find=0

    #         find_list.append(find)

    #     jump_df["find"]=find_list
    #     jump_df=jump_df[jump_df["find"]!=0]

    #     search_tbl=pd.DataFrame() 
    #     df_sim_result=pd.DataFrame()
    #     df_sim_result["pubchemid"]=jump_df.pubchemid
        
    #     for pubchemid  in list_pubchem:
    #         try:
    #             mol=Chem.MolFromSmiles(pcp.Compound.from_cid(pubchemid).canonical_smiles)
    #             fp=AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=3, nBits=bits)
    #             similarities = DataStructs.BulkTanimotoSimilarity(fp, mol_list_jump)  
    #             df_sim_result[pubchemid]=similarities 
    #         except:
    #             print("cant find this pubchemid", pubchemid)
    #     st.write(df_sim_result)       



 




