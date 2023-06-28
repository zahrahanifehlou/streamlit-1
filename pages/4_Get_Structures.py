import os
from urllib.parse import urlencode
from urllib.request import urlretrieve
import pandas as pd
import plotly.express as px
import psycopg2
import streamlit as st
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import tmap


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
        "Choose files with batchids if df null", accept_multiple_files=True
    )
    list_df = []
    for uploaded_file in uploaded_files:
        if ".csv" in uploaded_file.name:
            list_df.append(pd.read_csv(uploaded_file))

        if ".fth" in uploaded_file.name:
            list_df.append(pd.read_feather(uploaded_file))
    return list_df


if "df_cpds" not in st.session_state:
    list_df = upload_files()

    if len(list_df) > 0:
        df_keep = pd.concat(list_df)
        df_c = get_struc(df_keep)
        if len(df_c) > 0:
            st.session_state["df_kegg"] = df_c
else:
    df_cpd = st.session_state["df_cpds"]
    df_cpd = df_cpd[df_cpd["smile"] != "No result"].reset_index()
    if len(df_cpd) > 0:
        
        
        st.write("DATA frame",df_cpd)
        mol_list = [Chem.MolFromSmiles(smiles) for smiles in df_cpd["smile"].tolist()]
        
        # plot similarity-----------------------------------------------------------------------------------
        df_c = get_struc(mol_list)
        col1, col2 = st.columns(2)
        with col1:
            st.write("Structure Similarities", df_c)
        with col2:
            fig = px.imshow(df_c)
            st.plotly_chart(fig)

        # plot tmap-----------------------------------------------------------------------------------
        bits = 1024
        df_tmap = df_cpd.copy()
        morgan_list = [
            AllChem.GetMorganFingerprintAsBitVect(
                mol, useChirality=True, radius=3, nBits=bits
            )
            for mol in mol_list
        ]
        morgan_list = [tmap.VectorUchar(list(fp)) for fp in morgan_list]

        enc = tmap.Minhash(bits, seed=42)
        lf_morgan = tmap.LSHForest(bits)
        lf_morgan.batch_add(enc.batch_from_binary_array(morgan_list))
        lf_morgan.index()
        x, y, s, t, _ = tmap.layout_from_lsh_forest(lf_morgan)

        df_tmap["x"] = x
        df_tmap["y"] = y

        fig = px.scatter(df_tmap, x="x", y="y", hover_data=["name", "pubchemid"])
        st.write("TMAP")
        st.plotly_chart(fig, theme="streamlit", use_container_width=True)

        # plot smile-----------------------------------------------------------------------------------
        ik = 0
        cpt = 0
        cols = st.columns(5)
        for mol in mol_list:
            if cpt == 5:
                cpt = 0
            try:
                cols[cpt].image(Draw.MolToImage(mol), df_cpd["pubchemid"][ik])
                ik = ik + 1
                cpt = cpt + 1
            except:
                st.write("http error: ")
