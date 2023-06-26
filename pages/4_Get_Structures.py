# Import necessary libraries
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


@st.cache_resource
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])


def tanimoto_similarity(mol1, mol2):
    fp1 = AllChem.GetMorganFingerprint(mol1, 2)
    fp2 = AllChem.GetMorganFingerprint(mol2, 2)
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity


def get_pubchemPNG(df_cpd):
    mol_list = df_cpd["mol"].tolist()
    ik = 0
    cpt = 0
    cols = st.columns(5)
    for mol in mol_list:
        if cpt == 5:
            cpt = 0
        try:
            cols[cpt].image(
                Draw.MolToImage(mol), df_cpd["name"][ik] + ":" + df_cpd["pubchemid"][ik]
            )
            ik = ik + 1
            cpt = cpt + 1
        except:
            st.write("http error: ")


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


conn = init_connection()
if "df_cpds" not in st.session_state:
    list_df = upload_files()

    if len(list_df) > 0:
        df_keep = pd.concat(list_df)
        df_c = get_struc(df_keep)
        if len(df_c) > 0:
            st.session_state["df_kegg"] = df_c
else:
    df_cpd = st.session_state["df_cpds"]
    if len(df_cpd) > 0:
        st.write(df_cpd)
        smiles_list = df_cpd["smile"].tolist()

        mol_list = [Chem.MolFromSmiles(smiles) for smiles in df_cpd["smile"].tolist()]
        df_cpd["mol"] = mol_list
        df_c = get_struc(df_cpd["mol"])
        fig = px.imshow(df_c, title="Structure Similarities")
        st.plotly_chart(fig)
        st.write("Structure Similarities", df_c)
        get_pubchemPNG(df_cpd)

        if len(df_c) > 0:
            st.session_state["df_kegg"] = df_cpd
        else:
            st.write("No Data from pubchems")
            list_df_p = upload_files()
            if len(list_df_p) > 0:
                df3 = pd.concat(list_df_p)
                df_c = get_struc(df3)
            if len(df_c) > 0:
                st.session_state["df_kegg"] = df_c
