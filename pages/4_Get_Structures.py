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


@st.cache_resource
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])


def tanimoto_similarity(mol1, mol2):
    fp1 = AllChem.GetMorganFingerprint(mol1, 2)
    fp2 = AllChem.GetMorganFingerprint(mol2, 2)
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity

def get_struc(df):
    batch_list = df["batchid"].tolist()
    b_list = []



    df_cpd=pd.DataFrame()
    for b in batch_list:
        if "jcp2022_800" not in b:
            b_list.append("'" + b + "'")
            # b_list.append(b.upper())
    # st.write(b_list)
    if len(b_list)>0:
        sql_profile = "select * from batchs where batchid  in (" + ",".join(b_list) + ")"
        df_int = pd.read_sql(sql_profile, conn)

        b_list2 = df_int["pubchemid"].to_list()
        if len(b_list2)>0:
            bq = []
            for bs in b_list2:
                bq.append("'" + bs + "'")

            sql_batch = "select * from cpd where pubchemid in (" + ",".join(bq) + ")"
            df_cpd = pd.read_sql(sql_batch, conn)
            df_cpd.fillna('No result',inplace=True)
            df_cpd=df_cpd[df_cpd["smile"]!="No result"].reset_index()
            # Define the molecule ID
            st.write(df_cpd)
            ik = 0
            cpt = 0
            smiles_list = df_cpd['smile'].tolist()
            # Convert SMILES strings to RDKit molecules
            molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

            # Compute pairwise similarities between the N compounds
            similarity_matrix = []
            for i in range(len(molecules)):
                row = []
                for j in range(len(molecules)):
                    similarity = tanimoto_similarity(molecules[i], molecules[j])
                    row.append(similarity)
                similarity_matrix.append(row)

            sim_df = pd.DataFrame(similarity_matrix)
            fig = px.imshow(sim_df,title='Structure Similarities')
            st.plotly_chart(fig)
            # st.write('Structure Similarities', similarity_matrix)


            cols = st.columns(5)
            for id in df_cpd["pubchemid"]:
                mol_id = id
                if cpt==5:
                    cpt=0
                # Define the URL for the PubChem structure image API
                url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/PNG".format(mol_id)

                # Define the file name for the image
                file_name = "{}.png".format(mol_id)

                # Define the file path for the image
                file_path = os.path.join(os.getcwd(), file_name)

                # Download the image from the URL and save it to the file path
                urlretrieve(url, file_path)

                cols[cpt].image(file_path, df_cpd["name"][ik] + ":" + df_cpd["keggid"][ik])
                ik = ik + 1
                cpt =cpt+1
                os.remove(file_path)
    return df_cpd


def upload_files():
    uploaded_files = st.file_uploader("Choose files with batchids if df null", accept_multiple_files=True)

    # st.write(uploaded_file)
    list_df = []
    for uploaded_file in uploaded_files:
        if ".csv" in uploaded_file.name:
            list_df.append(pd.read_csv(uploaded_file))

        if ".fth" in uploaded_file.name:
            list_df.append(pd.read_feather(uploaded_file))
    return(list_df)

conn = init_connection()
if "df_keep" not in st.session_state:
    list_df=upload_files()

    if len(list_df) > 0:
        df2 = pd.concat(list_df)
        df_c = get_struc(df2)
        if len(df_c)>0:
            st.session_state["df_kegg"] = df_c
else:
    df2 = st.session_state["df_keep"]
    df_c = get_struc(df2)
    if len(df_c)>0:
        st.session_state["df_kegg"] = df_c
    else:
        st.write('No Data from pubchems')
        list_df_p = upload_files()
        if len(list_df_p) > 0:
            df3 = pd.concat(list_df_p)
            df_c = get_struc(df3)
        if len(df_c)>0:
            st.session_state["df_kegg"] = df_c
