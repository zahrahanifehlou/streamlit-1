import warnings

# import bbknn
import numpy as np
import pandas as pd
import plotly.express as px
import psycopg2

# names = []
# targets = []
# tar_hsa = []
# pubchems =[]
# effs = []
# cpd_ids = []
# def get_compound_info(compound_id):

#     cpd_ids.append(compound_id)
#     url = f'https://rest.kegg.jp/get/{compound_id}'

#     # Make the REST API call and get the response
#     response = requests.get(url)

#     # Parse the response using BeautifulSoup
#     soup = BeautifulSoup(response.content, 'html.parser')

#     # Extract the relevant information from the parsed response
#     result = soup.get_text()
#     if 'NAME' in result:
#         names.append(result.split("NAME")[1].split("\n")[0].strip().split(' ')[0].split(';')[0])
#         print("NAME: ", result.split("NAME")[1].split("\n")[0].strip().split(' ')[0].split(';')[0])
#     else:
#         names.append('None')
        
#     if 'TARGET' in result:
#         temp_str = result.split("TARGET")[1].split("\n")[0].strip().split(' ')[0]
#         targets.append(temp_str)
#         if len(result.split("TARGET")[1].split("\n")[0].strip().split('['))>1:
#             tar_hsa.append(result.split("TARGET")[1].split("\n")[0].strip().split('[')[1].split(']')[0].split(' ')[0])
#             print("TARGET: ", temp_str)
#         else:
#             tar_hsa.append('None')
#     else :
#         print("TARGET: ", "None")
#         targets.append('None')
#         tar_hsa.append('None')
    
#     if 'PubChem:' in result:
#         pubchems.append(result.split("PubChem:")[1].split("\n")[0].strip())
#         print("PubChem: ", result.split("PubChem:")[1].split("\n")[0].strip())
#     else :
#         pubchems.append('None')
#         print("PUBCHEM: ", "None")
    
#     if 'EFFICACY' in result:
#         effs.append(result.split("EFFICACY")[1].split("\n")[0].strip())
#         print("EFFICACY: ", result.split("EFFICACY")[1].split("\n")[0].strip())
#     else :
#         effs.append('None')
#         print("EFFICACY: ", "None")

# k=KEGG()
# response=k.list('disease')

# url = f"http://rest.kegg.jp/link/hsa/disease+ds:{code_dis}"
# # st.write(url)
# response = requests.get(url)
# list_disease=[]
# list_dis_id=[]
# if response.status_code == 200:
#         disease_data = response.text.split("\n")
#         for dis in disease_data:
#             if len(dis)>1:
#                 # list_dis_id.append(dis.split("\t")[0])
#                 list_disease.append(dis.split("\t")[1])
# # else:
# #     col2.write("Bad request")

# df_link_genes=pd.DataFrame()
# df_link_genes['genes']=list_disease

# url = f"http://rest.kegg.jp/link/drug/disease+ds:{code_dis}"
# # st.write(url)
# response = requests.get(url)
# list_disease=[]
# list_dis_id=[]
# if response.status_code == 200:
#         disease_data = response.text.split("\n")
#         for dis in disease_data:
#             if len(dis)>1:
#                 # list_dis_id.append(dis.split("\t")[0])
#                 list_disease.append(dis.split("\t")[1][3:])
# # else:
# #     col1.write("Bad request")

# df_link_drug=pd.DataFrame()
# df_link_drug['drugs']=list_disease


# import scanpy as sc
# import streamlit as st
# import umap

# conn_profileDB = psycopg2.connect(
#         host="192.168.2.131",
#         port="5432",
#         user="arno",
#         database="ksilink_cpds",
#         password="12345",
#     )

# col1, col2 = st.columns(2)
# list_sources = [
#     "Amgen" ,
#     "Astrazeneca",
#      "Bayer",
#      "Broad",
#      "Eisai",
#     "Janssen",
#      "Ksilink_25",
#      "Ksilink_625",
#      "Merck",
#     "Pfizer",
#      "Servier" ,
#     "Takeda" ,
#     "CRISPER"
# ]
# def init_connection():
#     return psycopg2.connect(**st.secrets["postgres"])
# conn = init_connection()

#  where batchid like 'BR%'
# choix_source = st.selectbox("Select the Source", list_sources)

# def clean_pubchemid(pubchemid):
#     return str(pubchemid).split('.')[0]

# def sel_code(option, conn):
#     # sql = f"SELECT cpdbatchs.pubchemid, cpdbatchs.batchid, platemap.plate, platemap.ctrltype, platemap.well \
#     # FROM platemap INNER JOIN cpdbatchs ON cpdbatchs.batchid = platemap.batchid " # WHERE platemap.source = '{option}' "
#     sql = f"SELECT cpdbatchs.pubchemid, cpdbatchs.batchid, platemap.plate, platemap.ctrltype, platemap.well FROM platemap \
#     INNER JOIN cpdbatchs ON cpdbatchs.batchid = platemap.batchid and platemap.ctrltype= '{option}' "
#     df = sql_df(sql, conn)
#     df['pubchemid'] = df['pubchemid'].apply(clean_pubchemid)
#     df=df.drop_duplicates(subset=[ "plate", "well"]).reset_index(drop=True)
#     return df

# def main():
#     df = sel_code('poscon_diverse',conn)
#     batch_list = [f"'{batchid}'" for batchid in df["batchid"].unique()]
#     sql_profile = f"SELECT * FROM aggcombatprofile WHERE metabatchid IN ({','.join(batch_list)})"
#     df_prof = sql_df(sql_profile, conn_profileDB)
#     df_prof=df_prof.drop_duplicates(subset=[ "metabatchid","metasource"]).reset_index(drop=True)

#     df_prof.reset_index(drop=True, inplace=True)
#     # st.write(df_prof[["metasource","metabatchid"]])
#     # st.write("profile")
#     num_col=[col for col in df_prof if 'meta' not in col]
#     adata=sc.AnnData(df_prof[num_col],obs=df_prof[["metasource"]])
#     # adata = sc.read('pancreas.h5ad', backup_url='https://www.dropbox.com/s/qj1jlm9w10wmt0u/pancreas.h5ad?dl=1')
#     sc.tl.pca(adata,svd_solver='arpack')
#     bbknn.bbknn(adata, batch_key='metasource')

#     df_prof=sc.AnnData.to_df(adata)
#     print(df_prof.describe())

def connect():
    import requests
    from unidecode import unidecode
    import json

    url2 ="http://labcollector.ksilink.int/webservice/index.php?v=2&module=strains"
    key="74a7d185e35ca33dc082f3e1605b914d1c6fc1c1add3ef4f96e6a284952199f2"

    # client = requests.Session()
    # client.timeout = 900000
    # request = requests.Request()
    # request.timeout = 900000
    # request.method = 'GET'
    # request.headers["X-LC-APP-Auth"] = key
    # request.headers["Accept"] = "text/xml"
    # request.headers["X-LC-APP-Charset"]=

    headers = {
        "X-LC-APP-Auth": key.encode("latin-1"),
        "Accept": "application/json"
    }
    response=requests.get(url2,headers=headers)
    # response = client.send(request.prepare(), url=url)
    r = response.text
    jj = json.loads(r)
    df = pd.json_normalize(jj['name'])
    

if __name__ == "__main__":
    # main()
    connect()