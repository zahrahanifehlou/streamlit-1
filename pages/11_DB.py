import os
import sys
from os import listdir
from pathlib import Path

import pandas as pd
import psycopg2
import requests
import streamlit as st

sys.path.append('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/')
from streamlib import sql_df

conn_meta = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksi_cpds",
    password="12345",
)

def to_str(pubchemid):
    return str(pubchemid)
# @st.cache_data
def get_cpds():
    # sql_kegg = "select cpdpath.pathid,keggcpdgene.geneid,gene.symbol, keggcpd.*, cpdbatchs.* from keggcpd\
    #             inner join keggcpdgene on keggcpd.keggid=keggcpdgene.keggid\
    #             inner join cpd on cpd.keggid=keggcpd.keggid \
    #             inner join cpdpath on cpdpath.pubchemid=cpd.pubchemid \
    #             inner join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid \
    #             inner join gene on keggcpdgene.geneid=gene.geneid"

    sql_kegg = "select keggcpdgene.geneid,gene.symbol, keggcpd.* from keggcpd\
                inner join keggcpdgene on keggcpd.keggid=keggcpdgene.keggid\
                inner join gene on keggcpdgene.geneid=gene.geneid"

    df_drug_meta = sql_df(sql_kegg, conn_meta)
    df_drug_meta = df_drug_meta.loc[:, ~df_drug_meta.columns.duplicated()].copy()
    # df_drug_meta = df_drug_meta.drop_duplicates(subset=["keggid", "source"]).reset_index(drop=True)

    df_drug_meta.dropna(subset="geneid", axis=0, inplace=True)
    return df_drug_meta

from pymongo import MongoClient

CONNECTION_STRING = "mongodb://192.168.2.127:27017/"

client = MongoClient(CONNECTION_STRING)

def get_plate(plate_names, collection='tags'):
    db=client['tags'][collection]

    if type(plate_names) is str:
        plate_names=[plate_names]

    res=[]
    for pl in plate_names:
        plates=[x for x  in db.find({'plate': pl})]
        if (len(plates) == 0) and collection == "tags":
            raise ValueError(f"plate map not found. please add plate map for the plate in checkout {pl}")
        res += plates
    return res


url2 ="http://labcollector.ksilink.int/webservice/index.php?v=2&module=strains"
key="74a7d185e35ca33dc082f3e1605b914d1c6fc1c1add3ef4f96e6a284952199f2"

headers = {
    "X-LC-APP-Auth": key.encode("latin-1"),
    "Accept": "application/json"
}
response=requests.get(url2,headers=headers)

r = response.json()

df = pd.DataFrame(r)
# st.write(df)
df.dropna(subset='project',inplace=True)
st.write('## From LabCollector:')
proj = st.selectbox('Select Project:',df['project'].unique())

df_proj=df[df["project"]==proj].reset_index(drop=True)

cells = st.multiselect(
    'Select Cell lines:', df_proj['name'].unique())
pos_ctrl=''
neg_ctrl=''
if len(cells)>1:
    col_pos,col_neg=st.columns(2)
    pos_ctrl = col_pos.selectbox('Positive Control',cells)
    neg_ctrl = col_neg.selectbox('Negative Control',cells)

# df_cells=df_proj[df_proj['name'].isin(cells)]

###################### SO MESSY HERE!!!!!!!!!!!!!!!!!!!!##############################################
cpd_inhib=''
df_inhib_sel=pd.DataFrame()
rad = st.sidebar.radio('## Any perturbations as controls?',['yes','no'])
# eb = st.empty()
if rad=='yes':
    orig=st.sidebar.radio('Registered in labcollector ?', ['yes','no'])
    if orig=='yes':
        url_cpd="http://labcollector.ksilink.int/webservice/index.php?v=2&module=chemicals"
        response=requests.get(url_cpd,headers=headers)
        cpds = response.json()
        df_cpds = pd.DataFrame(cpds)
        # col_name,col_selref=st.columns(2)
        st.multiselect('Chose Cpd Names',df_cpds['name'].unique())
        # st.write('Labc chemicals',df_cpds.head())
    if orig=='no':
        df_cpds_inhib = get_cpds()
        # st.write(df_cpds_inhib)
        inhib = st.sidebar.radio('Looking for  specific inhibitors?', ['yes','no'],key='inhib')
        if inhib=='yes':
            var_text = st.text_area("Enter your search",help='Name or ID separated by enter')
            var_t=var_text.split('\n')
            var_t = [t.strip().upper() for t in var_t]
            # dr_inhib = st.multiselect('select target(s):', df_cpds_inhib['symbol'].unique())
            # dr_inhib = [t.strip().upper() for t in dr_inhib]
            df_inhib_sel=df_cpds_inhib[df_cpds_inhib['symbol'].isin(var_t)]
            df_inhib_sel['ext_name']=df_inhib_sel['name']+'_'+df_inhib_sel['symbol']
            cpd_inhib=st.multiselect(f'select cpds targeting :{var_t}:',df_inhib_sel['ext_name'].unique())
            if cpd_inhib:
                st.write(df_inhib_sel[df_inhib_sel['ext_name'].isin(cpd_inhib)].reset_index(drop=True))

            # st.write(dr_inhib)
        else:
            df_cpds_inhib['ext_name']=df_cpds_inhib['name']+'_'+df_cpds_inhib['symbol']
            cpd_inhib = st.multiselect('Select Cpds_Target:', df_cpds_inhib['ext_name'].unique())
            df_inhib_sel=df_cpds_inhib[df_cpds_inhib['ext_name'].isin(cpd_inhib)]
            # df_inhib_sel['ext_name']=df_inhib_sel['name']+'_'+df_inhib_sel['symbol']
            st.write(df_inhib_sel)
if not df_inhib_sel.empty:
    list_pubchemid = [f"'{t}'" for t in df_inhib_sel[df_inhib_sel['ext_name'].isin(cpd_inhib)]["pubchemid"]]
    # st.write(list_pubchemid)
    if list_pubchemid:
        sql_jump=f"select * from cpdbatchs  WHERE cpdbatchs.pubchemid IN ({','.join(list_pubchemid)})"
        df_jump=sql_df(sql_jump,conn_meta)
        st.header('Existing Cpds in JUMP',divider='rainbow')
        st.write(df_jump.drop_duplicates(subset='batchid').reset_index(drop=True))
#############################################################################################################

# st.cache_data.clear()
st.header('',divider='rainbow')
# st.write(f'###  Project: {proj}')
# st.write(f'####  Cell Lines: {cells} ')
# st.write(f'####  Pertubations as ctrls: {cpd_inhib}')
st.latex(r'''\mathcal{\text{Summary}} ''')
st.latex(f'Project: {proj}')
for i in range(len(cells)):
    st.latex(f'Cells_{i}: {cells[i].split("#")[0].replace("_","+")}')
for i in range(len(cpd_inhib)):
    st.latex(f'Cpds_{i}: {cpd_inhib[i].split("_")[0]}')
st.divider()
st.header('Data', divider='rainbow')

on = st.sidebar.toggle('Activate Data Upload')
if on:
    #this is function to run
    # def get_folder_path():
    #     path = os.path.abspath('lib')
    #     p = subprocess.Popen(['python','tkDirSelector.py'], cwd=path,stdin=subprocess.PIPE,
    #                          stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    #     result,error = p.communicate()
    #     p.terminate()
    #     st.write(result)
    #     if isinstance(result,bytes):
    #         return result.decode('utf-8')
    #     if isinstance(result,str):
    #         return result
    #     else:
    #         return "error occured"
    # get_folder_path()
    try:
        paths = sorted(Path(f'/mnt/shares/L/Projects/{proj}/Checkout_Results/').iterdir(), key=os.path.getmtime, reverse=True)
    except:
        paths=Path('/mnt/shares/L/Projects/').iterdir()
    # st.write(paths)
    uploaded_files = st.multiselect("Choose result directories corresponding to this assay", paths)
    # st.write('Selected data directories to upload:')

    list_df=[]
    if uploaded_files:
        for item in uploaded_files:
            fth_file=[f for f in listdir(item) if f.startswith('ag')]
            if len(fth_file)>0:
                for fi in fth_file:
                    # st.write(item.as_posix())
                    list_df.append(pd.read_feather(item.as_posix()+'/'+fi))
        if len(list_df)>0:
            df_data =  pd.concat(list_df)
            st.write('Upload these Data to DB:',df_data)
            # crap = get_plate(df_data['Plate'])
            # # df_mongo = [pd.DataFrame(plate_data) for plate_data in crap]
            # # df_mo = pd.concat(df_mongo,axis=1)
            # df_mongo = pd.json_normalize(crap['Categories'],max_level=0)
            # st.write(df_mongo)
            # for item in crap:

        #     # # # Use json_normalize() to convert JSON to DataFrame
        #     #     # dict = json.loads(item)
        #     #     # df_mongo = pd.json_normalize(item)
        #     # # # print(df2)

        #     #     df_mongo=pd.DataFrame((item))
        #     #     st.write(df_mongo)
        else:
            st.write('No matching Data in these directories')

