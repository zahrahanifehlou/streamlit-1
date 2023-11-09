

import sys

import pandas as pd
import psycopg2
import streamlit as st

sys.path.append('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/')

@st.cache_resource
def sql_df(sql_str, _conn):
    cur = _conn.cursor()
    cur.execute(sql_str)
    rows = cur.fetchall()
    df_d = pd.DataFrame(rows, columns=[desc[0] for desc in cur.description])

    return df_d
# from streamlib import sql_df


def convert_df(df):
       return df.to_csv(index=False).encode('utf-8')


profile_conn = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksilink_cpds",
    password="12345",
)
st.set_page_config(
    layout="wide",
)
st.header("In house Dataset",divider='rainbow')
# sql_line = st.text_area("Enter your search",help="select * from cpd")
# if len(sql_line)>0:
#     st.write("your SQL is :",sql_line)
#     df_results= sql_df(sql_line, profile_conn)


   
# if len(df_results)>0:
#             st.write(df_results)
#             st.download_button(
#             label="Save results",
#             data=convert_df(df_results),
#             file_name="df_results.csv",
#             mime="csv",
#     )

sql_proj='select project from projectsprofile group by project'
res_proj=sql_df(sql_proj, profile_conn)

project_name = st.selectbox("Select Project name",res_proj,help='DM1')
sql_profile = f"SELECT assay from projectsprofile WHERE projectsprofile.project='{project_name}' group by assay"
res_assay=sql_df(sql_profile, profile_conn)
# list_assay = res_assay.assay.unique()
sel_proj = st.selectbox("Select Assay name",res_assay)


@st.cache_data
def get_data_once(sel_projet):
    sql_assay=f"SELECT * from projectsprofile WHERE projectsprofile.assay='{sel_projet}'"
    df_pro=sql_df(sql_assay, profile_conn)
    original_columns = ["project","assay", "name", "batchid","concentration", "tags", "plate", "well"]    
    pivot_df = pd.pivot_table(df_pro, index=original_columns, columns='feature', values='value').reset_index()
    return pivot_df


pivot_df=get_data_once(sel_proj)
tab1,tab2=st.tabs([f"Profiles in {sel_proj}", f"Summary in {sel_proj}"])
tab1.write(pivot_df)
st.download_button(
    label="Save",data=convert_df(pivot_df),file_name=f"{sel_proj}.csv",mime='csv',)
tab2.write(pivot_df.describe())
list_cols=pivot_df.columns
numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
cols_alpha = pivot_df.select_dtypes(exclude=numerics).columns
cols = st.columns(3)
sel_col = cols[0].selectbox('select column:', list_cols)
on = st.sidebar.toggle('Drop Duplicates')
if sel_col in cols_alpha:
    sel_sign = cols[1].selectbox('select:', '==')
    str_data=pivot_df[sel_col]
    # st.write(str_data)
    sel_value= cols[2].selectbox('select value',str_data)
    df_sel = pivot_df[pivot_df[sel_col] == sel_value]
else:
    sel_sign = cols[1].selectbox('select sign:', ['<','>','=='])
    val_data_min = pivot_df[sel_col].min()
    val_data_max =  pivot_df[sel_col].max()
    val_data_med =  pivot_df[sel_col].median()
    sel_value = cols[2].slider('Value',val_data_min,val_data_max,val_data_med )
    if sel_sign == '<':
        df_sel = pivot_df[pivot_df[sel_col] < sel_value]
    elif sel_sign == '>':
        df_sel = pivot_df[pivot_df[sel_col] > sel_value]
    elif sel_sign == '==':
        df_sel = pivot_df[pivot_df[sel_col] == sel_value]
if on:
    df_agg=df_sel.drop_duplicates(subset=['batchid','concentration']).reset_index()
else:
    df_agg=df_sel        
s2 = df_agg.style.highlight_min(subset=df_agg.select_dtypes(include=numerics).columns,props='color:white;background-color:darkred',axis=0)
s3 = s2.highlight_max(subset=df_agg.select_dtypes(include=numerics).columns,props='color:white;background-color:darkblue',axis=0)
st.dataframe(s3)
# ,column_config={sel_col:st.column_config.BarChartColumn("PlotSel",y_min=val_data_min,y_max=val_data_max),}




st.download_button(
        label="Save",data=convert_df(df_agg),file_name=f"{sel_col}+{sel_sign}+{sel_value}.csv",mime='csv',)
# cols1= st.columns(2)
# with cols1[0]:
#     to_find = st.text_input("Enter your search (Batchid)",help='FCCP')
# with cols1[1]:

    
    

    
    
# if project_name!="":
#     sql_profile = f"SELECT * from projectsprofile WHERE projectsprofile.project='{project_name}' \
#     AND ( batchid='{to_find}' OR UPPER(projectsprofile.name) LIKE UPPER('{to_find}%'))"
# else:
#     sql_profile = f"SELECT * from projectsprofile WHERE ( batchid='{to_find}'\
#     OR UPPER(projectsprofile.name) LIKE UPPER('{to_find}%'))"

# df_pro = sql_df(sql_profile, profile_conn)
# st.write(df_pro.assay.unique())
# st.write("--------------------------------------")

# disp = st.sidebar.toggle("Display Data")
# if len(df_pro)>0 and (to_find or project_name):


#     original_columns = ["project","assay", "name", "batchid","concentration", "tags", "plate", "well"]
#     for g, data in df_pro.groupby('assay'):
#         pivot_df = pd.pivot_table(data, index=original_columns, columns='feature', values='value').reset_index()
#         tab1,tab2=st.tabs([f"Profiles in {g}", f"Summary in {g}"])
#         if disp:
#             tab1.write(pivot_df)
#             st.download_button(
#                 label="Save",data=convert_df(pivot_df),file_name=f"{g}.csv",mime='csv',)
#             tab2.write(pivot_df.describe())


# profile_conn.close()

# for g in df_pro_gr.groups.keys:

#     pivot_df = pd.pivot_table(df_pro, index=original_columns, columns='feature', values='value').reset_index()

# # Print the pivoted DataFrame
# st.write(pivot_df)

# import datetime
# import os
# from glob import iglob
# import plotly.express as px
# import pymongo
# from pymongo import MongoClient
# from torch import concat
# CONNECTION_STRING = "mongodb://192.168.2.127:27017/"

# client = MongoClient(CONNECTION_STRING)
# db=client.tags.tags
# db2=client.tags.compounds
# @st.cache_data
# def getPlatesFromCompounds(query):
#     matcher=[
#         {
#             '$project': {
#                 'data': {
#                     '$objectToArray': '$map'
#                 },
#                 'plate': '$plate',
#                 'plateAcq': '$plateAcq',
#                 'proj': "$meta.project"
#             }
#         }, {
#             '$set': {
#                 'data': {
#                     '$filter': {
#                         'input': '$data',
#                         'as': 'p',
#                         'cond': {
#                             '$regexMatch': {
#                                 'input': '$$p.k',
#                                 'regex': query,
#                                 'options': 'i'
#                             }
#                         }
#                     }
#                 }
#             }
#         }, {
#             '$unwind': {
#                 'path': '$data',
#                 'preserveNullAndEmptyArrays': False
#             }
#         }
#     ]
#     col=db.aggregate(matcher)
#     col2 = db2.aggregate(matcher)
#     return pd.concat([pd.DataFrame(list(x)) for x in [col,col2]],ignore_index=True)

# @st.cache_data
# def getPluginsDataset(df):
#     commits={}
#     for (project, plate), data in df.groupby(["proj", 'plate']):
#         if project not in commits:
#             commits[project]=set()
#         basepath=f"/mnt/shares/L/PROJECTS/{project}/Checkout_Results/*/ag{plate}.fth"
#         for commit in iglob(basepath):
#             # st.write(commit)
#             commit = commit.replace('\\', '/').split('/')[7]

#             if commit != "BirdView":
#                 commits[project].add(commit)
#     d=[]
#     for (project, plate), data in df.groupby(["proj", 'plate']):
#         # st.write('PP', project)
#         # st.write('PP2', plate)
#         for commit in commits[project]:
#             # if commit =='efficientnet_run2_p1-4_positiveCpds':
#             #     st.write('Commit', commit)
#             basepath=f"/mnt/shares/L/PROJECTS/{project}/Checkout_Results/{commit}/*.json"
#             for x in iglob(basepath):
#                 # st.write('X', x)
#                 if os.path.exists(f"/mnt/shares/L/PROJECTS/{project}/Checkout_Results/{commit}/{plate}.fth"):
#                     jso=x.replace('\\', '/').split('/')[-1].replace('.json', '').split('_')
#                     if jso[-1].isnumeric() and jso[-2].isnumeric():
#                         d.append((project, "_".join(jso[:-2]), commit, plate, datetime.datetime.strptime(jso[-2],'%Y%m%d').date(), jso[-1]))
#                         break
#     return pd.DataFrame(d, columns=['Project', 'Plugin', 'Commit', 'Plate', 'Date', 'time'])


# if "choice_plugin" not in st.session_state:
#      st.session_state.choice_plugin = False

# var_text = st.text_input("Enter your search")
# if var_text:
#     st.cache_resource.clear()


#     df=getPlatesFromCompounds(var_text+'.*')
#     st.write(df)

#     experiments=getPluginsDataset(df)
#     st.write(experiments)

#     list_proj=experiments['Project'].unique().tolist()
#     Proj = st.selectbox('Chose Project', list_proj)
#     df_plug=experiments[experiments['Project']==Proj].reset_index(drop=True)
#     list_plugin = df_plug['Plugin'].unique().tolist()

#     # st.balloons()
#     Plug = st.selectbox('Chose Plugin', list_plugin)
#     st.session_state.choice_plugin=True
#     df_commit=df_plug[df_plug['Plugin']==Plug].reset_index(drop=True)

#     st.write(df_commit)

#     d = st.date_input(
#         "Chose a date",
#         value=[df_commit.Date.min(),df_commit.Date.max()],min_value=df_commit.Date.min())
#     # st.write('Chosen Dates:', d)

#     # df_commit['Date']= pd.to_datetime(df_commit['Date'], format='%Y%m%d').dt.date
#     if len(d)==2:
#         df_date=df_commit[df_commit['Date'].between(*d)].reset_index(drop=True)

#         st.write('Data for these dates', df_date)

#         list_commit = df_date['Commit'].unique()
#         list_df=[]
#         for commit in list_commit:
#             pat = f"/mnt/shares/L/PROJECTS/{Proj}/Checkout_Results/{commit}/"
#             list_ag=os.listdir(pat)
#             for ag in list_ag:
#                 if 'ag' in ag and '.fth' in ag:
#                     # st.write(pat+ag)
#                     df_temp= pd.read_feather(pat+ag)
#                     df_temp=df_temp[df_temp['tags'].str.contains(var_text)].reset_index(drop=True)
#                     list_df.append(df_temp)
#         df_all = pd.concat(list_df)
#         st.write('df_all', df_all)

#         # list_sel = st.multiselect('Select Commits',list_commit)

#         df_group=df_all.groupby('tags').median(numeric_only=True).reset_index()
#         filter_col = [
#         col for col in df_group.columns if "tags" not in col
#     ]
#         normalized_df=(df_group[filter_col]-df_group[filter_col].min())/(df_group[filter_col].max()-df_group[filter_col].min())
#         normalized_df['tags']=df_group['tags']
#         st.write("df_group", df_group)
#         # st.write('df/norm',len(df_group.columns),len(normalized_df.columns))
#         cpd_names = normalized_df.tags.values
#         df_plt = normalized_df.set_index("tags")
#         df_plt = df_plt[filter_col].T
#         # st.write('df/norm',len(df_plt.columns))
#         fig2 = px.line(df_plt, x=filter_col, y=cpd_names, width=1400, height=1000)
#         st.plotly_chart(fig2, theme="streamlit", use_container_width=True)

#         fig3=px.box(normalized_df,x='tags',y='Avg nuclei area')
#         st.plotly_chart(fig3, theme="streamlit", use_container_width=True)
# # list_df=[]
# # for com in df_commit['Commit']:
# #     pd.read_feather()