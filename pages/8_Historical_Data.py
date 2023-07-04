

import pandas as pd
import psycopg2
import streamlit as st

profile_conn = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksilink_cpds",
    password="12345",
)

cols1= st.columns(2)
with cols1[0]:
    to_find = st.text_input("Enter your search",help='FCCP')
with cols1[1]:
    project_name = st.text_input("Enter Project name",help='DM1')
if project_name!="":
    sql_profile = f"SELECT * from projectsprofile WHERE projectsprofile.project='{project_name}' AND ( batchid='{to_find}' OR UPPER(projectsprofile.name) LIKE UPPER('{to_find}%'))"
else:
    sql_profile = f"SELECT * from projectsprofile WHERE ( batchid='{to_find}' OR UPPER(projectsprofile.name) LIKE UPPER('{to_find}%'))"

df_pro = pd.read_sql(sql_profile, profile_conn)
st.write(df_pro)

original_columns = ['name', 'batchid', 'concentration', 'tags', 'plate','well','assay']
for g, data in df_pro.groupby('assay'):
    pivot_df = pd.pivot_table(data, index=original_columns, columns='feature', values='value').reset_index()
    tab1,tab2=st.tabs([f"Profiles in {g}", f"Summary in {g}"])
    tab1.write(pivot_df)
    tab2.write(pivot_df.describe())
    # st.write(pivot_df)


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