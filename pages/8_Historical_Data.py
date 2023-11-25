

import sys

import pandas as pd
import psycopg2
import streamlit as st
import polars as pl

sys.path.append('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/')

# @st.cache_resource
def sql_df(sql_str, _conn):
    # cur = _conn.cursor()
    # cur.execute(sql_str)
    # rows = cur.fetchall()
    # df_d = pd.DataFrame(rows, columns=[desc[0] for desc in cur.description])
    df_d =pl.read_database_uri(query=sql_str, uri=_conn)
    # df_e=pd.DataFrame(df_d)
    return df_d.to_pandas()
# from streamlib import sql_df
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])

# conn_meta = init_connection()
conn_meta = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"
def convert_df(df):
       return df.to_csv(index=True).encode('utf-8')

profile_conn = "postgres://arno:12345@192.168.2.131:5432/ksilink_cpds"
# profile_conn = psycopg2.connect(
#     host="192.168.2.131",
#     port="5432",
#     user="arno",
#     database="ksilink_cpds",
#     password="12345",
# )
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


st.write(res_proj)
project_name = st.selectbox("Select Project name",res_proj,help='DM1')
sql_profile = f"SELECT assay from projectsprofile WHERE projectsprofile.project='{project_name}' group by assay"
res_assay=sql_df(sql_profile, profile_conn)
# list_assay = res_assay.assay.unique()
sel_proj = st.selectbox("Select Assay name",res_assay)


@st.cache_data
def get_data_once(sel_projet):
    sql_assay=f"SELECT * from projectsprofile WHERE projectsprofile.assay='{sel_projet}'"
    df_pro=sql_df(sql_assay, profile_conn)
    original_columns = [ "name", "batchid", "tags", "plate", "well","project","assay","concentration"]  
 
    pivot_df = df_pro.pivot( index=original_columns, columns='feature', values='value').reset_index()
    return pivot_df


pivot_df=get_data_once(sel_proj)
tab1,tab2=st.tabs([f"Profiles in {sel_proj}", f"Summary in {sel_proj}"])
pivot_df=pivot_df.apply(pd.to_numeric, errors='ignore')
list_cols=pivot_df.columns
numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
cols_alpha = pivot_df.select_dtypes(exclude=numerics).columns
cols_num = pivot_df.select_dtypes(include=numerics).columns

for col in cols_num:
    pivot_df[col] = pivot_df[col].apply(lambda x: "{:.2e}".format(x))
# df_temp = pivot_df[cols_alpha]
# df_temp_num = pivot_df[cols_num].format('{:.4g}')
# pivot_df=pd.concat([df_temp,df_temp_num])
tab1.dataframe(pivot_df)
# st.dataframe(df_dcp[cols_num].style.format('{:.4g}'))
tab2.write(pivot_df.describe())

# df3=df2['EC_50 Nuclei_Tot'].apply(lambda x:"{:.2e}".format(float(x)))
# df4=pd.to_numeric(df3)
# # df5 = df4.applymap('{:.2f}'.format)
# # st.write(df4.apply(lambda x:"{:.2e}".format(float(x))))
# # pd.set_option('display.float_format', '{:.2g}'.format)
# st.dataframe(df4)
for col in cols_num:
    pivot_df[col] = pivot_df[col].astype(float)

cols = st.columns(3)
sel_col = cols[0].selectbox('select column:', list_cols,index=2)
on = st.sidebar.toggle('Drop Duplicates')
if sel_col in cols_alpha:
    sel_sign = cols[1].selectbox('select:', '==')
    str_data=pivot_df[sel_col].unique()
    # st.write(str_data)
    sel_value= cols[2].selectbox('select value',str_data)
    df_sel = pivot_df[pivot_df[sel_col] == sel_value]
else:
    sel_sign = cols[1].selectbox('select sign:', ['<','>','=='])
    val_data_min = pivot_df[sel_col].min()
    val_data_max =  pivot_df[sel_col].max()
    val_data_med =  pivot_df[sel_col].median()
    val_step=(val_data_max-val_data_min)/10
    sel_value = cols[2].slider('Value',val_data_min,val_data_max,val_data_med,step=val_step,format="%7.2e" )
    if sel_sign == '<':
        df_sel = pivot_df[pivot_df[sel_col] < sel_value]
    elif sel_sign == '>':
        df_sel = pivot_df[pivot_df[sel_col] > sel_value]
    elif sel_sign == '==':
        df_sel = pivot_df[pivot_df[sel_col] == sel_value]
if on:
    df_agg=df_sel.drop_duplicates(subset=['batchid','concentration']).reset_index()
    df_agg.set_index(['name','concentration','batchid'],inplace=True)   
else:
    df_agg=df_sel
# st.write(df_agg)  
s2 = df_agg.style.highlight_min(subset=df_agg.select_dtypes(include=numerics).columns,props='color:white;background-color:darkred',axis=0)
s3 = s2.highlight_max(subset=df_agg.select_dtypes(include=numerics).columns,props='color:white;background-color:darkblue',axis=0)
s4 =s3.set_sticky(axis="index")
# pd.set_option("styler.render.max_elements", 3000000)
st.dataframe(s4)
# ,column_config={sel_col:st.column_config.BarChartColumn("PlotSel",y_min=val_data_min,y_max=val_data_max),}
df_state=df_agg.reset_index()
b_list2 = df_state['batchid'].astype(str).to_list()
   
 
bq = []
for bs in b_list2:
    bq.append("'" + bs.strip().upper() + "'")


# sql_cpd="select cpd.*, cpdbatchs.batchid, keggcpdgene.geneid from cpd \
#         INNER join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid \
#         INNER join keggcpdgene on cpd.keggid=keggcpdgene.keggid \
#         "

# st.write(sql_cpd)         
# df_cpd_meta = sql_df(sql_cpd, conn_meta)
# df_cpd_meta = df_cpd_meta.loc[:, ~df_cpd_meta.columns.duplicated()]
# st.write(df_cpd_meta)
# st.write(df_cpd_infos)
# df_cpds= df_cpd_meta[df_cpd_meta['batchid'].astype(str).str.contains('|'.join(b_list2))].drop_duplicates(subset='pubchemid').reset_index(drop=True)
# st.write(df_cpds)
sql_genes = f"select * from cpdbatchs where batchid  in  (" + ",".join(bq) + ")"
df_cpd_infos=sql_df(sql_genes,conn_meta)
st.session_state['df_cpds'] = df_cpd_infos


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