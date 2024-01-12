

import sys

import pandas as pd
import psycopg2
import streamlit as st
import polars as pl
import streamlit.components.v1 as components
from pygwalker.api.streamlit import init_streamlit_comm, get_streamlit_html
from os import listdir
from pathlib import Path
import os
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import MinMaxScaler
import numpy as np
from pqdm.processes import pqdm
import glob
import tarfile
import io
import plotly.express as px

st.set_page_config(
    layout="wide",
)
sys.path.append('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/')
init_streamlit_comm()

@st.cache_resource
def get_pyg_html(df: pd.DataFrame) -> str:
    # When you need to publish your application, you need set `debug=False`,prevent other users to write your config file.
    # If you want to use feature of saving chart config, set `debug=True`
    html = get_streamlit_html(df,spec="./chart_meta_0.json", themeKey='vega',use_kernel_calc=True, debug=False)
    return html

@st.cache_resource
def sql_df(sql_str, _conn):

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


st.header("In house Dataset",divider='rainbow')

sql_proj='select project from projectsprofile group by project'
res_proj=sql_df(sql_proj, profile_conn)


# st.write(res_proj)
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
pivot_df=pivot_df.apply(pd.to_numeric, errors='ignore')
components.html(get_pyg_html(pivot_df), height=1000, scrolling=True)

st.header("Dataset from Phenolink not validated",divider='rainbow')

on = st.sidebar.toggle('Search Data')
dl = st.sidebar.toggle('Deep Learning Data')
cbc = st.sidebar.toggle('Cell by Cell Data')
tp = st.sidebar.toggle('Temporal Data')
list_proj=os.listdir('/mnt/shares/L/Projects/')
proj= st.selectbox('Choose your project',list_proj)

if on and tp:
   
    paths = sorted(Path(f'/mnt/shares/L/Projects/{proj}/Checkout_Results/').iterdir(), key=os.path.getmtime, reverse=True)
# st.write(paths)
    paths_clean=[f for f in paths if 'test' not in str((f)).lower()]
    uploaded_files = st.multiselect("Choose result directories corresponding to this assay", paths_clean)
    
    list_df=[]
    if uploaded_files:
        for item in uploaded_files:
            fth_file=[f for f in listdir(item) if f.startswith('ag')]
            if len(fth_file)>0:
                for fi in fth_file:
                    # st.write(item.as_posix())
                    if fi.endswith('.fth'):
                        list_df.append(pd.read_feather(item.as_posix()+'/'+fi))
                    else:
                        # st.write(item.as_posix()+'/'+fi)
                        list_df.append(pd.read_csv(item.as_posix()+'/'+fi, encoding='latin1'))
            else:
                st.warning('No data started with ag in this directory',icon="🚨" )
    if len(list_df)>0:
        df_data =  pd.concat(list_df)
        len_1=len(df_data)
        df_data.replace([np.inf, -np.inf], np.nan, inplace=True)

        df_data= df_data.dropna()
        len_2=len(df_data)
        df_data=df_data.apply(pd.to_numeric,errors='ignore')
        if len_1!=len_2:
            st.warning(f'{len_1-len_2} rows were removed due to NaN values')
        # Create correlation matrix
        
    else:
        exit(0)
    st.warning("In dev....")
    data=df_data.copy()
    numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
    
    col_sel = data.select_dtypes(include=numerics).columns.to_list()
    t2 = st.radio("Analyse", ["Contractility","Calcium"],0)
    dt=0.031
    df_agg=data.copy()
    # df_agg=data[col_sel]
    if t2=='Calcium':
        df_agg[col_sel] = (-df_agg[col_sel]).add(df_agg[col_sel].min(axis=1), axis = 0).add(df_agg[col_sel].max(axis=1), axis = 0)
    # df_agg.drop('tags',axis=1,inplace=True)
    time_cols = sorted([col for col in df_agg.columns if 'time' in col],key=lambda x: int(x.split("_")[-1]))
    sig_x=[i*dt for i in range(len(time_cols))]
    
    # df_agg['tags']=data['tags']
    # 
    listofexp=df_agg['Plate'].unique()
    sel_col_exp= st.selectbox("Chose exps",listofexp)
    df_agg2=df_agg[df_agg['Plate']==sel_col_exp]
    df_agg2=df_agg2.groupby(['tags'])
    for name, group in df_agg2:     
        title=group.tags.values[0]         
        cpd_names = group.Well.values
        # st.write(cpd_names)
        df_plt = group.set_index("Well")
        df_plt = df_plt[time_cols].T
        df_plt['t(s)']=sig_x
        st.plotly_chart(px.line(
            df_plt,
            x='t(s)',
            y=cpd_names,
            width=800,
            height=800,
            title=title,
            # line_shape='hv'
            
        ), theme="streamlit", use_container_width=True)
    # exit(0)


if on and not dl and not cbc and not tp:
    paths = sorted(Path(f'/mnt/shares/L/Projects/{proj}/Checkout_Results/').iterdir(), key=os.path.getmtime, reverse=True)
# st.write(paths)
    paths_clean=[f for f in paths if 'test' not in str((f)).lower()]
    uploaded_files = st.multiselect("Choose result directories corresponding to this assay", paths_clean)

    list_df=[]
    if uploaded_files:
        for item in uploaded_files:
            fth_file=[f for f in listdir(item) if f.startswith('ag')]
            if len(fth_file)>0:
                for fi in fth_file:
                    # st.write(item.as_posix())
                    if fi.endswith('.fth'):
                        list_df.append(pd.read_feather(item.as_posix()+'/'+fi))
                    else:
                        # st.write(item.as_posix()+'/'+fi)
                        list_df.append(pd.read_csv(item.as_posix()+'/'+fi, encoding='latin1'))
            else:
                st.warning('No data started with ag in this directory',icon="🚨" )
        if len(list_df)>0:
            df_data =  pd.concat(list_df)
            len_1=len(df_data)
            df_data.replace([np.inf, -np.inf], np.nan, inplace=True)

            df_data= df_data.dropna()
            len_2=len(df_data)
            df_data=df_data.apply(pd.to_numeric,errors='ignore')
            if len_1!=len_2:
                st.warning(f'{len_1-len_2} rows were removed due to NaN values',icon="🚨")
            len_col=len(df_data.columns)
            corr_matrix = df_data.corr().abs()

            # Select upper triangle of correlation matrix
            upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(np.bool))

            # Find index of feature columns with correlation greater than 0.95
            to_drop = [column for column in upper.columns if any(upper[column] > 0.95)]

            # Drop features 
            df_data = df_data.drop(df_data[to_drop], axis=1)
            len_col2=len(df_data.columns)
            if len_col!=len_col2:
                st.warning(f'{len_col-len_col2} out of {len_col} columns were removed due to correlation > |0.95| ',icon="🚨")
            components.html(get_pyg_html(df_data), height=1000, scrolling=True)
            umap_on=st.sidebar.toggle('UMAP')
            if umap_on:
                
                df_data.replace([np.inf, -np.inf], np.nan, inplace=True)
                df_data.dropna(inplace=True, axis=1)
                df_data.dropna(inplace=True)
                # st.write('Data after removing Nan',df_data.head(2))
                if len(df_data) > 2:
                    selector = VarianceThreshold()
                    numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
                    df_num = df_data.select_dtypes(include=numerics)
                    num_col = [x for x in df_num.columns if "Meta" not in x and "meta" not in x]
                    df_num = df_num[num_col]
                    dat_sel = selector.fit_transform(df_num)
                    col_sel = selector.get_feature_names_out()
                    data1 = pd.DataFrame(dat_sel, columns=col_sel).reset_index(drop=True)
                    # st.write('Data after removing Meta, Var=0',data1.head(2))
                    cols_alpha = df_data.select_dtypes(exclude=numerics).columns
                    data = pd.concat([data1, df_data[cols_alpha].reset_index(drop=True)], axis=1)
                    scaler = MinMaxScaler()
                    col_sel = data.select_dtypes(include=numerics).columns.to_list()
                    data_scaled = pd.DataFrame(
                        scaler.fit_transform(data.select_dtypes(include=numerics)),
                        columns=col_sel,
                    )
                    data_scaled = pd.concat([data_scaled, data[cols_alpha].reset_index(drop=True)], axis=1)
                    import umap

            #
                    model = umap.UMAP(random_state=42, verbose=False).fit(data_scaled[col_sel])
                    emb = model.transform(data_scaled[col_sel])
                    # st.write('UMAP VISU')
                    df_all_umap = pd.DataFrame()
                    df_all_umap["X_umap"] = emb[:, 0]
                    df_all_umap["Y_umap"] = emb[:, 1]
                    df_all_umap = pd.concat([df_all_umap, data[cols_alpha].reset_index(drop=True)], axis=1)
                    components.html(get_pyg_html(df_all_umap), height=1000, scrolling=True)

@st.cache_data
def loadDeepTar(files):
    with tarfile.open(files, 'r') as tar:
        list_df=[]
    # Replace 'your_file.feather' with the actual file name
        for member in tar.getmembers():
            # Check if the member is a file and has the '.feather' extension
            if member.isfile() and member.name.endswith('.fth'):
                # Extract the Feather file content
                l = member.name.replace("\\", "/").replace(".fth", "").split("/")[-1].split("_")
                # print(l)
                feather_content = tar.extractfile(member).read()
                try:
                    df = pd.read_feather(io.BytesIO(feather_content))
                    # df =cudf.read_feather
                    df = df.drop('tags',axis=1)
                    df["Well"] = l[1]
                    df["Plate"] = l[0]
                    list_df.append(df)
                except:
                    print("error")

    df2=pd.concat(list_df).groupby(["Plate", "Well"]).median()
    df2=df2.reset_index()
    
    return df2


if on and dl:
    sys.path.append("/mnt/shares/L/Code/KsilinkNotebooks/LIB/")
    import tools
    
    paths = sorted(Path(f'/mnt/shares/L/Projects/{proj}/Checkout_Results/').iterdir(), key=os.path.getmtime, reverse=True)
    # st.write(paths)
    paths_clean=[f for f in paths if 'test' not in str((f)).lower()]
    uploaded_file = st.selectbox("Choose result directory corresponding to this assay", paths_clean)
    list_df=[]
    if uploaded_file:
        files = glob.glob(f"{uploaded_file}/**/*deep_features.tar", 
                   recursive = True)
        
        if len(files)>0:
            with st.spinner(f'Wait for it... Loading {len(files)} files and computing UMAP'):
                result_deep = pqdm(files, loadDeepTar, n_jobs=20)
                alldata = pd.concat(result_deep).reset_index(drop=True)
                alldata = tools.setCategories(tools.retrieve_tags(alldata))
                alldata=tools.getScreenCategories(alldata)
                cols = [x for x in alldata.columns if 'Feature_' in x]
                cols_alpha = [x for x in alldata.columns if 'Feature_' not in x]
                st.write(alldata.sample(5))
                import umap

                emb = umap.UMAP(random_state=42, verbose=False).fit_transform(alldata[cols])
                df_all_umap = pd.DataFrame()
                df_all_umap["X"] = emb[:, 0]
                df_all_umap["Y"] = emb[:, 1]
                df_all_umap[cols_alpha]=alldata[cols_alpha]

                components.html(get_pyg_html(df_all_umap), height=1000, scrolling=True)
            st.success('Done!')


@st.cache_data
def loadCellbyCell(file):
    df2 = pd.read_feather(file)
    l = file.replace("\\", "/").replace(".fth", "").split("/")[-1].split("_")
    df2["Well"] = l[1]
    df2["Plate"] = l[0]
    return df2



if on and cbc:
    sys.path.append("/mnt/shares/L/Code/KsilinkNotebooks/LIB/")
    import tools    
    paths = sorted(Path(f'/mnt/shares/L/Projects/{proj}/Checkout_Results/').iterdir(), key=os.path.getmtime, reverse=True)
    # st.write(paths)
    paths_clean=[f for f in paths if 'test' not in str((f)).lower()]
    uploaded_file = st.selectbox("Choose result directory corresponding to this assay", paths_clean)
    list_df=[]
    if uploaded_file:
        files = glob.glob(f"{uploaded_file}/**/*cell_by_cell.fth", 
                   recursive = True)
        
        if len(files)>0:
            with st.spinner(f'Wait for it... Loading {len(files)} files'):
                result_deep=pqdm(files, loadCellbyCell,n_jobs=20)
                df_data = pd.concat(result_deep).reset_index(drop=True)
                df_data = tools.setCategories(tools.retrieve_tags(df_data))
                df_data=tools.getScreenCategories(df_data)
                # st.write('df_data',df_data)
                len_1=len(df_data)
                df_data.replace([np.inf, -np.inf], np.nan, inplace=True)

                df_data= df_data.dropna()
                len_2=len(df_data)
                df_data=df_data.apply(pd.to_numeric,errors='ignore')
                if len_1!=len_2:
                    st.warning(f'{len_1-len_2} rows were removed due to NaN values',icon="🚨")
                len_col=len(df_data.columns)
                corr_matrix = df_data.corr().abs()

                # Select upper triangle of correlation matrix
                upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(np.bool))

                # Find index of feature columns with correlation greater than 0.95
                to_drop = [column for column in upper.columns if any(upper[column] > 0.95)]

                # Drop features 
                df_data = df_data.drop(df_data[to_drop], axis=1)
                len_col2=len(df_data.columns)
                if len_col!=len_col2:
                    st.warning(f'{len_col-len_col2} out of {len_col} columns were removed due to correlation > |0.95| ',icon="🚨")
                    components.html(get_pyg_html(df_data), height=1000, scrolling=True)
            st.success('Done!')
