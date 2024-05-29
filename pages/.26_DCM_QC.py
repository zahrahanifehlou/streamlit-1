import streamlit as st

import sys  # To properly import the files !


import pandas as pd

import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np
import ipywidgets as widgets
from pqdm.processes import pqdm
import glob
import tarfile
import tools
from os import listdir
from pathlib import Path
import polars as pl
import os
import json
import multiprocessing
import io


nb_cpu = multiprocessing.cpu_count()


@st.cache_data
def loadDeepfth(files):
    l2 = files.replace("\\", "/").replace(".fth", "").split("/")[-1].split("_")
    # df = pd.read_feather(files)
    df = pd.read_feather(files)
    # df =cudf.read_feather
    # df = df.drop("tags", axis=1)
    # df["Well"] = l2[1]
    df["Plate"] = l2[0]
    # list_df.append(df)

    # df2 = df.groupby(["Plate"]).median(numeric_only=True)
    # df2 = df2.reset_index()

    return df


@st.cache_data
def loadDeepTar(files):
    with tarfile.open(files, "r") as tar:
        # Replace 'your_file.feather' with the actual file name
        list_df = []
        for member in tar.getmembers():
            # Check if the member is a file and has the '.feather' extension
            if member.isfile() and member.name.endswith(".fth"):
                # Extract the Feather file content
                l = (
                    member.name.replace("\\", "/")
                    .replace(".fth", "")
                    .split("/")[-1]
                    .split("_")
                )
                # print(l)
                feather_content = tar.extractfile(member).read()
                try:
                    df = pd.read_feather(feather_content)
                    # df =cudf.read_feather
                    df = df.drop("tags", axis=1)
                    df["Well"] = l[1]
                    df["Plate"] = l[0]
                    list_df.append(df)
                except:
                    print("error")

        df2 = pd.concat(list_df)
        df2 = df2.groupby(["Plate", "Well"]).median(numeric_only=True)
        df2 = df2.reset_index()

    return df2


sys.path.append("/mnt/shares/L/Code/KsilinkNotebooks/LIB/")
import tools

list_proj = os.listdir("/mnt/shares/L/PROJECTS/")
proj = st.selectbox("Choose your project", list_proj)


paths = sorted(
    Path(f"/mnt/shares/L/PROJECTS/{proj}/Checkout_Results/").iterdir(),
    key=os.path.getmtime,
    reverse=True,
)
# st.write(paths)
paths_clean = [f for f in paths if "test" not in str((f)).lower()]
uploaded_file = st.selectbox(
    "Choose result directory corresponding to this assay", paths_clean
)
list_df = []

files = glob.glob(f"{uploaded_file}/**/*roi.fth", recursive=True)
st.write(files)
for f in files:
    list_df.append(loadDeepfth(f))

# result_deep = pqdm(files, loadDeepfth, n_jobs=min(nb_cpu, 60))

alldata1 = pd.concat(list_df).reset_index(drop=True)
st.write(alldata1)
list_df = []
files2 = glob.glob(f"{uploaded_file}/**/*ag*.fth", recursive=True)
# st.write(files2)
for f in files2:
    list_df.append(loadDeepfth(f))

# result_deep2 = pqdm(files2, loadDeepfth, n_jobs=min(nb_cpu, 60))

alldata = pd.concat(list_df).reset_index(drop=True)

# alldata = pd.concat([alldata, alldata2], axis=1)

cols = [x for x in alldata.columns if "Feature_" in x]
cols_alpha = [x for x in alldata.columns if "Feature_" not in x]
st.write(alldata)
list_cols = alldata.columns
sel_col = st.selectbox("select column to display", list_cols)
