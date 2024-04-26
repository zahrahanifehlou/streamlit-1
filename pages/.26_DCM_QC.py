import streamlit as st

import sys  # To properly import the files !

# For tools modules
sys.path.append("/mnt/shares/L/Code/KsilinkNotebooks/LIB/")


# For exp_graphs & deep_loader modules
sys.path.append("/mnt/shares/L/Code/KsilinkNotebooks/DCM/helper_lib")

import pandas as pd
from glob import iglob as glob
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np
import ipywidgets as widgets

import tools

from exp_graphs import *
from deep_loader import *

import json


def load_deep_features(buffer, well, plate):  # Load all the features from the well
    tmp_df = pd.read_feather(buffer)

    tmp_df["Well"] = well
    tmp_df["Plate"] = plate

    return tmp_df


list_files = st.file_uploader("Load your files", accept_multiple_files=True)
if list_files:
    df = pd.concat([pd.read_feather(x) for x in list_files], ignore_index=True)
    cols = {x for x in df.columns}
    df = tools.getCategories(df)
    categories = [x for x in df.columns if x not in cols]
    st.write("Categories: ", df)
    df.loc[df.cell.isna(), "cell"] = df.CellLines[df.cell.isna()]
    df.loc[df.exp.isna(), "exp"] = df.Plate[df.exp.isna()].apply(
        lambda x: x.split("-")[1]
    )
    train_paths = list(set(["/".join(x.split("\\")[:-1]) for x in list_files]))
    print(train_paths)

    classes = {
        "All": {
            "Plate": lambda x: True
            if "16" in x or "21" in x or "20" in x or "19" in x
            else False
        },
    }

    train_df = load_features_from(
        train_paths, classes, file_sub="roi", loader=load_deep_features
    )
    train_df.loc[train_df.cell.isna(), "cell"] = train_df.loc[
        train_df.cell.isna(), "CellLines"
    ]
    train_df.loc[train_df.exp.isna(), "exp"] = train_df.Plate[
        train_df.exp.isna()
    ].apply(lambda x: x.split("-")[1])
    train_df["disp"] = train_df[list(["cell", "exp"])].apply(
        lambda row: " ".join(row.values.astype(str)), axis=1
    )
    train_df["disp"] = train_df[list(["cell", "exp"])].apply(
        lambda row: " ".join(row.values.astype(str)), axis=1
    )

    sdfa = train_df.groupby(list(["cell", "exp"]) + ["Well", "Plate"], as_index=False)[
        "Area"
    ].count()
    sdfa.rename(columns={"Area": "Nucleis"}, inplace=True)
    sdfa["disp"] = sdfa[list(["cell", "exp"])].apply(
        lambda row: " ".join(row.values.astype(str)), axis=1
    )
    st.write("sdfa", sdfa)
