import os
from glob import iglob

import pandas as pd
import plotly.express as px
import pymongo
import streamlit as st
from pymongo import MongoClient

CONNECTION_STRING = "mongodb://192.168.2.127:27017/"

client = MongoClient(CONNECTION_STRING)
db=client.tags.tags

def getPlatesFromCompounds(query):
    matcher=[
        {
            '$project': {
                'data': {
                    '$objectToArray': '$map'
                },
                'plate': '$plate',
                'plateAcq': '$plateAcq',
                'proj': "$meta.project"
            }
        }, {
            '$set': {
                'data': {
                    '$filter': {
                        'input': '$data',
                        'as': 'p',
                        'cond': {
                            '$regexMatch': {
                                'input': '$$p.k',
                                'regex': query,
                                'options': 'i'
                            }
                        }
                    }
                }
            }
        }, {
            '$unwind': {
                'path': '$data',
                'preserveNullAndEmptyArrays': False
            }
        }
    ]
    col=db.aggregate(matcher)
    return pd.DataFrame(list(col)).drop("_id", axis=1)

def getPluginsDataset(df):
    commits={}
    for (project, plate), data in df.groupby(["proj", 'plate']):
        if project not in commits:
            commits[project]=set()
        basepath=f"/mnt/shares/L/PROJECTS/{project}/Checkout_Results/*/ag{plate}.fth"
        for commit in iglob(basepath):
            commit = commit.replace('\\', '/').split('/')[4]
            if commit != "BirdView":
                commits[project].add(commit)
    d=[]
    for (project, plate), data in df.groupby(["proj", 'plate']):
        for commit in commits[project]:
            basepath=f"/mnt/shares/L/PROJECTS/{project}/Checkout_Results/{commit}/*.json"
            for x in iglob(basepath):
                jso=x.replace('\\', '/').split('/')[-1].replace('.json', '').split('_')
                if jso[-1].isnumeric() and jso[-2].isnumeric():
                    d.append((project, "_".join(jso[:-2]), commit, plate, jso[-2], jso[-1]))
                    break
    return pd.DataFrame(d, columns=['Project', 'Plugin', 'Commit', 'Plate', 'Date', 'time'])



var_text = st.text_input("Enter your search")

df=getPlatesFromCompounds(var_text+'.*')
st.write(df)

experiments=getPluginsDataset(df)
st.write(experiments)