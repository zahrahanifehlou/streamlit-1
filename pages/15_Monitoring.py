import numpy as np
import pandas as pd
import plotly.express as px
import pymongo
import streamlit as st
from pymongo import MongoClient

st.set_page_config(
    layout="wide",
)

st.header("Storage Monitoring",divider='rainbow')


CONNECTION_STRING = "mongodb://192.168.2.127:27017/"

client = MongoClient(CONNECTION_STRING)
db=client.tags.tags
pipeline=[
    {
        '$project': {
            'size': '$size',
            'project': "$meta.project",
            'plate':'$plate',
            'serverPath': {
                '$substr': [
                    '$serverPath', 0, 1
                ]
            },
            'acquisitionDate': '$acquisitionDate'
        }
    }
]
col=db.aggregate(pipeline)
dfSizes=pd.DataFrame(list(col)).sort_values("acquisitionDate").reset_index(drop=True)
fig = px.histogram(dfSizes, x="acquisitionDate", y="size",
                    title="Storage activities",color='project',template='plotly_dark')
fig.update_traces(xbins_size="M2")
fig.update_xaxes(showgrid=True, ticklabelmode="period",
                 dtick="M2", tickformat="%b\n%Y")
fig.update_layout(bargap=0.1)
st.plotly_chart(fig, theme="streamlit", use_container_width=True)

dfproj=dfSizes[["project", "size","serverPath"]].groupby(["project","serverPath"]).sum().reset_index()
dfproj =dfproj[(dfproj != 0).all(1)]

fig23 = px.treemap(dfproj,path=['serverPath','project'],values='size',color='project',
                   title='Projects TreeMap by storage'
                  )
fig23.update_traces(marker=dict(cornerradius=5))
st.plotly_chart(fig23, theme="streamlit", use_container_width=True)