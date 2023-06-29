
import warnings

import numpy as np
import pandas as pd
import psycopg2
import streamlit as st

warnings.filterwarnings("ignore")


def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])
conn = init_connection()

conn_profileDB = psycopg2.connect(
        host="192.168.2.131",
        port="5432",
        user="arno",
        database="ksilink_cpds",
        password="12345",
    )

data = pd.read_sql("select * from cpd", conn)
data.dropna(subset='keggid',inplace=True)
st.write(data.describe())