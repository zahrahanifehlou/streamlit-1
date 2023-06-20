import streamlit as st

st.write("# Welcome to Ksilink Databases")
# st.sidebar.success("Select a demo above.")

st.markdown(
    """
    The goal of this software is to be able to interrogate Jump Dataset.

    **👈 Select an applicaton from the left to see some examples
    of what it can do!

    ### Want to learn more?

    - [Load Data](http://192.168.2.131:8501/Load_Data) : just load fth/csv files, filter them and makes some basic displays
    - [Connect DB](http://192.168.2.131:8501/Connect_DB) will connect to DB and get data from DBs
    - [Plot Profiles](http://192.168.2.131:8501/Plot_Profiles) allow to plot profiles close to the ones selected in connect DB from a selected source
    - [Get structures](http://192.168.2.131:8501/Get_Structures) will display the structures and names of the closest profiles
    - [Get Bio Infos](http://192.168.2.131:8501/Get_Biological_Infos) will display the pathways and genes


"""
)
st.image("DB.png")
st.write(st.session_state)

if st.button('Clear Cache'):
    for key in st.session_state.keys():
        del st.session_state[key]