import streamlit as st

st.write("# Welcome to Ksilink Databases")
st.write('Can be access on : http://192.168.2.131:8501; http://192.168.2.130:8501; http://192.168.2.128:8501')
# st.sidebar.success("Select a demo above.")

st.markdown(
    """
    The goal of this software is to be able to interrogate Jump and historical Dataset.

    **👈 Select an applicaton from the left to see some examples
    of what it can do!

    ### Want to learn more?

    - [Load Data](http://192.168.2.131:8501/Load_Data) : just load fth/csv files, filter them and makes some basic displays
    - [Connect DB](http://192.168.2.131:8501/Connect_DB) connect to DB and get data from DBs
    - [Plot Profiles](http://192.168.2.131:8501/Plot_Profiles) plot profiles close to the ones selected in connect DB from a selected source
    - [Get structures](http://192.168.2.131:8501/Get_Structures) display the structures and names of the closest profiles
    - [Get Bio Infos](http://192.168.2.131:8501/Get_Biological_Infos) compute the overrepresation in pathways and connect to >50 Libraries
    - [Pathway Interaction](http://192.168.2.131:8501/Pathway_Interaction) display the involved genes in a specific pathway
    - [CRISPR & CPDS](http://192.168.2.131:8501/CRISPR_&_CPDS) compute the similarities between crispr and cpds.
    - [Historical Data](http://192.168.2.131:8501/Historical_Data) Gather all available Ksilink data

    ### Normal Process:
    - Load Data, Historical data and CRISPR & CPDS are independent
    - Otherwise you should start from ConnectDB...
    - You can also make it independetly and a csv or fth files will be necessary to be uploaded (with relevant column names)

"""
)
tab1,tab2,tab3=st.tabs(["Database", "Sources", "Cpds Per Sources"])
tab1.image("DB.png")
tab2.image("sources.png")
tab3.image('cpds_per_sources.png')
#st.write(st.session_state)

if st.button('Clear Cache'):
    for key in st.session_state.keys():
        del st.session_state[key]