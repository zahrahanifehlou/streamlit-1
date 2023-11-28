import streamlit as st

st.set_page_config(
    page_title="Welcome to Ksilink Databases",
    page_icon="🧊",
    layout="wide",
    # initial_sidebar_state="expanded",
    # menu_items={
    #     'Get Help': 'https://www.extremelycoolapp.com/help',
    #     'Report a bug': "https://www.extremelycoolapp.com/bug",
    #     'About': "# This is a header. This is an *extremely* cool app!"
    # }
)
st.image("Postgresql.png")
st.write("# Welcome to Ksilink Databases")
st.write('Can be access locally on : http://192.168.2.131:8501; http://192.168.2.130:8501')
# st.sidebar.success("Select a demo above.")
 #- [CRISPR & CPDS](http://192.168.2.131:8501/CRISPR_&_CPDS) compute the similarities between crispr and cpds.
#  - [Pathway Interaction](http://192.168.2.131:8501/Pathway_Interaction) display the involved genes in a specific pathway    
#      - [Clustering](http://192.168.2.131:8501/Clustering) Clustering Jump Profiles with known targets
#     - [Alignement](http://192.168.2.131:8501/Align_Data) Tentative of data alignement between sources (In progress)
st.markdown(
    """
    The goal of this software is to be able to interrogate Jump and historical Dataset.

    **👈 Select an applicaton from the left to see some examples
    of what it can do!

    ### Want to learn more?

    - [Load Data](http://192.168.2.131:8501/Load_Data) : just load fth/csv files, filter them and makes some basic displays
    - [Get Profiles](http://192.168.2.131:8501/Get_Profiles) connect to DB and Get Profiles
    - [Find Profiles](http://192.168.2.131:8501/Find_Similar_Profiles) find profiles close to the ones selected in connect DB from a selected source
    - [Find structures](http://192.168.2.131:8501/Find_Similar_Structures) find the structures and names of the closest profiles
    - [Get Bio Infos](http://192.168.2.131:8501/Get_Biological_Infos) compute the overrepresentation in pathways and connect to >50 Libraries    
    - [Historical Data](http://192.168.2.131:8501/Historical_Data) Retrieve Ksilink Data
    - [Protein-Protein Interactions](http://192.168.2.131:8501/StringDB) Found Protein-Protein Interactions and profile vizualisation
    - [Disease](http://192.168.2.131:8501/Disease) Retrieve Drugs/Genes disease specifics
    - [Convert Data](http://192.168.2.131:8501/Convert_Data) Retrieve Identifiers from different databases (pubchem,Kegg)
    - [Monitoring](http://192.168.2.131:8501/Monitoring) Storage Monitoring
    ### Normal Process:
    - Load Data, Historical data and CRISPR & CPDS are independent
    - Otherwise you should start from ConnectDB...
    - You can also make it independetly and a csv or fth files will be necessary to be uploaded (with relevant column names)

"""
)
tab1,tab2,tab3,tab4=st.tabs(["MetaDatabase", "Sources", "Cpds Per Sources","Profiles Database"])
tab1.image("DB.png")
tab2.image("sources.png")
tab3.image('cpds_per_sources.png')
tab4.image('profile_db.png')
#st.write(st.session_state)

if st.button('Clear Cache'):
    for key in st.session_state.keys():
        del st.session_state[key]