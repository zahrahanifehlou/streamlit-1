import pandas as pd
import streamlit as st

from Bio.KEGG.REST import kegg_get
from streamlib import sql_df


conn_meta = "postgres://arno:12345@192.168.2.131:5432/ksi_cpds"
conn_profileDB = "postgres://arno:12345@192.168.2.131:5432/ksilink_cpds"


# A bit of helper code to shorten long text
def head(text):
    """Print the first lines lines of the passed text."""
    stre = text.split("\n")
    return stre


def pub(sel_dis):
    keg_dis = head(kegg_get(sel_dis).read())
    list_ref = []
    list_auth = []
    list_title = []
    list_jour = []
    for item in keg_dis:
        if "REFERENCE" in item:
            list_ref.append(
                "https://pubmed.ncbi.nlm.nih.gov/"
                + item.split("REFERENCE")[1].split(":")[1]
            )
        if "AUTHORS" in item:
            list_auth.append(item.split("AUTHORS")[1])
        if "TITLE" in item:
            list_title.append(item.split("TITLE")[1])
        if "JOURNAL" in item:
            list_jour.append(item.split("JOURNAL")[1])
    df_pub = pd.DataFrame()
    df_pub["REFERENCE"] = list_ref
    df_pub["AUTHORS"] = list_auth
    df_pub["TITLE"] = list_title
    df_pub["JOURNAL"] = list_jour
    return df_pub


st.title("Diseases Information")
st.subheader("From Disease to CPDS/Genes", divider="rainbow")
sql_dis = "select * from disease"

df_dis = sql_df(sql_dis, conn_meta)

# st.write(df_dis)
sel_disease = st.selectbox("Chose the disease:", df_dis["disname"].unique())

code_dis = df_dis[df_dis["disname"] == sel_disease]["disid"].values[0]
st.write("---")
st.write("Kegg Code for selected disease: ", code_dis)
tog_ref = st.sidebar.toggle("Get References", help="retrieve main litterature")
if tog_ref:
    df_pubmed = pub(code_dis)
    st.write("REFERENCES:", df_pubmed)
st.write("---")

col1, col2 = st.columns(2)

sql_known_drugs = f"select * from keggcpddis where disid='{code_dis}'"
df_known_drug = sql_df(sql_known_drugs, conn_meta)
# st.write(df_known_drug)


with col1:
    st.write("---")
    st.write("Known Drugs for this disease", df_known_drug)
    st.write("---")


#################################GET Genes ##################

sql_known_genes = f"select * from genedis where disid='{code_dis}'"
df_link_genes = sql_df(sql_known_genes, conn_meta)


with col2:
    st.write("---")
    st.write("Known genes for this disease", df_link_genes)
    st.write("---")

#

on = st.sidebar.toggle("Retrieve Cpds Info", help="Will retrieve all infos")
if on and not df_known_drug.empty:
    # @st.cache_data
    def convert_df(df):
        # IMPORTANT: Cache the conversion to prevent computation on every rerun
        return df.to_csv().encode("utf-8")

    # @st.cache_data
    def get_cpds(list_drug):
        sql_ksi = (
            f"select * from platemap where platemap.batchid in ({','.join(list_drug)})"
        )
        df_ksi = sql_df(sql_ksi, conn_meta)

        return df_ksi

    list_drug = [f"'{t}'" for t in df_known_drug["keggid"]]
    # sql_drug="select * from keggcpd where keggid  in (" + ",".join(list_drug) + ")"
    sql_drug = f"SELECT keggcpd.*, keggcpdgene.geneid,gene.geneid, gene.symbol,cpdbatchs.batchid FROM keggcpd\
            LEFT JOIN keggcpdgene ON keggcpdgene.keggid=keggcpd.keggid\
            LEFT JOIN gene ON gene.geneid=keggcpdgene.geneid\
            LEFT JOIN cpdbatchs ON cpdbatchs.pubchemid=keggcpd.pubchemid\
            WHERE keggcpd.keggid IN ({','.join(list_drug)})"
    # sql_drug="select * from keggcpd"sql_df
    df_drug = sql_df(sql_drug, conn_meta)
    df_drug = df_drug.loc[:, ~df_drug.columns.duplicated()].copy()
    df_drug = df_drug.drop_duplicates(subset=["keggid", "batchid"]).reset_index(
        drop=True
    )
    st.subheader("Drug Infos from Kegg", divider="rainbow")
    st.write(df_drug)
    # on_save1 = st.sidebar.toggle('Save Drug Infos')

    list_drug_ksi = [f"'{t}'" for t in df_drug["keggid"]]
    sql_prof = f"select * from cpd where keggid IN ({','.join(list_drug_ksi)})"
    df_ksi = sql_df(sql_prof, conn_meta)
    # st.write("df_ksi", df_ksi)

    list_pub = [f"'{t}'" for t in df_ksi["pubchemid"]]
    sql_pub = f"select * from cpdbatchs where pubchemid IN ({','.join(list_pub)})"
    df_ksi2 = sql_df(sql_pub, conn_meta)

    df_ksi2 = df_ksi2.drop_duplicates()
    # st.write("df_ksi", df_ksi2)
    list_pub2 = [f"'{t}'" for t in df_ksi2["batchid"]]
    sql_pub3 = f"select * from platemap where batchid IN ({','.join(list_pub2)})"
    df_ksi_f = sql_df(sql_pub3, conn_meta)
    # st.write("df_ksi", df_ksi_f)

    list_source = df_ksi_f.assay.unique()
    st.subheader("Data from Jump", divider="rainbow")
    sel_source = st.selectbox("Chose the source", list_source)
    df_sel = df_ksi_f[df_ksi_f["assay"] == sel_source]
    st.write(f"Data from source: {sel_source}", df_sel)

    # TO DO : Get profiles and get from project profiles
    sql_pro_prof = (
        f"select * from projectsprofile where batchid IN ({','.join(list_pub2)})"
    )
    df_ksi_proj = sql_df(sql_pro_prof, conn_profileDB)
    list_proj = df_ksi_proj.project.unique()
    st.subheader("In House Projects", divider="rainbow")
    sel_project = st.selectbox("Chose the project", list_proj)
    df_ksi_proj = df_ksi_proj[df_ksi_proj["project"] == sel_project]

    st.write(df_ksi_proj)
    # on_save = st.sidebar.toggle('Save Data from source')


on_gene = st.sidebar.toggle("Retrieve Genes Info", help="Will retrieve all gene infos")
if on_gene and not df_link_genes.empty:
    list_genes = [f"'{t.upper()}'" for t in df_link_genes["geneid"]]
    sql_genes = "select * from gene where geneid  in (" + ",".join(list_genes) + ")"
    # sql_drug="select * from keggcpd"
    df_genes = sql_df(sql_genes, conn_meta)
    st.write("Genes Infos", df_genes)

st.subheader("From gene to diseases", divider="rainbow")
sql_help = "select gene.symbol from gene"
df_gene_help = sql_df(sql_help, conn_meta)
list_help_genes = df_gene_help["symbol"]
gene_of_int = st.selectbox("chose one of them", list_help_genes)
# gene_of_int= st.text_input('enter your gene of interest').upper()
sql_known_dis = f"select * from gene where symbol='{gene_of_int}'"
df_gene_of_int = sql_df(sql_known_dis, conn_meta)
if not df_gene_of_int.empty:
    geneid = df_gene_of_int[df_gene_of_int["symbol"] == gene_of_int]["geneid"].values[0]
    # st.write(symbol)
    sql_known_gene_in_dis = f"select * from genedis where geneid='{geneid}'"
    df_gene_in_dis = sql_df(sql_known_gene_in_dis, conn_meta)
    if not df_gene_in_dis.empty:
        disid = df_gene_in_dis["disid"].values
        df_list_dis = df_dis[df_dis["disid"].isin(disid)].reset_index(drop=True)
        st.write("known diseases involving this gene", df_list_dis)
        sel_dis = st.selectbox(
            "select the disease for references", df_list_dis["disid"].to_list()
        )
        if tog_ref:
            df_pub = pub(sel_dis)
            st.write("REFERENCES", df_pub)
else:
    st.warning("Unknown gene")


# from tdc.multi_pred import GDA
# data = GDA(name = 'DisGeNET')
# # data2  =data[data['Disease_ID']=='C1449563']
# st.write('TDC',data.get_data())
# C1449563
# st.write(split)
