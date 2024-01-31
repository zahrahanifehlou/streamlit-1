import sys
import numpy as np
import pandas as pd
import streamlit as st
import time
from tqdm import tqdm
import requests
import pubchempy as pcp

sys.path.append("/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit-1/lib/")
# from streamlib import sql_df
from query_pubchem import QueryPubChem

pubchem = QueryPubChem()


# def moa_annotation(field):
#     moa = ["" for _ in range(len(field))]

#     for x in range(len(field)):
#         compound = field[x].split(",")
#         # st.write(compound)
#         for y in range(len(list(compound))):
#             query_compound = list(compound)[y]
#             # st.write(query_compound)
#             response = requests.get(
#                 "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON?heading=FDA+Pharmacological+Classification".format(
#                     query_compound
#                 )
#             )
#             time.sleep(0.2)  # PubChem database allows only 5 queries per second
#             # st.write(response)
#             try:
#                 moa_section = response.json()["Record"]["Section"][0]["Section"][0][
#                     "Information"
#                 ]
#                 for i in range(len(moa_section)):
#                     temp = moa_section[i]["Value"]["StringWithMarkup"][0][
#                         "String"
#                     ].split(" - ")
#                     if temp[0] == "Mechanisms of Action [MoA]":
#                         if moa[x] == "":
#                             moa[x] = temp[-1]
#                         else:
#                             moa[x] = moa[x] + ", " + temp[-1]
#             except KeyError:
#                 continue

#     return moa


compounds = pd.DataFrame()
if "PubChem" not in st.session_state:
    exit(0)
list_cpd = st.session_state["PubChem"]
# list_cpd = [
#     "2979901",
#     "24772097",
#     "11691242",
#     "5154691",
#     "1179237",
#     "1104205",
#     "16070100",
#     "2915164",
#     "41018128",
#     "4996",
# ]
# compounds["pubchem_cid"] = list_cpd
list_pub = []
for item in list_cpd:
    list_pub.append(pcp.get_compounds(item, "cid", as_dataframe=True))

compounds = pd.concat(list_pub).reset_index()
compounds = compounds.astype({"cid": str})
# st.write(compounds)
compounds["pubchem_moa"] = pubchem.moa_annotation(compounds["cid"])

bioassay_df = pubchem.bioassay(compounds.cid)
with_target = bioassay_df.dropna(subset=["Target GeneID"]).reset_index(drop=True)
with_target["GeneName"] = pubchem.target(with_target.AID)
# st.write(with_target)
median_activity_value = (
    with_target[["InChIKey", "CID", "Activity Name", "Activity Value [uM]", "GeneName"]]
    .groupby(["CID", "InChIKey", "GeneName", "Activity Name"])
    .agg({"Activity Value [uM]": lambda x: np.median(x)})
    .reset_index()
)


# target_merge = (
#     median_activity_value[
#         ["InChIKey", "CID", "Activity Name", "Activity Value [uM]", "GeneName"]
#     ].groupby(["CID", "InChIKey", "Activity Name"])
#     # .agg(
#     #     {
#     #         "GeneName": lambda x: ",".join(list(set([i.upper() for i in x]))),
#     #         "Activity Value [uM]": lambda y: ",".join(
#     #             [str(np.around(j, 3)) for j in y]
#     #         ),
#     #     }
#     # )
#     # .reset_index()
# )
# target_merge = target_merge.reset_index()
median_activity_value.rename(
    columns={
        "CID": "pubchem_cid",
        "Activity Name": "Activity Type",
        "GeneName": "Gene",
        "Activity Value [uM]": "Median Activity Value [uM]",
    },
    inplace=True,
)

target_merge = median_activity_value.astype({"pubchem_cid": str})

target_merge = target_merge.merge(
    compounds[["iupac_name", "cid", "canonical_smiles", "pubchem_moa"]],
    left_on="pubchem_cid",
    right_on="cid",
)

# final_table = target_merge[
#     [
#         "Name",
#         "CAS",
#         "pubchem_cid",
#         "Canonical_Smiles",
#         "InChIKey",
#         "Target",
#         "Gene",
#         "Activity Type",
#         "Median Activity Value [uM]",
#     ]
# ]

st.write(target_merge)
