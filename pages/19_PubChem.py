import numpy as np
import pandas as pd
import streamlit as st

import pubchempy as pcp

# from streamlib import sql_df
from query_pubchem import QueryPubChem

pubchem = QueryPubChem()

st.title("Pubchem Information")
st.subheader("fetching compounds information from PubChem", divider="rainbow")

compounds = pd.DataFrame()
if "PubChem" not in st.session_state:
    exit(0)
list_cpd = st.session_state["PubChem"]

list_pub = []
for item in list_cpd:
    list_pub.append(pcp.get_compounds(item, "cid", as_dataframe=True))

compounds = pd.concat(list_pub).reset_index()
compounds = compounds.astype({"cid": str})
# st.write(compounds)
st.warning("Getting MOA")
compounds["pubchem_moa"] = pubchem.moa_annotation(compounds["cid"])

st.warning("Getting bioassay")
bioassay_df = pubchem.bioassay(compounds.cid)
with_target = bioassay_df.dropna(subset=["Target GeneID"]).reset_index(drop=True)
st.warning("Getting target, can be long...")
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

st.write(f"Info from Pubchem for your {len(list_cpd)} compounds", target_merge)
