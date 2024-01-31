import pubchempy as pcp
import time
from tqdm import tqdm
import requests
import pandas as pd
import io
from stqdm import stqdm


class QueryPubChem(object):
    def __init__(self):
        pass

    @staticmethod
    def cas_to_cid(field):
        cid_array = []

        for x in stqdm(range(len(field))):
            cas_no = field[x]
            cid_array.append([])
            if not str(cas_no) == "nan":
                temp_cid = pcp.get_cids(cas_no)

                if not temp_cid:
                    temp_cid = pcp.get_cids("CAS-" + cas_no)

                cid_array[-1] += temp_cid
                time.sleep(0.2)  # PubChem database allows only 5 queries per second

            cid_array[-1] = ", ".join([str(x) for x in cid_array[-1]])

        return cid_array

    @staticmethod
    def moa_annotation(field):
        moa = ["" for _ in range(len(field))]

        for x in stqdm(range(len(field))):
            compound = field[x].split(",")

            for y in range(len(list(compound))):
                query_compound = list(compound)[y]

                response = requests.get(
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON?heading=FDA+Pharmacological+Classification".format(
                        query_compound
                    )
                )
                time.sleep(0.2)  # PubChem database allows only 5 queries per second

                try:
                    moa_section = response.json()["Record"]["Section"][0]["Section"][0][
                        "Information"
                    ]
                    for i in range(len(moa_section)):
                        temp = moa_section[i]["Value"]["StringWithMarkup"][0][
                            "String"
                        ].split(" - ")
                        if temp[0] == "Mechanisms of Action [MoA]":
                            if moa[x] == "":
                                moa[x] = temp[-1]
                            else:
                                moa[x] = moa[x] + ", " + temp[-1]
                except KeyError:
                    continue

        return moa

    @staticmethod
    def bioassay(field):
        df = pd.DataFrame()

        for x in stqdm(range(len(field))):
            compound = field[x].split(",")

            for y in range(len(list(compound))):
                query_compound = list(compound)[y]

                response = requests.get(
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/assaysummary/CSV".format(
                        query_compound
                    )
                )
                time.sleep(0.2)  # PubChem database allows only 5 queries per second

                assay = pd.read_csv(
                    io.StringIO(response.content.decode("iso-8859-1")),
                    dtype={"Target GeneID": str},
                )

                try:
                    assay = (
                        assay[assay["Activity Name"].isin(["IC50", "EC50", "AC50"])]
                        .dropna(subset=["Activity Value [uM]"])
                        .reset_index(drop=True)
                    )
                    assay = assay[
                        [
                            "Target GeneID",
                            "AID",
                            "SID",
                            "CID",
                            "Activity Name",
                            "Activity Value [uM]",
                        ]
                    ]
                    assay = assay[assay["Activity Value [uM]"] > 0.001].reset_index(
                        drop=True
                    )
                    assay["InChIKey"] = pcp.get_properties("InChIKey", query_compound)[
                        0
                    ]["InChIKey"]
                    df = df.append(assay).drop_duplicates().reset_index(drop=True)
                except KeyError:
                    continue

        return df

    @staticmethod
    def target(field):
        gene_name = []

        for x in stqdm(range(len(field))):
            aid = field[x]

            response = requests.get(
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{}/targets/GeneSymbol/JSON".format(
                    aid
                )
            )
            time.sleep(0.2)
            gene_name.append(
                response.json()["InformationList"]["Information"][0]["GeneSymbol"][0]
            )

        return gene_name
