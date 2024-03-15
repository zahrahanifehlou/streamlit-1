import pandas as pd
import json
import requests
import streamlib as st
from sklearn.metrics.pairwise import cosine_similarity
import plotly.express as px
import umap
import polars as pl


################################### CONNECTION DBs ##############################################


def int_to_str(ints):
    return str(ints)


def str_to_float(strs):
    return float(strs)


def sql_df(sql_str, conn):
    df_d = pl.read_database_uri(query=sql_str, uri=conn)
    return df_d.to_pandas()


################################### STRINGDB PPI ##############################################


def get_stringDB(df_all_umap, thresh=0.7, genecol="target"):
    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    params = {
        "identifiers": "\r".join(df_all_umap[genecol].to_list()),
        # "identifiers" : "\r".join(["p53", "BRCA1", "cdk2", "Q99835"]), # your protein list
        "species": 9606,  # species NCBI identifier
        "limit": 1,  # only one (best) identifier per input protein
        "echo_query": 1,  # see your input identifiers in the output
        "caller_identity": "www.awesome_app.org",  # your app name
    }

    request_url = "/".join([string_api_url, output_format, method])
    st.write(request_url)
    results = requests.post(request_url, data=params)

    list_id0 = []
    list_id1 = []
    list_inter = []
    list_edges = []
    for line in results.text.strip().split("\n"):
        l2 = line.strip().split("\t")

        if len(l2) > 4:
            p1, p2 = l2[2], l2[3]

            # filter the interaction according to experimental score
            experimental_score = float(l2[10])
            if experimental_score >= thresh:
                # print
                # print("\t".join([p1, p2, "experimentally confirmed (prob. %.3f)" % experimental_score]))
                list_id0.append(p1)
                list_id1.append(p2)
                list_edges.append((p1, p2))
                list_inter.append(experimental_score)
    return list(set(list_edges))


################################### STRING DB ENRICHMENT ##############################################


def get_stringDB_enr(df_all_umap, genecol="target", cat="KEGG", fdra=0.01):
    if not isinstance(df_all_umap, pd.DataFrame):
        list_genes = df_all_umap
        df_all_umap = pd.DataFrame()
        df_all_umap[genecol] = list_genes

    string_api_url = "https://string-db.org/api"
    output_format = "json"
    method = "enrichment"
    params = {
        "identifiers": "\r".join(df_all_umap[genecol].to_list()),
        # "identifiers" : "\r".join(["p53", "BRCA1", "cdk2", "Q99835"]), # your protein list
        "species": 9606,  # species NCBI identifier
        "limit": 1,  # only one (best) identifier per input protein
        "echo_query": 1,  # see your input identifiers in the output
        "caller_identity": "www.awesome_app.org",  # your app name
    }
    request_url = "/".join([string_api_url, output_format, method])

    response = requests.post(request_url, data=params)
    data = json.loads(response.text)
    list_go = []
    list_names = []
    list_dec = []
    list_fdr = []
    list_category = []
    list_p_value = []
    for row in data:
        term = row["term"]
        preferred_names = ",".join(row["preferredNames"])
        fdr = float(row["fdr"])
        description = row["description"]
        category = row["category"]
        p_value = row["p_value"]

        if category == cat and fdr < fdra:
            list_go.append(term)
            list_names.append(preferred_names)
            list_fdr.append(format(fdr, "5.2e"))
            list_p_value.append(format(p_value, "5.2e"))
            list_dec.append(description)
            list_category.append(category)
            # st.write("\t".join([term,preferred_names, str(fdr), description]))
    df_go = pd.DataFrame()
    df_go["GO Term"] = list_go
    df_go["Preferred Names"] = list_names
    df_go["Description"] = list_dec
    df_go["fdr"] = list_fdr
    df_go["p_val"] = list_p_value
    df_go["category"] = list_category
    return df_go


def get_list_category(df_all_umap, genecol="target"):
    string_api_url = "https://string-db.org/api"
    output_format = "json"
    method = "enrichment"
    params = {
        "identifiers": "\r".join(df_all_umap[genecol].to_list()),
        # "identifiers" : "\r".join(["p53", "BRCA1", "cdk2", "Q99835"]), # your protein list
        "species": 9606,  # species NCBI identifier
        "limit": 1,  # only one (best) identifier per input protein
        "echo_query": 1,  # see your input identifiers in the output
        "caller_identity": "www.awesome_app.org",  # your app name
    }
    request_url = "/".join([string_api_url, output_format, method])

    response = requests.post(request_url, data=params)
    data = json.loads(response.text)
    list_category = []
    for row in data:
        category = row["category"]
        list_category.append(category)
        # st.write("\t".join([term,preferred_names, str(fdr), description]))

    return list(set(list_category))


def get_cpds(conn_meta):
    sql_query = """
    SELECT cpd.*, gene.*, keggcpdgene.*, cpdgene.pubchemid
    FROM cpd
    INNER JOIN keggcpdgene ON cpd.keggid = keggcpdgene.keggid
    INNER JOIN gene ON keggcpdgene.geneid = gene.geneid
    INNER JOIN cpdgene ON cpd.pubchemid = cpdgene.pubchemid
    WHERE cpd.keggid IS NOT NULL
    """
    df_drug_meta = sql_df(sql_query, conn_meta)
    df_drug_meta = df_drug_meta.loc[:, ~df_drug_meta.columns.duplicated()].copy()
    df_drug_meta = df_drug_meta.drop_duplicates(subset=["pubchemid"]).reset_index(
        drop=True
    )

    return df_drug_meta


###########################################################################


def convert_df(df):
    return df.to_csv(index=False).encode("utf-8")


def get_col_colors(df):
    feat_cols = [col for col in df.columns if not col.startswith("Meta")]
    df_out = df[feat_cols]
    prefix_colors = {
        "ER": "red",
        "DNA": "blue",
        "RNA": "green",
        "AGP": "orange",
        "Mito": "pink",
        "mito": "pink",
    }
    col_colors = [
        prefix_colors.get(col.split("_")[0], "white") for col in df_out.columns
    ]
    if "name" not in df_out.columns:
        df_out["name"] = df["metacpdname"] + "_" + df["metasource"]
    df_plt = df_out.set_index("name")
    return df_plt, col_colors


def get_sql_kegg(table_name="keggcpdgene", col_name="geneid", list_geneid=["hdac6"]):
    where_clause = f" WHERE UPPER({table_name}.{col_name}) IN ({','.join(list_geneid)})"

    if table_name == "keggcpd":
        sql = f"SELECT * FROM keggcpd{where_clause}"
    else:
        sql = f"SELECT keggcpd.keggid, keggcpd.keggcpdname, {table_name}.{col_name} FROM keggcpd\
                INNER JOIN {table_name} ON {table_name}.keggid = keggcpd.keggid{where_clause}\
                GROUP BY keggcpd.keggcpdname, keggcpd.keggid, {table_name}.{col_name}"

    return sql


def get_sql_jump(table_name="cpdgene", col_name="geneid", list_geneid=["hdac6"]):
    select_clause = "SELECT cpd.pubchemid, cpdbatchs.batchid, cpd.synonyms, cpd.keggid, cpd.cpdname, cpd.smile"

    if table_name == "keggcpddis":
        sql = f"{select_clause}, keggcpddis.disid FROM cpdbatchs\
                INNER JOIN cpd ON cpdbatchs.pubchemid=cpd.pubchemid\
                INNER JOIN keggcpddis ON keggcpddis.keggid=cpd.keggid WHERE UPPER({table_name}.{col_name}) IN ({','.join(list_geneid)})"
    elif table_name in ["cpd", "cpdbatchs"]:
        sql = f"{select_clause} FROM cpdbatchs RIGHT JOIN cpd ON cpd.pubchemid=cpdbatchs.pubchemid WHERE UPPER({table_name}.{col_name}) IN ({','.join(list_geneid)})"
    else:
        sql = f"{select_clause}, {table_name}.{col_name} FROM cpdbatchs\
                INNER JOIN {table_name} ON {table_name}.pubchemid=cpdbatchs.pubchemid\
                INNER JOIN cpd ON cpdbatchs.pubchemid=cpd.pubchemid WHERE UPPER({table_name}.{col_name}) IN ({','.join(list_geneid)})\
                GROUP BY cpd.pubchemid, cpdbatchs.batchid, cpd.synonyms, cpd.keggid, cpd.cpdname, cpd.smile, {table_name}.{col_name}"

    return sql


########################################################################################


def find_sim_cpds(df1, df2):
    meta_cols = [
        "geneid",
        "name",
        "pubchemid",
        "keggid",
        "efficay",
        "smile",
        "source",
        "cpdname",
        "batchid",
    ]
    filter_col1 = [col for col in df1.columns if col not in meta_cols]
    filter_col2 = [col for col in df2.columns if col not in meta_cols]

    filter_col = list(set(filter_col1) & set(filter_col2))
    simi = cosine_similarity(df1[filter_col], df2[filter_col])
    return simi


def find_umap(df, title):
    meta_cols = [
        "geneid",
        "name",
        "pubchemid",
        "keggid",
        "efficay",
        "smile",
        "source",
        "cpdname",
        "batchid",
    ]
    filter_cols = [col for col in df.columns if col not in meta_cols]

    reducer = umap.UMAP(densmap=True, random_state=42, verbose=True)
    embedding = reducer.fit_transform(df[filter_cols])
    df_emb = pd.DataFrame({"x": embedding[:, 0], "y": embedding[:, 1]})
    df_emb["batchid"] = df["batchid"]
    fig = px.scatter(
        df_emb, x="x", y="y", hover_data=["batchid"], color="metasource", title=title
    )
    st.plotly_chart(fig, theme="streamlit", use_container_width=True)
