import psycopg2
import pandas as pd
import json
import requests
################################### CONNECTION DBs ##############################################
conn_meta = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksi_cpds",
    password="12345",
)

conn_prof = psycopg2.connect(
    host="192.168.2.131",
    port="5432",
    user="arno",
    database="ksilink_cpds",
    password="12345",
)

def int_to_str(ints):
    return str(ints)
def str_to_float(strs):
    return float(strs)

def sql_df(sql_str,conn=conn_meta):
    cur = conn.cursor()
    cur.execute(sql_str)
    rows = cur.fetchall()
    df_d = pd.DataFrame(rows, columns=[desc[0] for desc in cur.description])
    return df_d

################################### STRINGDB PPI ##############################################
def get_stringDB(df_all_umap,thresh=0.7,genecol='target'):    

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    params = {
        "identifiers" : "\r".join(df_all_umap[genecol].to_list()),
        # "identifiers" : "\r".join(["p53", "BRCA1", "cdk2", "Q99835"]), # your protein list
        "species" : 9606, # species NCBI identifier 
        "limit" : 1, # only one (best) identifier per input protein
        "echo_query" : 1, # see your input identifiers in the output
        "caller_identity" : "www.awesome_app.org" # your app name

    }

    request_url = "/".join([string_api_url, output_format, method])

    results = requests.post(request_url, data=params)

    list_id0=[]
    list_id1=[]
    list_inter=[]
    list_edges=[]
    for line in results.text.strip().split("\n"):
        l = line.strip().split("\t")
        p1, p2 = l[2], l[3]

        ## filter the interaction according to experimental score
        experimental_score = float(l[10])
        if experimental_score >=thresh:
            ## print 
            # print("\t".join([p1, p2, "experimentally confirmed (prob. %.3f)" % experimental_score]))
            list_id0.append(p1)
            list_id1.append(p2)
            list_edges.append((p1,p2))
            list_inter.append(experimental_score)
    return list(set(list_edges))

################################### STRING DB ENRICHMENT ##############################################
def get_stringDB_enr(df_all_umap,genecol='target',cat='KEGG',fdra=0.01):    

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "json"
    method = "enrichment" 
    params = {
        "identifiers" : "\r".join(df_all_umap[genecol].to_list()),
        # "identifiers" : "\r".join(["p53", "BRCA1", "cdk2", "Q99835"]), # your protein list
        "species" : 9606, # species NCBI identifier 
        "limit" : 1, # only one (best) identifier per input protein
        "echo_query" : 1, # see your input identifiers in the output
        "caller_identity" : "www.awesome_app.org" # your app name
    }
    request_url = "/".join([string_api_url, output_format, method])
    
    response = requests.post(request_url, data=params)
    data = json.loads(response.text)
    list_go=[]
    list_names=[]
    list_dec=[]
    list_fdr=[]
    list_category=[]
    list_p_value=[]
    for row in data:
        term = row["term"]
        preferred_names = ",".join(row["preferredNames"])
        fdr = float(row["fdr"])
        description = row["description"]
        category = row["category"]
        p_value  = row['p_value']

        if category == cat and fdr < fdra:
            list_go.append(term)
            list_names.append(preferred_names)
            list_fdr.append(format(fdr,"5.2e"))
            list_p_value.append(format(p_value,"5.2e"))
            list_dec.append(description)
            list_category.append(category)
            # st.write("\t".join([term,preferred_names, str(fdr), description]))
    df_go=pd.DataFrame()
    df_go['GO Term']= list_go
    df_go['Preferred Names']= list_names
    df_go['Description']=list_dec
    df_go['fdr']=list_fdr
    df_go['p_val']=list_p_value
    df_go['category']=list_category
    return df_go       


def get_list_category(df_all_umap,genecol='target'):
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "json"
    method = "enrichment" 
    params = {
        "identifiers" : "\r".join(df_all_umap[genecol].to_list()),
        # "identifiers" : "\r".join(["p53", "BRCA1", "cdk2", "Q99835"]), # your protein list
        "species" : 9606, # species NCBI identifier 
        "limit" : 1, # only one (best) identifier per input protein
        "echo_query" : 1, # see your input identifiers in the output
        "caller_identity" : "www.awesome_app.org" # your app name
    }
    request_url = "/".join([string_api_url, output_format, method])
    
    response = requests.post(request_url, data=params)
    data = json.loads(response.text)   
    list_category=[]   
    for row in data:       
        category = row["category"]       
        list_category.append(category)
            # st.write("\t".join([term,preferred_names, str(fdr), description]))

    return list(set(list_category))       
    
def get_cpds():
   
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
    df_drug_meta = df_drug_meta.drop_duplicates(subset=["pubchemid"]).reset_index(drop=True)
         
   
    return df_drug_meta  

