import networkx as nx
import igraph as ig
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
import requests
from pyvis.network import Network
import matplotlib
from streamlit_agraph import agraph, Node, Edge, Config
# from lib import streamlib

from streamlib import (
    get_list_category,
    get_stringDB_enr,
    int_to_str,
    sql_df,
    str_to_float,
)


st.set_page_config(layout="wide")


conn_meta = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"
conn_prof = "postgres://arno:12345@192.168.2.131:5432/ksilink_cpds"


def get_stringDB(df_all_umap, thresh=0.7, genecol="target"):
    string_api_url = "https://version-11-5.string-db.org/api"
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
    # st.write(request_url)
    results = requests.post(request_url, data=params)

    list_id0 = []
    list_id1 = []
    list_inter = []
    list_edges = []
    for line in results.text.strip().split("\n"):
        l = line.strip().split("\t")
        # st.write(l)
        if len(l) > 4:
            p1, p2 = l[2], l[3]

            # filter the interaction according to experimental score
            experimental_score = float(l[10])
            if experimental_score >= thresh:
                # print
                # print("\t".join([p1, p2, "experimentally confirmed (prob. %.3f)" % experimental_score]))
                list_id0.append(p1)
                list_id1.append(p2)
                if (p1, p2) not in list_edges:
                    list_edges.append((p1, p2))
                    list_inter.append(experimental_score)
    return list_edges, list_inter


def get_relation(df_genes):
    if not df_genes.empty:
        thres = 0
        list_edges = get_stringDB(df_genes, 0, "symbol")
        st.write(f"retrieved : {len(list_edges)} Interaction with {thres} threshold")

        st.write("## Computing Network")

        import networkx as nx

        H = nx.Graph(list_edges)
        import igraph as ig
        import leidenalg as la

        G = ig.Graph.from_networkx(H)
        partition = la.find_partition(G, la.ModularityVertexPartition)
        subg = partition.subgraphs()

        # thres_db = col_b.slider("cluster Thresholds", 2, 100,9,1)
        st.write(f"Total Number of clusters before Threshold: {len(subg)}")

        ################################### ENRICHMENT ##############################################
        list_cat = get_list_category(df_genes, "symbol")
        categ = st.selectbox("Select Category", list_cat)
        df_go_ento = get_stringDB_enr(df_genes, "symbol", categ)
        st.write(f"{categ} Enrichment", df_go_ento)
        df_go_ento["log_p_val"] = np.log10(df_go_ento["p_val"].apply(str_to_float))
        fig_bar = px.bar(df_go_ento, x="Description", y="log_p_val")
        st.plotly_chart(fig_bar)

        # subax1 = plt.subplot(121)
        st.write("### Displaying the full network/graph clustered with Leiden approach")
        # col1, col2,col3 = st.columns(3)

        vert_size = st.sidebar.slider(
            "vertex2 size", min_value=0.2, max_value=20.0, step=0.1, value=1.0
        )
        lab_size = st.sidebar.slider(
            "label2 size", min_value=0.2, max_value=20.0, step=0.1, value=1.0
        )
        box_size = st.sidebar.slider(
            "box2 size", min_value=400, max_value=1600, step=50, value=600
        )
        ig.config["plotting.backend"] = "matplotlib"
        subax1 = ig.plot(
            partition,
            vertex_label=partition.graph.vs["_nx_name"],
            vertex_label_size=lab_size,
            vertex_size=vert_size,
            margin=10,
            bbox=(0, 0, box_size, box_size),
        )
        # nx.draw(GG, with_labels=True,node_size=10,font_size=4)
        st.pyplot(subax1.figure, use_container_width=False)


@st.cache_data
def get_sample(df):
    st.write("onky 300 genes will be taken")
    df_genes = df.sample(300)
    return df_genes


################################### LOAD GENE CSV DATA ##############################################
for key in st.session_state.keys():
    del st.session_state[key]
# df_genes=get_cpds(conn_meta)
# unique_gene=df_genes['symbol'].unique()
# gene_col='symbol'
# st.write("TOTOT: ",len(df_genes['keggid']))
on = st.sidebar.toggle("Activate Data Upload")
if on:
    df_genes = pd.DataFrame()
    uploaded_file = st.file_uploader(
        "Choose csv file with a list of gene symbols  otherwise the reference 300  GeneMOA will be loaded",
        accept_multiple_files=False,
    )
    if uploaded_file:
        df_genes = pd.read_csv(uploaded_file)
        st.write("Loaded csv file:", df_genes)
        gene_col = st.text_input(
            "Write the column name with gene symbol", value="symbol"
        )
        if gene_col:
            if gene_col in df_genes.columns:
                df_genes.dropna(subset=gene_col, inplace=True)
                unique_gene = df_genes[gene_col].unique()
                st.write(f"you entered : {len(unique_gene)} genes")
                if len(unique_gene) > 300:
                    df_genes = get_sample(df_genes)
                    unique_gene = df_genes[gene_col].unique()
            else:
                st.warning(
                    "No columns symbol (i.e HDAC6,DHFR...) in your dataframe, check your file",
                    icon="⚠️",
                )
                exit(0)
            # get_relation(df_genes)
        # st.write(df_genes)
else:
    df_genes = pd.read_csv("Supp Data 3 - 300-Gene MoA sgRNA Library.csv")
    st.write("Loaded csv file:", df_genes.head())

    gene_col = st.text_input("Write the column name with gene symbols", value="symbol")
    unique_gene = df_genes[gene_col].unique()
    st.write(f"### You entered : {len(unique_gene)} genes")


# unique_gene=df_genes[gene_col].unique()
if not df_genes.empty:
    disp = st.sidebar.toggle("Display Data")
    if disp:
        st.write(df_genes)

    ################################### LOAD CPDS ##############################################
    choice = st.sidebar.radio(
        "Chose between cpds or crispr data", options=["Cpds", "CRISPR"]
    )
    df_inter = pd.DataFrame()
    if choice == "Cpds":
        st.warning("## Loading CPDS in development! ")
        # exit(0)

        b_list2 = unique_gene
        bq = []
        for bs in b_list2:
            if bs != "":
                bq.append("'" + bs + "'")

        sql_kegg = (
            "select cpdpath.pathid,keggcpdgene.geneid,gene.symbol, keggcpd.*, cpdbatchs.* from keggcpd\
                    left join keggcpdgene on keggcpd.keggid=keggcpdgene.keggid\
                    left join cpd on cpd.keggid=keggcpd.keggid \
                    left join cpdpath on cpdpath.pubchemid=cpd.pubchemid \
                    left join cpdbatchs on cpd.pubchemid=cpdbatchs.pubchemid \
                    left join gene on keggcpdgene.geneid=gene.geneid"
        )

        df_drug_meta = sql_df(sql_kegg, conn_meta)
        df_drug_meta = df_drug_meta.loc[:, ~df_drug_meta.columns.duplicated()].copy()
        # df_drug_meta = df_drug_meta.drop_duplicates(subset=["keggid"]).reset_index(
        # drop=True
        # )
        df_drug_meta = df_drug_meta[df_drug_meta["symbol"].isin(b_list2)]
        df_drug_meta.dropna(subset="symbol", axis=0, inplace=True)
        # disp2= st.toggle('Display Data')
        if disp:
            st.write("df_drug_meta", df_drug_meta)
            st.write(df_drug_meta.symbol.unique())

        # df_cpd_src = df_drug_meta[df_drug_meta["source"] == choix_source].reset_index(
        #     drop=True
        # )
        b_list2 = df_drug_meta["batchid"].to_list()
        bq = []
        for bs in b_list2:
            if bs is not None:
                bq.append("'" + bs + "'")

        sql_profile = (
            "select * from aggprofile where metabatchid  in (" + ",".join(bq) + ")"
        )

        df_cpd_prof = sql_df(sql_profile, conn_prof)
        list_sources = df_cpd_prof["metasource"].unique().tolist()
        choix_source = st.selectbox("Select the Source", list_sources)
        df_cpd_prof = df_cpd_prof[
            df_cpd_prof["metasource"] == choix_source
        ].reset_index(drop=True)
        if disp:
            st.write("df_cpd_prof", df_cpd_prof)

        # # st.write(df_cpd_src.sample(5))

        df_cpds_merge = df_cpd_prof.merge(
            df_drug_meta, left_on="metabatchid", right_on="batchid"
        ).reset_index(drop=True)

        # st.write("TOTOT: ",unique_gene)
        # df_inter = df_cpds_merge[df_cpds_merge[gene_col].isin(unique_gene)].reset_index(
        #     drop=True
        # )
        # df_inter = df_inter.drop_duplicates(subset=["keggid"]).reset_index(drop=True)
        df_inter = df_cpds_merge.drop_duplicates(subset=["keggid"]).reset_index(
            drop=True
        )
        if disp:
            st.write("Merged", df_cpds_merge)
            st.write("CDPS_GENES_COMMON", df_inter)
            st.write(" Unique GENES", df_inter["geneid"].unique())

    thres_gene = st.sidebar.slider("Min number of genes", 0, 100, 20, 1)
    ################################### LOAD CRISPR ##############################################
    if choice == "CRISPR":
        st.write("## Loading Crispr")

        bq = []
        for bs in unique_gene:
            bq.append("'" + bs + "'")
        sql_genes = (
            "select gene.*,genebatchs.batchid from gene  INNER JOIN genebatchs ON  genebatchs.geneid=gene.geneid where symbol  in ("
            + ",".join(bq)
            + ")  "
        )
        df_inter = sql_df(sql_genes, conn_meta)
        if not df_inter.empty:
            df_inter.drop_duplicates(inplace=True)
            df_inter["chromosome"] = df_inter["chromosome"].apply(int_to_str)
            st.write("GeneInfos with crispr profiles in Jump", df_inter)
            df_inter.reset_index(drop=True, inplace=True)

        else:
            st.warning(
                "Maybe you did not enter enough data or they are not symbols(i.e: HDAC6,DHFR)"
            )
            exit(0)
    if len(df_inter) < thres_gene:
        st.warning(
            "the intersection between your list and our data is pretty small, results will not be accurate or useful,\
                    you can fixed it by changing min number of genes or selecting CRISPR if not done yet",
            icon="⚠️",
        )
        exit(0)
    ################################### NETWORK ##############################################
    st.write("## Loading StringDB PPI")
    # col_a,col_b=st.columns(2)
    thres = st.sidebar.slider("Interaction Thresholds", 0.0, 1.0, 0.8, 0.02)
    if not df_inter.empty:
        list_edges, list_inters = get_stringDB(df_inter, thres, "symbol")
        st.write(f"retrieved : {len(list_edges)} Interaction with {thres} threshold")
        list_cat = get_list_category(df_inter, "symbol")
        categ = st.selectbox("Select Category", list_cat)
        df_go_ento = get_stringDB_enr(df_inter, "symbol", categ)

        st.write(f"{categ} Enrichment", df_go_ento)
        df_go_ento["log_p_val"] = -np.log10(df_go_ento["p_val"].apply(str_to_float))
        fig_bar = px.bar(df_go_ento, x="Description", y="log_p_val")
        st.plotly_chart(fig_bar)
        st.write("## Computing Network")
        H = nx.Graph(list_edges)
        G = ig.Graph.from_networkx(H)
        G = G.simplify()
        # partition = la.find_partition(
        #     G, la.ModularityVertexPartition, n_iterations=-1, weights=list_inters
        # )

        # subg = partition.subgraphs()
        communities = G.community_edge_betweenness()
        communities = communities.as_clustering()
        num_communities = len(communities)
        palette = ig.RainbowPalette(n=num_communities)
        list_gene = []
        list_clust = []
        cluster = 96
        list_df = []
        list_desc = []
        nodes = []
        # thres_db = col_b.slider("cluster Thresholds", 2, 100,9,1)
        net = Network(notebook=False)
        st.write(f"Total Number of clusters: {num_communities}")
        for i, g in enumerate(communities):
            cluster = cluster + 1
            G.vs[g]["color"] = i
            # st.write(G.vs[g].attributes())
            community_edges = G.es.select(_within=g)
            community_edges["color"] = i
            df_temp = get_stringDB_enr(G.vs[g]["_nx_name"], cat=categ)
            if not df_temp.empty:
                list_desc.append(str(df_temp["Description"].iloc[0]))
            else:
                list_desc.append("None")
            couleur = matplotlib.colors.to_hex(palette.get(i), keep_alpha=True)
            for u in G.vs[g]:
                # st.write(couleur)
                nodes.append(
                    Node(
                        id=u.index,
                        label=u["_nx_name"],
                        color=couleur,
                        title=list_desc[i],
                        shape="circle",
                        size=200,
                    )
                )
                # net.add_node(
                #     u.index,
                #     label=u["_nx_name"],
                #     color=couleur,
                #     shape="circle",
                #     title=u["_nx_name"],
                #     value=20,
                # )
                if not df_temp.empty:
                    list_df.append(str(df_temp["Description"].iloc[0]))
                    list_gene.append(u["_nx_name"])

                else:
                    list_df.append("Null")
                    list_gene.append(u["_nx_name"])
            # print(g.vs['_nx_name'])
            # if len(g.vs['_nx_name'])>thres_db:
        #     for name in G.vs[g]["_nx_name"]:
        #         list_clust.append(chr(cluster))
        #         list_gene.append(name)

        df_clust = pd.DataFrame()
        df_clust["symbol"] = list_gene
        df_clust["description"] = list_df
        # if df_clust.empty:
        #     st.warning("Please reduce the threshold to get more interaction")
        #     exit()
        # if disp:
        #     st.write("DF_CLUST", df_clust)
        df_umap_cluster = df_inter.merge(
            df_clust, left_on="symbol", right_on="symbol"
        ).reset_index(drop=True)
        # # df_umap_cluster['chromosome']=df_umap_cluster['chromosome'].apply(int_to_str)
        # if disp:
        #     st.write(df_umap_cluster)

        ################################### ENRICHMENT ##############################################
        A = G.get_edgelist()
        edges = [Edge(source=i, target=j, type="CURVE_SMOOTH") for (i, j) in A]
        config = Config(
            width=1200,
            height=1200,
            directed=False,
            nodeHighlightBehavior=True,
            highlightColor="#F7A7A6",
            collapsible=True,
            node={"labelProperty": "label"},
            physics=True,
            # link={"labelProperty": "label", "renderLabel": True},
        )
        st.write("### Displaying the full network/graph clustered with Leiden approach")
        return_value = agraph(nodes=nodes, edges=edges, config=config)
        # net.add_edges(A)
        # net.show("/mnt/shares/L/Temp/titi.html", notebook=False)

        # HtmlFile = open("/mnt/shares/L/Temp/titi.html", "r", encoding="utf-8")
        # source_code = HtmlFile.read()

        # components.html(source_code, height=1200, width=1200)

        # vert_size = st.sidebar.slider(
        #     "vertex size", min_value=0.2, max_value=20.0, step=0.1, value=1.0
        # )
        # lab_size = st.sidebar.slider(
        #     "label size", min_value=0.2, max_value=20.0, step=0.1, value=2.4
        # )
        # box_size = st.sidebar.slider(
        #     "box size", min_value=400, max_value=1600, step=50, value=600
        # )

        # fig, ax = plt.subplots()

        # ig.config["plotting.backend"] = "matplotlib"
        # ig.plot(
        #     communities,
        #     palette=palette,
        #     vertex_label=G.vs["_nx_name"],
        #     edge_width=1,
        #     target=ax,
        #     vertex_size=20,
        #     vertex_label_size=lab_size,
        # )
        # legend_handles = []
        # for i in range(num_communities):
        #     handle = ax.scatter(
        #         [],
        #         [],
        #         s=100,
        #         facecolor=palette.get(i),
        #         edgecolor="k",
        #         label=list_desc[i],
        #     )
        #     legend_handles.append(handle)
        # ax.legend(
        #     handles=legend_handles,
        #     title="Community:",
        #     bbox_to_anchor=(0, 1.0),
        #     bbox_transform=ax.transAxes,
        # )

        # st.pyplot(ax.figure, use_container_width=False)

        # nbr_of_cluster = ord(df_clust["cluster"].max()) - 96

        # comment for nothing

        st.write(
            f'## Computing UMAP with the main {num_communities} clusters and {len(df_umap_cluster["symbol"].unique())} genes'
        )
        st.write(
            "To increase or decrease number of clusters please change cluster threshold above"
        )

        numerics = ["float16", "float32", "float64"]
        # sql_dot="select * from umapemd"
        # df_temp=sql_df(sql_dot,conn_prof)
        # st.write(df_temp)
        # choix_source = "Ksilink_25"
        if choice == "Cpds":
            umap_sql = f"select * from umapemd where metasource='{choix_source}'"
        else:
            umap_sql = "select * from umapemd where metasource='CRISPER'"
        df_umap = sql_df(umap_sql, conn_prof)

        if choice == "Cpds":
            # df_umap = df_umap.dropna(subset="keggid", axis=0)
            # df_umap = df_umap.set_index("keggid")
            dict1 = df_umap_cluster.set_index("batchid").to_dict()["description"]
            dict2 = df_umap_cluster.set_index("batchid").to_dict()["symbol"]
            # st.write("Dict", dict1)
            df_umap["description"] = df_umap["metabatchid"].map(dict1)
            # df_umap["description"] = df_umap_cluster["description"]
            df_umap["target"] = df_umap["metabatchid"].map(dict2)
            df_umap["target"] = df_umap["target"].fillna("")
            df_umap_cluster = df_umap_cluster.reset_index()
            df_umap = df_umap.reset_index()
            df_umap["size"] = 5
            df_umap["size"] = df_umap["metabatchid"].apply(
                lambda x: 0.5 if x not in df_umap_cluster["batchid"].to_list() else 5
            )

            # dict1 = df_umap_cluster.set_index('keggid').to_dict()['cluster']
            # df_umap["cluster"] = df_umap['metakeggid'].map(dict1)
        else:
            dict1 = df_umap_cluster.set_index("symbol").to_dict()["description"]
            df_umap["description"] = df_umap["metaname"].map(dict1)
            df_umap["size"] = 5
            df_umap["size"] = df_umap["metaname"].apply(
                lambda x: 0.5 if x not in df_umap_cluster["symbol"].to_list() else 5
            )
            df_umap["target"] = df_umap["metaname"]
            df_umap["target"] = df_umap["target"].apply(
                lambda x: x if x in df_umap_cluster["symbol"].to_list() else ""
            )

        # df_umap["cluster"] = df_umap["cluster"].replace(list_enr)
        df_umap["description"] = df_umap["description"].fillna("others")
        # st.write(df_umap)
        ################################### SIMILARITY ##############################################
        # st.write("## Cosine Similarity")

        # cluster = df_umap
        # cluster = cluster[~cluster["cluster"].isin(["others", "Null"])]
        # cluster.reset_index(inplace=True, drop=True)
        # cluster["cluster"] = cluster["cluster"].str.split(" ").str[0]
        # col_dict = dict(
        #     zip(
        #         list(cluster.cluster.unique()),
        #         [
        #             tuple(int(c * 255) for c in cs)
        #             for cs in sns.color_palette("husl", len(cluster.cluster.unique()))
        #         ],
        #     )
        # )
        # cluster["color"] = cluster["cluster"].map(col_dict)
        # cluster["color"] = cluster["color"].fillna("yellow")
        # gene_color = cluster.set_index("metaname").to_dict()["color"]
        # sql_sim = "select * from crisprcos"
        # df_sim = sql_df(sql_sim, conn_prof)

        # df_sim.drop_duplicates(inplace=True)
        # df_sim.rename(
        #     columns={
        #         "symbol1": "Gene_1",
        #         "symbol2": "Gene_2",
        #         "sim": "Cosine_similarity",
        #     },
        #     inplace=True,
        # )
        # df_sim = df_sim[df_sim["Gene_1"].isin(cluster.metaname)]
        # df_sim = df_sim[df_sim["Gene_2"].isin(cluster.metaname)]
        # # df_sim=df_sim[df_sim['Cosine_similarity']>0.70]
        # df_sim.reset_index(inplace=True, drop=True)
        # final_df = (
        #     df_sim.merge(cluster, left_on="Gene_1", right_on="metaname", how="inner")
        #     .merge(cluster, left_on="Gene_2", right_on="metaname", how="inner")
        #     .drop(columns=["metaname_y", "metaname_x"])
        #     .rename(columns={"gene_group_name": "Gene_2_gene_group_name"})
        # )
        # # st.write(final_df)
        # interactions = final_df[["Gene_1", "Gene_2", "Cosine_similarity"]]
        # inter = np.array(final_df[["Gene_1", "Gene_2"]])
        # tup_inter = list(map(tuple, inter))
        # # st.write(tup_inter)
        # list_sim = final_df["Cosine_similarity"].to_list()
        # # G2 = nx.Graph()
        # G2 = nx.Graph(tup_inter)
        # IG = ig.Graph.from_networkx(G2)
        # partition2 = la.find_partition(
        #     IG, la.ModularityVertexPartition, n_iterations=-1
        # )
        # subax2 = ig.plot(
        #     partition2,
        #     vertex_label=partition2.graph.vs["_nx_name"],
        #     vertex_label_size=lab_size,
        #     vertex_size=vert_size,
        #     margin=10,
        #     bbox=(0, 0, box_size, box_size),
        # )
        # # nx.draw(GG, with_labels=True,node_size=10,font_size=4)
        # st.pyplot(subax2.figure, use_container_width=False)

        # G = nx.Graph(name="Protein Interaction Graph")
        # interactions = np.array(interactions)
        # # st.write(interactions)
        # for i in range(len(interactions)):
        #     interaction = interactions[i]
        #     a = interaction[0]  # protein a node
        #     b = interaction[1]  # protein b node
        #     w = float(
        #         interaction[2]
        #     )  # score as weighted edge where high scores = low weight
        #     G.add_weighted_edges_from([(a, b, w)])  # add weighted edge to graph

        # pos = nx.spring_layout(G)

        # # edges trace
        # edge_x = []
        # edge_y = []
        # xtext = []
        # ytext = []
        # str_sim = []

        # for edge in G.edges():
        #     x0, y0 = pos[edge[0]]
        #     x1, y1 = pos[edge[1]]
        #     edge_x.append(x0)
        #     edge_x.append(x1)
        #     edge_x.append(None)
        #     edge_y.append(y0)
        #     edge_y.append(y1)
        #     edge_y.append(None)
        #     tmp = final_df[
        #         (final_df["Gene_1"] == edge[0]) & (final_df["Gene_2"] == edge[1])
        #     ]
        #     if len(tmp) > 0:
        #         sim2 = tmp.Cosine_similarity.values[0]
        #         sim2 = float("{:.2f}".format(sim2))
        #         xtext.append((x0 + x1) / 2)
        #         ytext.append((y0 + y1) / 2)
        #         str_sim.append(sim2)

        # edge_trace = go.Scatter(
        #     x=edge_x,
        #     y=edge_y,
        #     line=dict(color="black", width=2),
        #     hoverinfo="none",
        #     showlegend=False,
        #     mode="lines",
        # )

        # # nodes trace
        # node_x = []
        # node_y = []
        # text = []
        # colors = []
        # for node in G.nodes():
        #     x, y = pos[node]
        #     node_x.append(x)
        #     node_y.append(y)

        #     text.append(node)
        #     colors.append(gene_color[node])

        # node_trace = go.Scatter(
        #     x=node_x,
        #     y=node_y,
        #     text=text,
        #     mode="markers+text",
        #     showlegend=False,
        #     hoverinfo="none",
        #     marker=dict(
        #         color=["rgb({}, {}, {})".format(r, g, b) for r, g, b in colors],
        #         size=50,
        #         line=dict(color="black", width=1),
        #     ),
        # )
        # # layout

        # # figure
        # structer_trace = go.Scatter(
        #     x=xtext,
        #     y=ytext,
        #     mode="text",
        #     hoverinfo="none",
        #     text=str_sim,
        #     textposition="top right",
        #     textfont=dict(color="blue", size=10),
        # )
        # layout = dict(
        #     autosize=True,
        #     template="plotly_white",
        #     width=1200,
        #     height=800,
        #     hovermode="closest",
        #     xaxis=dict(
        #         linecolor="white", showgrid=False, showticklabels=False, mirror=True
        #     ),
        #     yaxis=dict(
        #         linecolor="black", showgrid=False, showticklabels=False, mirror=True
        #     ),
        # )

        # fig = go.Figure(data=[edge_trace, node_trace, structer_trace], layout=layout)

        # st.plotly_chart(fig)

        if choice == "Cpds":
            # df_umap["keggid"] = df_umap_cluster['keggid']

            fig1 = px.scatter(
                df_umap,
                x="umap1",
                y="umap2",
                color="description",
                text="target",
                size="size",
                width=800,
                height=800,
                hover_data=["target", "metaname"],
            )
        else:
            fig1 = px.scatter(
                df_umap,
                x="umap1",
                y="umap2",
                color="description",
                text="target",
                size="size",
                width=800,
                height=800,
                hover_data=["target"],
            )
        st.plotly_chart(fig1, theme="streamlit", use_container_width=True)
        # st.write(df_umap)

    ################################### SIMILARITY ##############################################
    # st.write("## Similarity")
    # import seaborn as sns
    # df_umap_cluster.set_index('symbol',inplace=True)
    # lut = dict(zip(set(df_umap_cluster['cluster'].values), sns.hls_palette(len(set(df_umap_cluster['cluster'])), l=0.5, s=0.8)))
    # row_colors = pd.DataFrame(df_umap_cluster['cluster'].values)[0].map(lut)

    # sns.set(font_scale=0.6)

    # figsize = [20, 20]
    # fig_clusmap, ax1 = plt.subplots()
    # fig_clusmap = sns.clustermap(
    #             df_umap_cluster.select_dtypes(include=numerics),
    #             # metric="cosine",
    #             figsize=figsize,
    #             row_colors=row_colors.values,
    #         method="ward",
    #             xticklabels=False,
    #             yticklabels=True,
    #             col_cluster=False,
    #             cmap="vlag",
    #             center=0,
    #             vmin=-5,
    #             vmax=5,
    #             # ax=ax1
    #         )
