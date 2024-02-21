import streamlit as st

# from streamlit_plotly_events import plotly_events
import numpy as np
import plotly.express as px

from streamlib import sql_df
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import cosine_similarity
import streamlit.components.v1 as components
import networkx as nx
import igraph as ig
import leidenalg as la
import matplotlib.pyplot as plt
import psycopg2
# import webbrowser
# webbrowser.open_new_tab("http://librenms.ksilink.int/overview?dashboard=1")

# components.iframe("http://librenms.ksilink.int/",height=1600)
# https://www.google.com/

conn_meta = "postgres://arno:123456@192.168.2.131:5432/ksi_cpds"
conn_prof = "postgres://arno:12345@192.168.2.131:5432/ksilink_cpds"

# st.title('Experimental: works only on 131')
# sql_rep='select symbol1 from crisprcos'
# df_rep= sql_df(sql_rep,conn_prof)
# df_rep = df_rep.drop_duplicates().reset_index(drop=True)

# df_rep=pd.read_csv('Supp Data 3 - 300-Gene MoA sgRNA Library.csv',usecols=['symbol'])
# df_rep= df_rep.append({'symbol':'DMSO'},ignore_index=True)
# # st.write(df_rep)
# # bq = []
# # for bs in df_rep['symbol1']:
# #     bq.append("'" + bs.upper() + "'")


# st.subheader("added metasource to cpdbatch", divider="rainbow")

# sql_meta = "select * from cpdbatchs"
# df_meta = sql_df(sql_meta, conn_meta)
# list_cpd = [f for f in df_meta.batchid.unique() if "JCP" not in f]
# st.write(
#     "CPD_BATCHS",
#     df_meta[df_meta["batchid"].isin(list_cpd)].reset_index(drop=True),
# )
# st.write("df_CPD_BATCHS", df_meta)


st.subheader(
    "Confusion when creating DB in 8_create_DB... between Metadata_JCP2022 and Metadata_broad_sample",
    divider="rainbow",
)

col1, col2, col3 = st.columns(3)
source = "source_7_625"
df = pd.read_feather(
    f"/mnt/shares/L/PROJECTS/JUMP-CP/jumpAWS/concatSources/profiles/combat_{source}.fth"
)
# st.write(df.columns)
st.write(df[df["Metadata_broad_sample"] == "S8808-01"])

# # st.header(
# #     "Inconsistency batchid between Umap.batchid and addprofile.batchid",
# #     divider="rainbow",
# # )
# col1.write("from the feather file")
# col1.write(df)


# # # sql_umqpemd =  f"SELECT * FROM aggprofile where metasource='CRISPER' and metabatchid  in (" + ",".join(bq) + ") "
sql_umqpemd = "SELECT * FROM aggprofile where metasource='Ksilink_25'"
df_src_agg = sql_df(sql_umqpemd, conn_prof)
col2.write("Aggregation profiles")
col2.write(df_src_agg.metabatchid)

sql_umqpemd = "SELECT * FROM umapemd where metasource='Ksilink_25'"
df_src_emd = sql_df(sql_umqpemd, conn_prof)
col3.write("Umap")
col3.write(df_src_emd.metabatchid)

# df_keep = pd.read_csv("df_keep_prof.csv")
# st.write(df_keep)

df_inter = df_src_emd[~df_src_emd["metabatchid"].isin(df_src_agg["metabatchid"])]
st.subheader(
    "In total more than 5000 ids are not common between agg and umap for Ksilink",
    divider="rainbow",
)
st.write("Not in agg", df_inter.metabatchid)

st.subheader(
    "Test tables in meta, print unique batchid and total batchids", divider="rainbow"
)
conn = psycopg2.connect(conn_meta)
cursor = conn.cursor()
cursor.execute(
    "select table_name from information_schema.tables where table_schema='public'"
)
tables = [i[0] for i in cursor.fetchall()]
# st.write(tables)
for table in tables:
    sql_req = f"select * from {table}"
    # st.write(table)
    df_tab = sql_df(sql_req, conn_meta)
    if "batchid" in df_tab.columns:
        list_cpd = [f for f in df_tab.batchid.unique() if "JCP" not in f]
        st.write(
            f"{table} -> batchid:",
            (len(df_tab.batchid.unique()), len(df_tab.batchid), len(list_cpd)),
        )
    if "pubchemid" in df_tab.columns:
        st.write(
            f"{table} -> pubchemid:",
            (len(df_tab.pubchemid.unique()), len(df_tab.pubchemid)),
        )
    if "pubchemsid" in df_tab.columns:
        st.write(
            f"{table} -> pubchemsid:",
            (len(df_tab.pubchemsid.unique()), len(df_tab.pubchemsid)),
        )


st.subheader(
    "Test tables in prof, print unique batchid and total batchids", divider="rainbow"
)
conn = psycopg2.connect(conn_prof)
cursor = conn.cursor()
cursor.execute(
    "select table_name from information_schema.tables where table_schema='public'"
)
tables = [i[0] for i in cursor.fetchall()]
# st.write(tables)
for table in tables:
    sql_req = f"select * from {table}"
    # st.write(table)
    df_tab = sql_df(sql_req, conn_prof)
    # if "batchid" in df_tab.columns:
    #     list_cpd = [f for f in df_tab.batchid.unique() if "JCP" not in f]
    #     st.write(
    #         f"{table} -> batchid:",
    #         (len(df_tab.batchid.unique()), len(df_tab.batchid), len(list_cpd)),
    #     )
    if "metabatchid" in df_tab.columns:
        list_cpd = [f for f in df_tab.metabatchid.unique() if "JCP" not in f]
        st.write(
            f"{table} -> metabatchid:",
            (len(df_tab.metabatchid.unique()), len(df_tab.metabatchid), len(list_cpd)),
        )
    if "pubchemid" in df_tab.columns:
        st.write(
            f"{table} -> pubchemid:",
            (len(df_tab.pubchemid.unique()), len(df_tab.pubchemid)),
        )


st.subheader("Removed pubchemid column", divider="rainbow")
sql_kegg_cpd = "select metabatchid from aggprofile where metasource='Ksilink_625'"
df_agg = sql_df(sql_kegg_cpd, conn_prof)
# st.write("agg_batchid", df_keggcpd)
sql_kegg_cpd = "select batchid from cpdbatchs where metasource='Ksilink_625'"
df_batch = sql_df(sql_kegg_cpd, conn_meta)

cola, colb = st.columns(2)
df = df_agg[~df_agg["metabatchid"].isin(df_batch.batchid)]
cola.write(df.reset_index(drop=True))
df2 = df_batch[~df_batch["batchid"].isin(df_agg.metabatchid)]
colb.write(df2.reset_index(drop=True))
exit(0)
# df_sel = df_src_emd[df_src_emd["metabatchid"]=='DMSO']
# sim_crispr = find_sim_cpds(df_src_emd, df_sel)
# df_hist_cpd = pd.DataFrame(
#     {"sim": sim_crispr.flatten().tolist(
#     ), "metabatchid": df_src_emd["metabatchid"]}
# )
# df_hist_cpd['metabatchid']=df_hist_cpd['metabatchid'].str.split('_').str[0]
# df_hist_cpd=df_hist_cpd[df_hist_cpd['metabatchid'].isin(df_rep['symbol'])]

# df_rna=pd.read_csv('Recursion_U2OS_expression_data.csv')
# dict_rna = df_rna.set_index('gene').to_dict()['zfpkm']

# df_hist_cpd['zfpkm']=df_hist_cpd['metabatchid'].map(dict_rna)
# df_hist_cpd.dropna(subset='zfpkm', inplace=True)
# df_hist_cpd=df_hist_cpd[df_hist_cpd['zfpkm']>-30]

# df_hist_cpd['expressed']='No'
# df_hist_cpd.loc[df_hist_cpd['zfpkm'] > -3, 'expressed'] = 'Yes'


# #
# # df_hist_cpd=df_hist_cpd[df_hist_cpd['zfpkm']<-3]
# # st.write(df_hist_cpd)
# fig = px.scatter(df_hist_cpd,x='sim',y='zfpkm',hover_data=['metabatchid'],color='expressed')
# st.plotly_chart(fig)

# df_rna2=pd.read_csv('src13_u2os.csv')
# df_rna2['zfpkm']=df_rna2['symbol'].map(dict_rna)

# df_rna2 = df_rna2.dropna(axis=1)
# df_rna2=df_rna2[df_rna2['zfpkm']>-30]
# df_rna2=df_rna2[df_rna2['TPM']>3]
# df_rna2['log10TPM']=np.log10(df_rna2['TPM'])
# fig2= px.scatter(df_rna2,x='log10TPM',y='zfpkm',hover_data=['symbol'],trendline='ols', trendline_color_override = 'red')
# st.plotly_chart(fig2)
# model = px.get_trendline_results(fig2)
# alpha = model.iloc[0]["px_fit_results"].params[0]
# beta = model.iloc[0]["px_fit_results"].params[1]
# st.write(f'alpha={alpha}, beta={beta}')
# from tdc.resource import PrimeKG

# data = PrimeKG(path="../data")
# drug_feature = data.get_features(feature_type="drug")
# data.to_nx()
# st.write(data.get_node_list("disease"))
st.set_page_config(
    layout="wide",
)


@st.cache_data
def load_data():
    primekg = pd.read_csv("../kg.csv", low_memory=False)
    return primekg


vert_size = st.sidebar.slider(
    "vertex size", min_value=0.2, max_value=20.0, step=0.1, value=1.0
)
lab_size = st.sidebar.slider(
    "label size", min_value=0.2, max_value=20.0, step=0.1, value=2.4
)

KG = load_data()
# sel_relation = st.selectbox("Select Relation", KG.relation.unique())
# b = KG[KG["relation"] == sel_relation]
list_gene = ["HDAC6", "TTN", "CALM3"]


var_text = st.text_area("Enter your genes list", help="Name or ID separated by enter")
var_t = var_text.split("\n")
list_gene = [t.strip().upper() for t in var_t]

# KG = KG.query('x_name == "HDAC6"|y_name=="TTN"')
KG = KG.query("x_name == @list_gene | y_name==@list_gene")
# st.write(b.sample(20))
# b = primekg.sample(200)
# st.write(b[b.relation == "protein_protein"][["x_name", "y_name"]].values)
graph_list = []
# ig.config["plotting.backend"] = "matplotlib"
options = {
    "node_color": "blue",
    "node_size": vert_size,
    "width": 1,
    "font_size": lab_size,
    "edge_color": "green",
    "alpha": 0.6,
}
for i in KG.relation.unique():
    # if i == "protein_protein":
    G = nx.Graph()
    G.add_edges_from(KG[KG.relation == i][["x_name", "y_name"]].values, relation=i)
    # H = ig.Graph.from_networkx(G)
    # st.write("Relation", i)
    # st.write("Rel", KG[KG.relation == i].reset_index())
    # list_path = [p for p in nx.all_shortest_paths(G, source="HDAC6", target="TTN")]
    # st.write(list(nx.connected_components(G)))
    # GG = nx.Graph()
    # for p in list_path:
    #     nx.add_path(GG, p)

    # H2 = ig.Graph.from_networkx(GG)
    # st.pyplot(
    #     ig.plot(
    #         H,
    #         vertex_label=H.vs["_nx_name"],
    #         vertex_label_size=lab_size,
    #         vertex_size=vert_size,
    #         edge_width=0.1,
    #     ).figure,
    #     use_container_width=True,
    #     clear_figure=True,
    # )
    fig, ax = plt.subplots()
    ax.set_title(i)
    nx.draw(G, with_labels=True, **options)
    st.pyplot(
        fig,
        use_container_width=True,
        clear_figure=True,
    )
    G.clear()
    # H.clear()
# subax1 =
# nx.draw(GG, with_labels=True,node_size=10,font_size=4)

# st.plotly_chart(subax1)
# a = primekg.query('y_name=="HDAC6"|x_name=="HDAC6"')
# st.write(a)
t = st.toggle("Compute")
if t:
    import os
    import cv2
    import numpy as np
    from joblib import Parallel, delayed
    import torch
    from torchvision import transforms
    from torch.autograd import Variable
    from torch import nn
    import torchvision.transforms as T
    import multiprocessing
    from torchvision.models import efficientnet_b0, EfficientNet_B0_Weights
    import gc

    m_cpu = multiprocessing.cpu_count()

    torch.cuda.empty_cache()
    # exit(0)
    # Load the model
    # model = torch.load('/mnt/shares/L/PROJECTS/JUMP-CRISPR/Code/streamlit_arno/streamlit-1/dinov2_vits14_linear4_head.pth')
    # 'cuda' if torch.cuda.is_available() else
    dinov2_vits14 = torch.hub.load("facebookresearch/dinov2", "dinov2_vits14")
    num_gpus = torch.cuda.device_count()
    st.write(num_gpus)
    if m_cpu > 30:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        # model=nn.DataParallel(dinov2_vits14,device_ids=list(range(num_gpus)))
        # st.write(dir(model))
        # st.write(help(model.prediction))
        weights = EfficientNet_B0_Weights.DEFAULT
        model = nn.DataParallel(efficientnet_b0(weights=weights))
        # preprocess = weights.transforms()

        # model=efficientnet_b0(EfficientNet_B0_Weights)
    else:
        device = torch.device("cpu")
        model = dinov2_vits14

    model.to(device)
    model.eval()
    # st.write(help(model.features))
    # exit(0)
    transform_image = T.Compose(
        [T.ToTensor(), T.Resize(244), T.CenterCrop(224), T.Normalize([0.5], [0.5])]
    )
    # model.eval()

    # Create a new model that outputs the features from the second to last layer
    # feature_model = nn.Sequential(*list(model.children())[:-1])

    def process_image(file1):
        # Construct the file names for the other channels
        l = file1.replace("\\", "/").replace(".tif", "").split("/")[-1].split("_")
        Plate = l[0]
        Well = l[1]
        # st.write(Well)
        file2 = file1.replace("C01.tif", "C02.tif")
        file2 = file2.replace("A01", "A02")
        file3 = file1.replace("C01.tif", "C03.tif")
        file3 = file2.replace("A01", "A03")
        # st.write(file1)

        # Read the images
        img1 = cv2.imread(file1, -1)
        img2 = cv2.imread(file2, -1)
        img3 = cv2.imread(file3, -1)
        # st.write(img2)
        # Convert to grayscale for processing
        # gray = cv2.cvtColor(img1, cv2.COLOR_BGR2GRAY)

        # Thresholding and morphological operations to isolate nuclei
        _, thresh = cv2.threshold(img1, 250, 65535, cv2.THRESH_BINARY)
        thresh = cv2.convertScaleAbs(thresh)
        nb_blobs, im_with_separated_blobs, stats, _ = cv2.connectedComponentsWithStats(
            thresh
        )
        sizes = stats[:, -1]

        sizes = sizes[1:]
        nb_blobs -= 1
        min_size = 2000

        # output image with only the kept components
        im_result = np.zeros_like(im_with_separated_blobs)
        # for every component in the image, keep it only if it's above min_size
        for blob in range(nb_blobs):
            if sizes[blob] >= min_size:
                # see description of im_with_separated_blobs above
                im_result[im_with_separated_blobs == blob + 1] = 255
        im_result = cv2.convertScaleAbs(im_result)
        contours, _ = cv2.findContours(
            im_result, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
        )
        # st.image(im_result)
        # exit(0)
        # all_embeddings ={}
        # emb=[]
        # list_plate=[]
        # list_well=[]
        # list_score=[]
        # list_cat=[]
        nuclei = []
        for i, contour in enumerate(contours):
            # Get bounding box coordinates
            x, y, w, h = cv2.boundingRect(contour)
            h2 = 22
            w2 = 22
            # Crop the nucleus from each channel
            if (
                y - h > 0
                and x - w > 0
                and y + h < img1.shape[0]
                and x + w < img1.shape[1]
            ):
                nucleus1 = img1[y - h2 : y + h + h2, x - w2 : x + w + w2]
                nucleus2 = img2[y - h2 : y + h + h2, x - w2 : x + w + w2]
                nucleus3 = img3[y - h2 : y + h + h2, x - w2 : x + w + w2]
                # st.write((h,w))
                # Stack nuclei along the third axis to create a multi-channel nucleusq
                nucleus = np.dstack([nucleus3, nucleus3, nucleus3])

                # Resize to match the input size expected by the model
                nucleus = cv2.resize(
                    nucleus, (224, 224)
                )  # Replace with your input size
                # st.image(nucleus,clamp=True)
                # Normalize the image
                nucleus = np.float32(
                    nucleus / np.max(nucleus)
                )  # Use 65535 for 16-bit images

                # Convert to PyTorch tensor and add an extra dimension
                transform = transforms.Compose([transforms.ToTensor()])
                # nucleus  = transform_image(np.float32(nucleus)).unsqueeze(0)
                # st.write(nucleus)
                nuclei.append(transform(np.float32(nucleus)).unsqueeze(0))
        #         # batch = preprocess((nucleus)).unsqueeze(0)
        #         # st.write(batch)
        #         # temp=model(batch.to(device)).squeeze(0)
        #         # emb.append(temp)
        #         prediction = model(nucleus.to(device)).squeeze(0).softmax(0)
        #         class_id = prediction.argmax().item()
        #         score = prediction[class_id].item()
        #         category_name = weights.meta["categories"][class_id]
        #         list_score.append(score)
        #         list_cat.append(category_name)
        #         # # nucleus = transform(nucleus)
        #         # nucleus = nucleus.unsqueeze(0)
        #         list_plate.append(Plate)
        #         list_well.append(Well)
        #         # Extract features
        #         # features = model(Variable(nucleus))
        #         with torch.no_grad():
        #             embeddings = model.features(nucleus.to(device))
        #             embed = nn.AvgPool2d(7)
        #             toto =embed(embeddings)
        #             # all_embeddings[file1+'_'+str(i)] = np.array(embeddings[0].cpu().numpy()).reshape(1, -1).tolist()
        #             temp = np.squeeze(np.array(toto[0].cpu().numpy())).tolist()
        #             # temp2 = np.squeeze(np.array(temp).reshape(-1, 384))
        #             # st.write(toto.shape)

        #             # flat_list = [item for sublist in temp2 for item in sublist]
        #             emb.append(temp)
        #     # Process the features as needed
        # ...q
        # flat_list = [item for sublist in emb for item in sublist]q
        # df2=pd.DataFrame()
        # df2['Plate']=list_plate
        # df2['Well']=list_well
        # df2['Score']=list_score
        # df2['Category']=list_cat
        # del model

        return nuclei

    def gpu_comp(nuc):
        # Convert to PyTorch tensor and add an extra dimension
        # transform = transforms.Compose([transforms.ToTensor()])
        # # nucleus  = transform_image(np.float32(nucleus)).unsqueeze(0)
        # # st.write(nucleus)
        # nucleus=transform(np.float32(nuc))
        # batch = preprocess((nucleus)).unsqueeze(0)
        # st.write(batch)
        # temp=model(batch.to(device)).squeeze(0)
        # emb.append(temp)
        prediction = model(nuc.to(device)).softmax(0)
        class_id = prediction.argmax(dim=1)
        # st.write(class_id.shape)
        # score = prediction[class_id].item()
        category_name = [weights.meta["categories"][w] for w in class_id]
        # category_name = weights.meta["categories"][class_id]
        return category_name

    def embed(nuc):
        with torch.no_grad():
            embeddings = model.module.features(nuc.to(device))
            embed = nn.AvgPool2d(7)
            toto = embed(embeddings).squeeze()
            temp = np.array(toto.cpu().numpy())
            return temp

    import time

    start = time.time()
    # Get list of all C01.tif files (nuclei)
    path_dcm = "/mnt/shares/O/BTSData/MeasurementData/DCM Exp31 LDC 4doses DRC/DCM-LDC-PL03c_20231119_041305/DCM-LDC-PL03c"
    path_test = "/mnt/shares/L/Temp/DC/"
    # tif_files = [f for f in os.listdir(path_dcm) if f.endswith('C01.tif')]
    tif_files = [f for f in os.listdir(path_dcm) if f.endswith("C01.tif")]
    # img1 = cv2.imread(tif_files[0], -1)
    # st.write(tif_files)
    list_df = Parallel(n_jobs=-1, prefer="threads", max_nbytes=None)(
        delayed(process_image)(f"{path_dcm}/{file}") for file in tif_files
    )
    del tif_files
    # st.write(list_df[0])
    flat_list = [item for sublist in list_df for item in sublist]
    # st.write(len(flat_list))
    # del list_df
    b = np.squeeze(torch.stack(flat_list))

    # st.write(b.shape)
    del flat_list
    g = np.array_split(b, 1000)
    del b
    list_pred = []
    list_emb = []
    # c=torch.tensor(0)
    for i, item in enumerate(g):
        list_pred.extend(gpu_comp(item))
        list_emb.append(embed(item))
        # st.write('remaining: ', len(g)-i)
    end = time.time()
    st.write("Total Exec Time= ", (end - start) / 60)

    out = np.concatenate(list_emb, axis=0)

    # st.write(len(list_pred))
    # st.write(out.shape)
    # b=[]
    # for i,item in enumerate(list_pred):
    #     b.append(np.squeeze(list_pred[i]))
    #     # st.write(b.shape)
    #     # list_pred = gpu_comp(b)
    #     # st.write(len(list_pred))
    # st.write(b)
    # feature=[process_image(f"{path_test}/{file}") for file in tif_files]
    # list_df=[]
    # for i,file in enumerate(tif_files):
    #     feat=process_image(f"{path_test}/{file}")
    #     list_df.append(feat)
    #     st.write("Remaining Images:",len(tif_files)-i)
    # # embedding_list = list(feature.values())
    # flat_list = [item for sublist in feature for item in sublist]
    # exit(0)
    import umap

    # embed = np.array(embedding_list).reshape(-1, 384)
    # df=pd.concat(list_df,ignore_index=True)
    # df['Plate']
    # st.write(df)
    # st.write(df['Category'].value_counts())
    # exit(0)
    numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
    emb = umap.UMAP(random_state=42, verbose=False).fit_transform(out)
    df_all_umap = pd.DataFrame()
    df_all_umap["X"] = emb[:, 0]
    df_all_umap["Y"] = emb[:, 1]
    # df_all_umap['Plate']=df['Plate']
    # df_all_umap['Well']=df['Well']
    df_all_umap["Category"] = list_pred
    fig = px.scatter(df_all_umap, x="X", y="Y", color="Category")
    st.plotly_chart(fig, theme="streamlit", use_container_width=True)
    del model
    gc.collect()
    # torch.cuda.
# st.write(df)
# exit(0)

# def toto():
#     sql_umqpemd =  f"SELECT * FROM aggprofile where metasource='CRISPER'"
#     df_src_emd = sql_df(sql_umqpemd, conn_prof)
#     return df_src_emd


# df_dcp= pd.read_csv('/mnt/shares/L/Temp/Test/test.csv')
# # numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
# # cols_num = df_dcp.select_dtypes(include=numerics).columns
# # st.dataframe(df_dcp[cols_num].style.format('{:.4g}'))
# # st.write(df_dcp.dtypes)
# st.dataframe(df_dcp,column_config={"EC_50_Nuclei_Tot":st.column_config.NumberColumn("EC_50_Nuclei_Tot", format="%7.2g",)})
# numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
# cols_num = df_dcp.select_dtypes(include=numerics).columns
# for col in cols_num:
#     df_dcp[col] = df_dcp[col].apply(lambda x: "{:.2e}".format(x))


# st.dataframe(df_dcp)

# df=toto()
# all_cpds=["Staurosporine","Saccharin","Sorbitol",   "CA-074","BAYK8644",
#     "Lys05","Cucurbitacin","FCCP","Rapamycin","Cladribine","Cytarabine",
#     "Etoposide","Berberine","Fluphenazine","Docetaxel","Oxibendazole",
#     "Ethoxyquin","Rotenone","GSK2879552","BAY K8644","NMDA","Tetrandrine",
#     'Aloxistatin','Dexamethasone','Quinidine','LY2109761','AMG900','NVS-PAK1-1','Mirk-IN-1','Daporinad']
# df_dcps=pd.DataFrame()
# df_dcps['cpds']=all_cpds
# df_dcps.to_csv('/mnt/shares/L/Temp/cpds_tox.csv',index=None)
# cols = [c for c in df.columns if  not c.startswith("meta")]
# meta_cols = [c for c in df.columns if  c.startswith("meta")]
# X = df[cols]
# # KNN
# nb_cluster=st.slider('Number of clusters',min_value=2,max_value=30,value=15,step=1)
# #kmeans = KMeans(n_clusters=nb_cluster, random_state=0).fit(X)


# ## PCA
# # from sklearn.decomposition import PCA
# from cuml import PCA
# pca_2d = PCA(n_components=2)
# projection_2d = pca_2d.fit_transform(X)
# # st.write(projection_2d[0])
# emd_pca = pd.DataFrame()
# emd_pca["pca1"] = projection_2d[0]
# emd_pca["pca2"] = projection_2d[1]
# kmeans = KMeans(n_clusters=nb_cluster, random_state=0).fit(emd_pca)
# # emd_pca[meta_cols] = df[meta_cols]

# emd_pca['Cluster']=kmeans.labels_
# emd_pca['Cluster'] = emd_pca['Cluster'].astype(str)
# fig1 = px.scatter(
#                 emd_pca,
#                 x="pca1",
#                 y="pca2",
#                 color="Cluster",


#                 title=f" PCA ",
#                 # hover_data=["metabatchid"],
#             )

# st.plotly_chart(fig1, theme="streamlit",
#                                 use_container_width=True)

# ## UMAP
# # import umap
# from cuml import UMAP
# umap_2d = UMAP(n_components=2, n_neighbors=30,
#                 min_dist=0, verbose=True, random_state=42)
# mapper=umap_2d.fit(X)
# projection_2d = mapper.transform(X)
# emd_umap = pd.DataFrame()
# emd_umap["umap1"] = projection_2d[0]
# emd_umap["umap2"] = projection_2d[1]
# # emd_umap[meta_cols] = df[meta_cols]
# emd_umap['Cluster']=kmeans.labels_
# emd_umap['Cluster'] = emd_umap['Cluster'].astype(str)
# fig = px.scatter(
#                 emd_umap,
#                 x="umap1",
#                 y="umap2",
#                 color="Cluster",


#                 title=f" UMAP ",
#                 # hover_data=["metabatchid"],
#             )

# st.plotly_chart(fig, theme="streamlit",
#                                 use_container_width=True)
# from cuml.metrics import trustworthiness
# cu_score=trustworthiness(X,projection_2d)
# st.write(" cuml's trustworthiness score : ",cu_score)

# st.write(df)
# import matplotlib.pyplot as plt
# import plotly.figure_factory as ff
# df_g1=emd_pca.groupby("Cluster").mean()
# df_g2=emd_umap.groupby("Cluster").mean()
# st.write(df_g1)
# st.write(df_g2)
# fig3=ff.create_quiver( df_g1["pca1"],df_g1["pca2"],(df_g2["umap1"]-df_g1["pca1"]),df_g2["umap2"]-df_g1["pca2"],scale=1)
# st.plotly_chart(fig3)
# import umap.plot
# fig6 = umap.plot.diagnostic(mapper, diagnostic_type='pca')
# # st.write(fig6)
# st.pyplot(fig6.axes.figure)

# fig7=umap.plot.diagnostic(mapper, diagnostic_type='vq')
# st.pyplot(fig7.axes.figure)

# fig8=umap.plot.diagnostic(mapper, diagnostic_type='local_dim')
# st.pyplot(fig8.axes.figure)

# fig9=umap.plot.diagnostic(mapper, diagnostic_type='neighborhood')
# st.pyplot(fig9.axes.figure)

# fig10 = umap.plot.connectivity(mapper, edge_bundling='hammer')
# st.pyplot(fig10.axes.figure)
# X=df_g1.values

# grid_unit_size = .02
# margin = 0.25

# x_min = X[:, 0].min() - margin
# x_max = X[:, 0].max() + margin
# y_min = X[:, 1].min() - margin
# y_max = X[:, 1].max() + margin

# x_range = np.arange(start = x_min, stop = x_max, step = grid_unit_size )
# y_range = np.arange(start = y_min, stop = y_max, step = grid_unit_size )

# x_gridvalues, y_gridvalues = np.meshgrid(x_range, y_range)

# # COMPUTE PROBABILITIES ON GRID
# gridvalues_combined_tidy = np.vstack([x_gridvalues.flatten(), y_gridvalues.flatten()]).T
# # knn_class_probabilities = knn_classifier.predict_proba(gridvalues_combined_tidy)

# probability_postive_class = knn_class_probabilities[:,1]

# fig = go.Figure(data=[
#     go.Contour(
#         x = x_range
#         ,y = y_range
#         ,z = probability_postive_class.reshape(x_gridvalues.shape)
#         ,colorscale = 'RdBu'
#         #,alpha = .3
#     )
# ])

# x = np.linspace(-4, 4, 9)
# y = np.linspace(-5, 5, 11)
# random_data = np.random.random((11, 9))
# x_1, y_1 = np.meshgrid(df_g1["pca1"], df_g1["pca2"])
# data1,data2 = np.meshgrid(df_g2["umap1"],df_g2["umap2"])

# # data= np.concatenate(df_g2["umap1"].values,df_g2["umap2"].values)
# st.write(data1)
# fig6 = plt.plot_surface(x_1,y_1,data1, cmap = 'jet')

# st.pyplot(fig6.axes.figure)
