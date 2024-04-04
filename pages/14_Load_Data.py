import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import MinMaxScaler
from pygwalker.api.streamlit import init_streamlit_comm, get_streamlit_html
import streamlit.components.v1 as components

st.header("Get a quick overview of your dataset", divider="rainbow")

uploaded_files = st.file_uploader("Choose files", accept_multiple_files=True)
numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
# data=pd.DataFrame()

init_streamlit_comm()


# test
#
@st.cache_resource
def get_pyg_html(df: pd.DataFrame) -> str:
    # When you need to publish your application, you need set `debug=False`,prevent other users to write your config file.
    # If you want to use feature of saving chart config, set `debug=True`
    html = get_streamlit_html(
        df,
        spec="./chart_meta_0.json",
        themeKey="vega",
        use_kernel_calc=True,
        debug=False,
    )
    return html


def load_files(uploaded_files):
    # st.write(uploaded_file)
    list_df = []
    for uploaded_file in uploaded_files:
        if ".csv" in uploaded_file.name:
            list_df.append(pd.read_csv(uploaded_file))

        if ".fth" in uploaded_file.name:
            list_df.append(pd.read_feather(uploaded_file))
    return list_df


def filter_data(list_dfs):
    data = pd.concat(list_dfs)
    st.write("Raw Data", data.head(2))
    data.replace([np.inf, -np.inf, ""], np.nan, inplace=True)
    data.dropna(inplace=True, axis=1)
    data.dropna(inplace=True)
    # st.write('Data after removing Nan',data.head(2))
    if len(data) > 2:
        selector = VarianceThreshold()

        df_num = data.select_dtypes(include=numerics)
        num_col = [x for x in df_num.columns if "Meta" not in x and "meta" not in x]
        df_num = df_num[num_col]
        dat_sel = selector.fit_transform(df_num)
        col_sel = selector.get_feature_names_out()
        data1 = pd.DataFrame(dat_sel, columns=col_sel).reset_index(drop=True)
        # st.write('Data after removing Meta, Var=0',data1.head(2))
        cols_alpha = data.select_dtypes(exclude=numerics).columns
        data = pd.concat([data1, data[cols_alpha].reset_index(drop=True)], axis=1)
        data = data.dropna(axis=1)
        tab1, tab2, tab3 = st.tabs(["Dataframe", "Samples", "Summary"])
        tab1.write(data)
        tab2.write(data.sample(5))
        tab3.write(data.describe())
    return data


list_df = load_files(uploaded_files)

if len(list_df) > 0:
    data = filter_data(list_df)
    if "fieldId" in data.columns.tolist():
        data["fieldId"] = data["fieldId"].astype(str)
    if "tags" in data.columns.tolist():
        new = data["tags"].str.split(";", n=1, expand=True)
        for i in range(len(new.columns)):
            str_tag = "tag" + "_" + str(i)
            data[str_tag] = new[i]
    else:
        st.warning(
            "No column tags in your dataset or empty tags, taking Plate+Well as tags"
        )
        data["tags"] = data["Plate"] + "_" + data["Well"]
    # st.write("Tags", data.tags)
    col1, col2, col3 = st.columns(3)
    t = col1.radio("Time Series", ["yes", "no"], 1)

    col_sel = data.select_dtypes(include=numerics).columns.to_list()
    if t == "yes":
        st.warning("In dev....")
        t2 = col2.radio("Analyse", ["Contractility", "Calcium"], 0)
        dt = 0.031
        df_agg = data.copy()
        # df_agg=data[col_sel]
        if t2 == "Calcium":
            df_agg[col_sel] = (
                (-df_agg[col_sel])
                .add(df_agg[col_sel].min(axis=1), axis=0)
                .add(df_agg[col_sel].max(axis=1), axis=0)
            )
        # df_agg.drop('tags',axis=1,inplace=True)
        time_cols = sorted(
            [col for col in df_agg.columns if "time" in col],
            key=lambda x: int(x.split("_")[-1]),
        )
        sig_x = [i * dt for i in range(len(time_cols))]

        # df_agg['tags']=data['tags']
        #
        listofexp = df_agg["Plate"].unique()
        sel_col_exp = st.selectbox("Chose exps", listofexp)
        df_agg2 = df_agg[df_agg["Plate"] == sel_col_exp]
        df_agg2 = df_agg2.groupby(["tags"])
        for name, group in df_agg2:
            title = group.tags.values[0]
            cpd_names = group.Well.values
            # st.write(cpd_names)
            df_plt = group.set_index("Well")
            df_plt = df_plt[time_cols].T
            df_plt["t(s)"] = sig_x
            st.plotly_chart(
                px.line(
                    df_plt,
                    x="t(s)",
                    y=cpd_names,
                    width=800,
                    height=800,
                    title=title,
                    # line_shape='hv'
                ),
                theme="streamlit",
                use_container_width=True,
            )
            T = 1.0 / 30.0
            N = len(time_cols)
            import scipy
            from scipy import signal
            import matplotlib.pyplot as plt
            import pywt

            sel_col_fft = st.selectbox("Chose exp", df_plt.columns)
            yy = df_plt[sel_col_fft].values
            # yf = scipy.fftpack.fft(yy)
            widths = np.arange(1, N // 8)
            widths = np.geomspace(1, len(yy), num=10)
            time = np.linspace(0, 1, len(yy))
            sampling_period = np.diff(time).mean()
            # xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
            # cwtmatr = signal.cwt(yy, signal.ricker, widths)
            # cwtmatr_yflip = np.flipud(cwtmatr)
            cwtmatr, freqs = pywt.cwt(
                yy, widths, "mexh", sampling_period=sampling_period
            )
            fig_cwt, ax = plt.subplots()
            #     ax.imshow(cwtmatr, cmap='PRGn', aspect='auto',
            #    vmax=abs(cwtmatr).max(), vmin=-abs(cwtmatr).max())

            ax.pcolormesh(sig_x, freqs, cwtmatr)

            # ax.imshow(cwtmatr_yflip,  cmap='PRGn', aspect='auto',
            #     vmax=abs(cwtmatr).max(), vmin=-abs(cwtmatr).max())
            st.pyplot(fig_cwt, use_container_width=True)
    else:
        # if "tags" in data.columns.tolist():
        # data['tags']=data['Well']

        # st.write(df_plt)
        # N=len(time_cols)

        # df_fft=pd.DataFrame()
        # df_fft['xf']=xf
        # df_fft['yf']=2.0/N * np.abs(yf[:N//2])
        # components.html(get_pyg_html(df_fft), height=1000, scrolling=True)
        # # st.write('df_fft',df_fft)
        # fig_fft = px.scatter(df_fft,xf,yf)
        # st.plotly_chart(fig_fft, theme="streamlit", use_container_width=True)

        #

        # ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))

        # # plt.show()
        # st.pyplot(fig_fft, use_container_width=True)

        g = col3.radio("MinMax", ["yes", "no"])
        col_sel = data.select_dtypes(include=numerics).columns.to_list()
        if g == "no" and t == "no":
            components.html(get_pyg_html(data), height=1600, scrolling=True)
        if g == "yes" and t == "no":
            scaler = MinMaxScaler()
            col_sel = data.select_dtypes(include=numerics).columns.to_list()
            cols_alpha = data.select_dtypes(exclude=numerics).columns
            data_scaled = pd.DataFrame(
                scaler.fit_transform(data.select_dtypes(include=numerics)),
                columns=col_sel,
            )
            data_scaled["tags"] = data["tags"]
            data_scaled["fieldId"] = data["fieldId"]
            components.html(get_pyg_html(data_scaled), height=1000, scrolling=True)
            df_agg = data_scaled.groupby("tags").median().reset_index()
            # t= st.radio("Time Series", ["yes","no"],1)
            title = "MinMax profiles"
            # if t=="yes":
            #     st.warning("In dev....")
            #     df_agg=data[col_sel]
            #     # df_agg.drop('tags',axis=1,inplace=True)
            #     df_agg.columns = sorted([col for col in df_agg.columns if 'time' in col],key=lambda x: int(x.split("_")[-1]))
            #     df_agg['tags']=data['tags']
            #     col_sel = df_agg.select_dtypes(include=numerics).columns.to_list()
            #     title='Temporal Profiles'

            st.write(df_agg)
            cpd_names = df_agg.tags.values
            df_plt = df_agg.set_index("tags")
            df_plt = df_plt[col_sel].T
            fig4 = px.line(
                df_plt,
                x=col_sel,
                y=cpd_names,
                width=800,
                height=800,
                title=title,
            )
            st.plotly_chart(fig4, theme="streamlit", use_container_width=True)
            # components.html(get_pyg_html(df_plt), height=1000, scrolling=True)

            import umap

            #
            model = umap.UMAP(random_state=42, verbose=False).fit(data_scaled[col_sel])
            emb = model.transform(data_scaled[col_sel])

            df_all_umap = pd.DataFrame()
            df_all_umap["X_umap"] = emb[:, 0]
            df_all_umap["Y_umap"] = emb[:, 1]
            df_all_umap[cols_alpha] = data[cols_alpha]
            df_all_umap[col_sel] = data_scaled[col_sel]

            components.html(get_pyg_html(df_all_umap), height=1000, scrolling=False)
