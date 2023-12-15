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

# test
#
@st.cache_resource
def get_pyg_html(df: pd.DataFrame) -> str:
    # When you need to publish your application, you need set `debug=False`,prevent other users to write your config file.
    # If you want to use feature of saving chart config, set `debug=True`
    html = get_streamlit_html(df,spec="./chart_meta_0.json", themeKey='vega',use_kernel_calc=True, debug=False)
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
    data.replace([np.inf, -np.inf], np.nan, inplace=True)
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

        tab1, tab2, tab3 = st.tabs(["Dataframe", "Samples", "Summary"])
        tab1.write(data.head(2))
        tab2.write(data.sample(5))
        tab3.write(data.describe())
    return data


list_df = load_files(uploaded_files)

if len(list_df) > 0:
    data = filter_data(list_df)
    if "tags" in data.columns.tolist():
        new=data['tags'].str.split(";",n=1,expand=True)
        # st.write(len(new.columns))
        for i in range(len(new.columns)):
            str_tag='tag'+'_'+str(i)
            data[str_tag]=new[i]
    # pivot_df=pivot_df.apply(pd.to_numeric, errors='ignore')
    components.html(get_pyg_html(data), height=1000, scrolling=True)
    if "tags" in data.columns.tolist():
        g = st.radio("MinMax", ["yes", "no"])
        col_sel = data.select_dtypes(include=numerics).columns.to_list()
        if g == "no":
            desc = st.multiselect("Select some descriptors", col_sel)
            # st.write(desc)
            fig1 = px.box(data, x="tags", y=desc, notched=True)
            st.plotly_chart(fig1, theme="streamlit", use_container_width=False)

            col1, col2 = st.columns(2)
            desc1 = col1.multiselect("Select descriptors", col_sel)
            # desc2 = col2.selectbox("Select descriptor 2", df_num.columns)
            if len(desc1) < 2:
                st.write("Please select 2 descriptors")
            if len(desc1) == 2:
                fig2 = px.scatter(
                    data,
                    x=desc1[0],
                    y=desc1[1],
                    hover_data=["Plate", "Well"],
                    color="tags",
                    title="Raw Scatter",
                )
                st.plotly_chart(fig2, theme="streamlit", use_container_width=False)

        scaler = MinMaxScaler()
        data_scaled = pd.DataFrame(
            scaler.fit_transform(data.select_dtypes(include=numerics)), columns=col_sel
        )

        data_scaled["tags"] = data["tags"]
        df_agg = data_scaled.groupby("tags").median().reset_index()
        st.write(df_agg)
        cpd_names = df_agg.tags.values
        df_plt = df_agg.set_index("tags")
        df_plt = df_plt[col_sel].T
        fig3 = px.line(
            df_plt,
            x=col_sel,
            y=cpd_names,
            width=800,
            height=800,
            title="MinMax profiles",
        )
        st.plotly_chart(fig3, theme="streamlit", use_container_width=False)

        if g == "yes":
            scaler = MinMaxScaler()
            col_sel = data.select_dtypes(include=numerics).columns.to_list()
            data_scaled = pd.DataFrame(
                scaler.fit_transform(data.select_dtypes(include=numerics)),
                columns=col_sel,
            )

            data_scaled["tags"] = data["tags"]
            desc = st.multiselect("Select some descriptors", col_sel)
            fig1 = px.box(data_scaled, x="tags", y=desc, notched=True)
            st.plotly_chart(fig1, theme="streamlit", use_container_width=False)

            col1, col2 = st.columns(2)
            desc1 = col1.multiselect("Select descriptors", col_sel)
            # desc2 = col2.selectbox("Select descriptor 2", df_num.columns)
            if len(desc1) < 2:
                st.write("Please select 2 descriptors")
            if len(desc1) == 2:
                fig2 = px.scatter(
                    data_scaled,
                    x=desc1[0],
                    y=desc1[1],
                    color="tags",
                    title="MinMax Scatter",
                )
                st.plotly_chart(fig2, theme="streamlit", use_container_width=False)
            import umap

            #
            model = umap.UMAP(random_state=42, verbose=False).fit(data_scaled[col_sel])
            emb = model.transform(data_scaled[col_sel])

            df_all_umap = pd.DataFrame()
            df_all_umap["X_umap"] = emb[:, 0]
            df_all_umap["Y_umap"] = emb[:, 1]
            df_all_umap["tags"] = data_scaled["tags"]

            fig3 = px.scatter(
                df_all_umap,
                x="X_umap",
                y="Y_umap",
                title="umap",
                hover_data=["tags"],
                color="tags",
            )
            st.plotly_chart(fig3, theme="streamlit", use_container_width=False)  #
    else:
        st.warning("No column tags in your dataset")
        g = st.radio("MinMax", ["yes", "no"])
        if g == "no":
            col_sel = data.select_dtypes(include=numerics).columns.to_list()
            desc = st.selectbox("Select descriptor", col_sel)
            fig1 = px.box(data, y=desc, notched=True)
            st.plotly_chart(fig1, theme="streamlit", use_container_width=False)

            col1, col2 = st.columns(2)
            desc1 = col1.multiselect("Select descriptors", col_sel)
            # desc2 = col2.selectbox("Select descriptor 2", df_num.columns)
            if len(desc1) < 2:
                st.write("Please select 2 descriptors")
            if len(desc1) == 2:
                fig2 = px.scatter(
                    data,
                    x=desc1[0],
                    y=desc1[1],
                    title="Raw Scatter",
                )
                st.plotly_chart(fig2, theme="streamlit", use_container_width=False)

        if g == "yes":
            col_sel = data.select_dtypes(include=numerics).columns.to_list()
            scaler = MinMaxScaler()
            data_scaled = pd.DataFrame(
                scaler.fit_transform(data.select_dtypes(include=numerics)),
                columns=col_sel,
            )

            desc = st.selectbox("Select descriptor", col_sel)
            fig1 = px.box(data_scaled, y=desc, notched=True)
            st.plotly_chart(fig1, theme="streamlit", use_container_width=False)

            col1, col2 = st.columns(2)
            desc1 = col1.multiselect("Select descriptors", col_sel)
            # desc2 = col2.selectbox("Select descriptor 2", df_num.columns)
            if len(desc1) < 2:
                st.write("Please select 2 descriptors")
            if len(desc1) == 2:
                fig2 = px.scatter(
                    data_scaled,
                    x=desc1[0],
                    y=desc1[1],
                    title="MinMax Scatter",
                )
                st.plotly_chart(fig2, theme="streamlit", use_container_width=False)
            import umap

            model = umap.UMAP(random_state=42, verbose=False).fit(data_scaled[col_sel])
            emb = model.transform(data_scaled[col_sel])
            cols_alpha = data.select_dtypes(exclude=numerics).columns
            df_all_umap = pd.DataFrame()
            df_all_umap = data[cols_alpha]
            df_all_umap["X_umap"] = emb[:, 0]
            df_all_umap["Y_umap"] = emb[:, 1]

            fig3 = px.scatter(
                df_all_umap, x="X_umap", y="Y_umap", title="umap", hover_data=cols_alpha
            )
            st.plotly_chart(fig3, theme="streamlit", use_container_width=False)
