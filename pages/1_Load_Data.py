import pandas as pd
import plotly.express as px
import streamlit as st

uploaded_files = st.file_uploader("Choose files", accept_multiple_files=True)
# st.write(uploaded_file)
list_df = []
for uploaded_file in uploaded_files:
    if ".csv" in uploaded_file.name:
        list_df.append(pd.read_csv(uploaded_file))

    if ".fth" in uploaded_file.name:
        list_df.append(pd.read_feather(uploaded_file))


if len(list_df) > 0:
    import numpy as np
    from sklearn.feature_selection import VarianceThreshold

    data = pd.concat(list_df)
    st.write('Raw Data', data)
    data.replace([np.inf, -np.inf], np.nan, inplace=True)
    data.dropna(inplace=True,axis=1)
    data.dropna(inplace=True)
    st.write('Data after removing Nan',data)
    if len(data) > 2:
        selector = VarianceThreshold()
        numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
        df_num = data.select_dtypes(include=numerics)
        dat_sel = selector.fit_transform(df_num)
        col_sel = selector.get_feature_names_out()
        data1 = pd.DataFrame(dat_sel, columns=col_sel).reset_index(drop=True)
        cols_alpha = data.select_dtypes(exclude=numerics).columns
        data = pd.concat([data1, data[cols_alpha].reset_index(drop=True)], axis=1)

        tab1, tab2, tab3 = st.tabs(["Dataframe", "Samples", "Summary"])
        tab1.write(data)
        tab2.write(data.sample(5))
        tab3.write(data.describe())

        from sklearn.preprocessing import MinMaxScaler

        g = st.radio("MinMax", ["yes", "no"])
        if g == "no":
            desc = st.selectbox("Select descriptor", col_sel)
            fig1 = px.box(data, x="tags", y=desc, notched=True)
            st.plotly_chart(fig1, theme="streamlit", use_container_width=True)

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
                st.plotly_chart(fig2, theme="streamlit", use_container_width=True)

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
            df_plt, x=col_sel, y=cpd_names, width=1400, height=1000, title="MinMax profiles"
        )
        st.plotly_chart(fig3, theme="streamlit", use_container_width=True)

        if g == "yes":
            scaler = MinMaxScaler()
            data_scaled = pd.DataFrame(
                scaler.fit_transform(data.select_dtypes(include=numerics)), columns=col_sel
            )
            data_scaled["tags"] = data["tags"]
            desc = st.selectbox("Select descriptor", col_sel)
            fig1 = px.box(data_scaled, x="tags", y=desc, notched=True)
            st.plotly_chart(fig1, theme="streamlit", use_container_width=True)

            col1, col2 = st.columns(2)
            desc1 = col1.multiselect("Select descriptors", col_sel)
            # desc2 = col2.selectbox("Select descriptor 2", df_num.columns)
            if len(desc1) < 2:
                st.write("Please select at 2 descriptors")
            if len(desc1) == 2:
                fig2 = px.scatter(
                    data_scaled,
                    x=desc1[0],
                    y=desc1[1],
                    color="tags",
                    title="MinMax Scatter",
                )
                st.plotly_chart(fig2, theme="streamlit", use_container_width=True)
            import umap

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
            st.plotly_chart(fig3, theme="streamlit", use_container_width=True)
