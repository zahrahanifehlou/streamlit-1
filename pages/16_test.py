import streamlit as st
from streamlit_plotly_events import plotly_events
import numpy as np
import plotly.express as px

x = np.random.rand(100)
y = np.random.rand(100)
# Writes a component similar to st.write()
fig = px.scatter(x=x, y=y)
selected_points = plotly_events(fig,click_event=True)
st.write(selected_points)