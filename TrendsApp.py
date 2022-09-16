# %%
import os
import pandas as pd
import plotly.express as px
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu, kruskal, wilcoxon
import streamlit as st
from PIL import Image
from src.check_password import check_password
from io import BytesIO
from itertools import combinations, combinations
from streamlit_extras.switch_page_button import switch_page
from src.functions import reset_state_to_default

# %%

st.set_page_config(
    layout="wide", page_title="TrendAnalysis", page_icon=Image.open("fav.ico")
)
hide_streamlit_style = """
              <style>
              #MainMenu {visibility: hidden;}
              footer {visibility: hidden;}

            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

# %%


def st_main_header():
    global sample_data_message, message_box
    st.title("Trend Analysis")
    st.sidebar.image("Mybiotics_LOGO - Large.png", width=350)
    st.sidebar.markdown("""---""")


def st_upload_file():
    upload_column = st.container()
    upload_column.subheader("File Upload")
    upload_column.text("Upload file to start working")
    upload_data_widget = upload_column.file_uploader(label="Upload File", type=["csv"])
    # message=None
    if upload_data_widget:
        open_file(upload_data_widget.getvalue())
        st.success("Custom Data Loaded")
    sample_data = upload_column.checkbox("Use Sample Data")
    if sample_data:
        st.session_state["df"] = pd.read_csv("sample_data.csv")
        st.success("Using test data")
    upload_column.markdown("""---""")

    if "df" in st.session_state:
        # if st.button("Proceed"):
        st.session_state["df"] = process_data(st.session_state["df"])
        switch_page("Filter Data")


@st.experimental_memo(show_spinner=False,)  # allow_output_mutation=True)
def open_file(file_uploader_data):
    st.session_state["df"] = pd.read_csv(BytesIO(file_uploader_data), engine="python")


def process_data(df):

    # Metadata columns astype str
    meta_columns = [
        x
        for x in st.session_state["df"].columns.tolist()
        if x not in ["OTU", "level_1", "RA", "Notes"]
    ]

    st.session_state["df"][meta_columns] = st.session_state["df"][meta_columns].astype(
        str
    )
    st.session_state["meta_columns"] = meta_columns

    return df


def main():

    if "df" in st.session_state:
        st.button("Reset and upload new file", on_click=reset_state_to_default)

    if "df" not in st.session_state:
        st_upload_file()


# st.write(list(st.session_state.keys()))

if __name__ == "__main__":
    main()
