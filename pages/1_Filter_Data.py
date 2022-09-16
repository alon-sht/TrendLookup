# %%
import pandas as pd
import streamlit as st
from PIL import Image
from streamlit_extras.switch_page_button import switch_page
from src.functions import filter_dataframe, update_session_state
from src.data_functions import sort_samples

update_session_state(update_all=True)
# %%

st.set_page_config(
    layout="wide", page_title="Filter Data", page_icon=Image.open("fav.ico")
)
hide_streamlit_style = """
              <style>
              #MainMenu {visibility: hidden;}
              footer {visibility: hidden;}

            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True)


# st.write(st.session_state)
def main():
    st.session_state["filtered_df"] = filter_dataframe(st.session_state["df"], None)
    st.session_state["sorter_mean"], st.session_state["sorter_median"] = sort_samples(
        st.session_state["filtered_df"]
    )
    st.sidebar.metric("DataFrame Length", st.session_state.df.shape[0])
    st.sidebar.metric("Filtered DataFrame", st.session_state.filtered_df.shape[0])


# if 'filtered_df' in st.session_state:
#     st.session_state['filtered_df']=sort_samples(st.session_state['filtered_df'])

if __name__ == "__main__":
    if "df" in st.session_state:
        main()
    else:
        import time

        seconds_to_wait = 4
        for i in range(seconds_to_wait):
            st.error(
                f"# Data not uploaded. Returning to homepage in {seconds_to_wait-i} seconds"
            )
            time.sleep(1)
        switch_page("TrendsApp")
