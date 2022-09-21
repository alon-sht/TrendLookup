# %%
from unicodedata import category
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
    st.session_state["filtered_metadata"] = filter_dataframe(st.session_state["metadata"], None)
    st.session_state["sorter_mean"], st.session_state["sorter_median"] = sort_samples(
        st.session_state["data"]
    )
    st.session_state['filtered_df']=pd.merge(st.session_state['filtered_metadata'],st.session_state['data'],on='SampleID',how='inner')
    st.session_state['meta_columns']=st.session_state['filtered_metadata'].columns.tolist()
    for col in st.session_state['meta_columns']:
        st.session_state['filtered_df'][col]=st.session_state['filtered_df'][col].astype('category')
    st.sidebar.metric("Total No of Samples ", st.session_state.metadata.shape[0])
    st.sidebar.metric("Filtered Samples", st.session_state.filtered_metadata.shape[0])
    def memory_usage(df):
        return(round(df.memory_usage(deep=True).sum() / 1024 ** 2, 2))
    
    with st.sidebar.expander("Memory Usage Monitor"):
        st.write(f"filtered_metadata : {memory_usage(st.session_state['filtered_metadata'])}")
        st.write(f"metadata : {memory_usage(st.session_state['metadata'])}")
        st.write(f"data : {memory_usage(st.session_state['data'])}")
        st.write(f"df : {memory_usage(st.session_state['filtered_df'])}")

# if 'filtered_df' in st.session_state:
#     st.session_state['filtered_df']=sort_samples(st.session_state['filtered_df'])

if __name__ == "__main__":
    if 'authentication_status' in st.session_state:
        if st.session_state["authentication_status"]:
            if "data" in st.session_state:
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
    else:
        switch_page("TrendsApp")
    
