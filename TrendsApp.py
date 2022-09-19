# %%
import os
import pandas as pd
import streamlit as st
from PIL import Image
# from src.check_password import check_password
from io import BytesIO
from itertools import combinations, combinations
from streamlit_extras.switch_page_button import switch_page
from src.functions import reset_state_to_default
import streamlit_authenticator as stauth


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
    # upload_data_widget = upload_column.file_uploader(label="Upload File", type=["csv"])
    sample_data = upload_column.checkbox("Use Sample Data",key='sample_data')
    
    upload_metadata = upload_column.file_uploader(label="Upload Metadata", type=["csv"])
    if upload_metadata:
        st.session_state['metadata']=open_file(upload_metadata.getvalue())
    upload_data = upload_column.file_uploader(label="Upload Data", type=["csv"])
    if upload_data:
        st.session_state['data']=open_file(upload_data.getvalue())
            
            # st.session_state['df']=pd.merge(st.session_state['metadata'],st.session_state['data'],on='SampleID')

    # if upload_data_widget:
    #     st.session_state['df']=open_file(upload_data_widget.getvalue())
    #     st.success("Custom Data Loaded")
    
    if sample_data:
        temp=pd.read_csv("sample_data.csv")
        st.session_state["metadata"] = temp[['SampleID','Notes','ExperimentName','Day/TimePoint','ReplicateGroup','DonorName']].drop_duplicates()
        st.session_state["data"] = temp[['OTU','RA','SampleID']]
        del temp
    upload_column.markdown("""---""")

    if "data" and 'metadata' in st.session_state:
        if st.button("Proceed"):
        # st.session_state["df"] = process_data(st.session_state["df"])
            switch_page("Filter Data")


@st.experimental_memo(show_spinner=False,)  # allow_output_mutation=True)
def open_file(file_uploader_data):
    df = pd.read_csv(BytesIO(file_uploader_data), engine="python",dtype='object')
    df.replace(np.nan,'None')
    return df

# def process_data(df):

#     # Metadata columns astype str
#     meta_columns = [
#         x
#         for x in st.session_state["df"].columns.tolist()
#         if x not in ["OTU", "level_1", "RA", "Notes"]
#     ]

#     # st.session_state["df"][meta_columns] = st.session_state["df"][meta_columns].astype(
#     #     str
#     # )
#     st.session_state["meta_columns"] = meta_columns

#     return df


def main():

    if "df" in st.session_state:
        st.button("Reset and upload new file", on_click=reset_state_to_default)

    if "df" not in st.session_state:
        st_upload_file()


# st.write(list(st.session_state.keys()))

if __name__ == "__main__":
    import yaml
    with open('secrets/config.yaml') as file:
        config = yaml.load(file, Loader=stauth.SafeLoader)

    authenticator = stauth.Authenticate(
        config['credentials'],
        config['cookie']['name'],
        config['cookie']['key'],
        config['cookie']['expiry_days'],
        config['preauthorized']
    )
    name, authentication_status, username = authenticator.login('Login', 'main')

    # if authentication_status:
    #     authenticator.logout('Logout', 'main')
    #     st.write(f'Welcome *{name}*')
    #     # st.title('Some content')
    # elif authentication_status == False:
    #     st.error('Username/password is incorrect')
    # elif authentication_status == None:
    #     st.warning('Please enter your username and password')
        
    if st.session_state["authentication_status"]:
        st.write(f'Welcome *{st.session_state["name"]}*')
        authenticator.logout('Logout', 'main')
        # if st.button("Change Password"):
        # with st.expander("Change Password"):
        #     try:
        #         if authenticator.reset_password(username, 'Reset password'):
        #             st.success('Password modified successfully')
        #     except Exception as e:
        #         st.error(e)
        main()
        
    elif st.session_state["authentication_status"] == False:
        st.error('Username/password is incorrect')
    elif st.session_state["authentication_status"] == None:
        st.warning('Please enter your username and password')
    
