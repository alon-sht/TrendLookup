# %%
import streamlit as st
import pandas as pd


#%%
def filter_dataframe(df, columns=None, allow_single_value_widgets=False):
    # Parse the df and get filter widgets based for provided columns
    count_unique = df.nunique()
    chosen_filters = st.multiselect(
        "Which columns to filter by?",
        options=count_unique[count_unique > 1].index.tolist(),
        key="chosen_filter",
    )
    columns = chosen_filters
    if not columns:  # if columns not provided, use all columns to create widgets
        return df
    if allow_single_value_widgets:
        threshold = 0
    else:
        threshold = 1
    widget_dict = {}
    filter_widgets = st.container()
    filter_widgets.warning(
        "After selecting filters press the 'Apply Filters' button at the bottom."
    )
    if not allow_single_value_widgets:
        filter_widgets.markdown(
            "Only showing columns that contain more than 1 unique value."
        )
    with filter_widgets.form(key="data_filters"):
        not_showing = []
        for y in df[columns]:
            if str(y) in st.session_state:  # update value from session state if exists
                selected_opts = st.session_state[str(y)]
            else:  # if doesnt exist use all values as defaults
                selected_opts = df[y].unique().tolist()
            if len(df[y].unique().tolist()) > threshold:  # checks if above threshold
                widget_dict[y] = st.multiselect(
                    label=str(y),
                    options=df[y].unique().tolist(),
                    default=selected_opts,
                    key=str(y),
                )
                df = df[df[y].isin(widget_dict[y])]
            else:  # if doesnt pass threshold
                not_showing.append(y)
        if not_showing:  # if the list is not empty, show this warning
            st.warning(
                f"Not showing filters for {' '.join(not_showing)} since they only contain one unique value."
            )
        submit_button = st.form_submit_button("Apply Filters")

    # reset_button = st.sidebar.button(
    #     "Reset All Filters",
    #     key="reset_buttons",
    #     on_click=reset_filter_widgets_to_default,
    #     args=(df, columns),
    # )
    filter_widgets.warning(
        "Dont forget to apply filters by pressing 'Apply Filters' at the bottom."
    )
    return df


def reset_filter_widgets_to_default(df, columns):
    for y in df[columns]:
        if str(y) in st.session_state:
            del st.session_state[y]


def get_session_state_dict():
    session_state_dict = {
        k: v
        for k, v in st.session_state.items()
        if "button" not in k and "file_uploader" not in k and "FormSubmitter" not in k
    }
    return session_state_dict


def update_session_state(
    update_all: bool = False, manual_session_state_dict: dict = None
):
    if update_all:
        session_state_dict = get_session_state_dict()
    if manual_session_state_dict:
        session_state_dict = manual_session_state_dict
    # st.session_state.update(session_state_dict)

    failed = []
    succeeded = []
    for k, v in session_state_dict.items():
        try:
            st.session_state[k] = v
            succeeded.append(k)
        except:
            failed.append(k)


def reset_state_to_default(custom_exceptions=None):
    default_exceptions = ["logged_in"]
    if type(custom_exceptions) == list:
        exceptions = default_exceptions + custom_exceptions
    else:
        exceptions = default_exceptions
    for y in st.session_state:
        if y not in exceptions:
            del st.session_state[y]
