# %%
import pandas as pd
import streamlit as st

# %%


def sort_samples(df):
    df['RA']=df['RA'].astype(float)
    sorter_mean = (
        df.groupby("OTU")
        .agg({"RA": "mean"})
        .sort_values("RA", ascending=False)
        .index.tolist()
    )
    sorter_median = (
        df.groupby("OTU")
        .agg({"RA": "mean"})
        .sort_values("RA", ascending=False)
        .index.tolist()
    )
    # sorter_choose = st.radio(
    #     "Sort samples by Mean or Median",
    #     options=["Mean", "Median"],
    #     index=0,
    #     key="mean_or_median",
    # )  # ,horizontal=True)
    # if sorter_choose == "Mean":
    #     sorter = sorter_mean
    # elif sorter_choose == "Median":
    #     sorter = sorter_median
    # sorterIndex = dict(zip(sorter, range(len(sorter))))
    # df["ind"] = df["OTU"].map(sorterIndex)
    # df.sort_values(by="ind", inplace=True)
    return sorter_mean, sorter_median
    # return df,sorter
