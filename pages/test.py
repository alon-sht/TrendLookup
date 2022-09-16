#%%.
import streamlit as st
import pandas as pd
import plotly.express as px

#%%
df = pd.read_csv(
    r"H:\My Drive\MyBiotics\Experiments and Results\RD628_zr4871.3.16S_211025.zymo\zr4871.3.16S_211025.zymo\00...AllSamples.Bac16Sv34\Composition_Barplots\Composition_Summary_L3.txt",
    delimiter="\t",
    skiprows=1,
)
df = df.set_index("#OTU ID")
df = df.drop(
    [
        "RD645.1",
        "RD645.4",
        "RD645.5",
        "RD645.6",
        "RD645.9",
        "RD645.10",
        "RD645.11",
        "RD645.14",
        "RD645.15",
        "RD645.16",
        "RD645.19",
        "RD645.20",
        "RD624.2",
        "RD624.5",
        "RD624.6",
        "RD624.9",
        "RD624.10",
        "RD624.13",
        "RD624.14",
        "RD646.1",
        "RD646.2",
        "RD646.3",
        "RD646.4",
        "RD646.5",
        "RD646.6",
        "RD646.7",
        "RD646.8",
        "RD646.9",
        "RD646.10",
        "RD646.11",
        "RD646.12",
        "RD646.13",
        "RD646.14",
        "RD646.15",
        "RD646.16",
    ],
    axis="columns",
)
df_sub = df.subtract(df["SDD005"], axis="index")
dfm = df.melt(ignore_index=False)
dfm_sub = df_sub.melt(ignore_index=False)
# %%
# %%
fig = px.density_heatmap(
    dfm, y=dfm.index, x="variable", z="value", width=1200, height=1000
)
st.plotly_chart(fig, use_container_width=True)

fig = px.density_heatmap(
    dfm_sub,
    y=dfm.index,
    x="variable",
    z="value",
    width=1200,
    height=1000,
    color_continuous_scale=px.colors.sequential.RdBu,
    color_continuous_midpoint=0,
)
st.plotly_chart(fig, use_container_width=True)
# %%
