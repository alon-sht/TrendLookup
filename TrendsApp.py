# %%
import os
import pandas as pd
import plotly.express as px
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu,kruskal,wilcoxon
import streamlit as st
from PIL import Image
from src.check_password import check_password
from io import BytesIO
from itertools import combinations,combinations

st.set_page_config(
    layout="wide", page_title="TrendAnalysis", page_icon=Image.open("fav.ico")
)
hide_streamlit_style = """
              <style>
              #MainMenu {visibility: hidden;}
              footer {visibility: hidden;}
              [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {
                     width: 500px;
                     }
              [data-testid="stSidebar"][aria-expanded="false"] > div:first-child {
                     width: 500px;
                     margin-left: -500px;
                     }
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True)


def st_main_header():
    global sample_data_message, message_box
    st.title("Trend Analysis")
    st.sidebar.image("Mybiotics_LOGO - Large.png", width=350)
    st.sidebar.markdown("""---""")
    message_box = st.container()
    sample_data_message = message_box.error(
        "#### Upload your data using the 'File Upload' widget OR use sample data by clicking 'Use Sample Data' checkbox."
    )


def st_sidebar_pick_file():
    # Sidebar dropdown to pick file to look at
    global level, df

    file_select = st.sidebar.container()
    file_options = {
        "Level 3 (Class) - Stool Only": "all_donors L3.csv",
        "Level 4 (Order) - Stool Only": "all_donors L4.csv",
        "Level 5 (Family) - Stool Only": "all_donors L5.csv",
        "Level 6 (Genus) - Stool Only": "all_donors L6.csv",
        # 'Level 7 (Species) - Stool Only':"all_donors L7.csv",
        "Level 3 (Class) - All Samples": "all_samples L3.csv",
        "Level 4 (Order) - All Samples": "all_samples L4.csv",
        "Level 5 (Family) - All Samples": "all_samples L5.csv",
    }

    # level=file_select.selectbox(label='Pick file to use', options=[None]+list(file_options.keys()), index=0)#'6 - All Samples'
    file_select.markdown("""---""")
    wd = ""
    df = pd.DataFrame(columns=["OTU", "RA"])
    # if level is not None:
    try:
        file = file_options[level]
        df = pd.read_csv(os.path.join(wd, file), engine="python")
    except:
        st.warning("This file doesn't currently work")
        df = pd.DataFrame()


def st_sidebar_upload_file():
    global upload_data_widget, df, sample_data
    upload_column = st.sidebar.container()
    upload_column.subheader("File Upload")
    upload_column.text("Upload file to start working")
    upload_data_widget = upload_column.file_uploader(label="Upload File", type=["csv"])
    if upload_data_widget:
        df = pd.read_csv(BytesIO(upload_data_widget.getvalue()), engine="python")
        sample_data_message.success("Custom Data Loaded")
    sample_data = upload_column.checkbox("Use Sample Data")
    if sample_data:
        df = pd.read_csv("sample_data.csv")
    upload_column.markdown("""---""")


def st_main_raw_data():
    # Show raw data table in main container

    raw_data = st.container()

    raw_data.subheader("Raw Data")
    show_raw = raw_data.checkbox("Show Raw Data", value=False)
    if show_raw:
        raw_data.write(df.astype(str), use_container_width=True)
    # raw_data.markdown("""---""")
    # data=st.container()


def st_sidebar_sort_samples():
    global meta_columns, sorter, data
    meta_columns = [
        x for x in df.columns.tolist() if x not in ["OTU", "level_1", "RA", "Notes"]
    ]
    df[meta_columns] = df[meta_columns].astype(str)
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
    sorter_choose = st.radio(
        "Sort samples by Mean or Median",
        options=["Mean", "Median"],
        index=0,
        key="mean_or_median",
    )  # ,horizontal=True)
    st.markdown("""---""")
    if sorter_choose == "Mean":
        sorter = sorter_mean
    elif sorter_choose == "Median":
        sorter = sorter_median
    sorterIndex = dict(zip(sorter, range(len(sorter))))
    df["ind"] = df["OTU"].map(sorterIndex)
    df.sort_values(by="ind", inplace=True)


def show_ra_of_all():
    ra_all = st.container()
    ra_all.subheader("Relative Abundance of All Samples")
    ra_all.text("For many samples (>30) this might be very slow")
    ra_all.text(
        f"You currently have {len(df_filtered['SampleID'].unique())} samples after applying filters with {len(df_filtered['OTU'].unique())} unique OTUs"
    )
    load_plot = ra_all.checkbox("Load Plot", value=False)
    ra_all_plot = ra_all.container()
    if load_plot:
        ra_all_legend = ra_all_plot.checkbox("Show Legend on Plot", value=False)
        fig = px.bar(
            df,
            x="SampleID",
            y="RA",
            color="OTU",
            barmode="stack",
        )
        fig.update_xaxes(
            dtick=1,
        )
        # fig.show(renderer="svg")
        fig.update_layout(showlegend=ra_all_legend)
        ra_all_plot.plotly_chart(fig, use_container_width=True)
        # li=[]
        fig.data[0]["showlegend"] == False
        # for x in range(10):
        # if fig.data[x]['showlegend']==True:
        # li.append(fig.data[x]['name'])
        # ra_all_plot.text(li)

    ra_all.markdown("""---""")


def st_sidebar_data_filters():
    # Show data filters in the sidebar inside an expander
    # Filters are updated with the press of a button (and not automatically)
    # Also number of samples before and after filtering is shown

    global df_filtered

    filters = st.sidebar.container()
    filters.subheader("Data Filters")
    filter_widgets = filters.expander(
        "Filter Widgets (Click to Expand). After selecting filters click the button at the bottom."
    )
    query = f""
    widget_dict = {}
    form = filter_widgets.form("form1")
    for col in meta_columns:
        if col not in ["SampleID", "Replicates", "ReplicateGroup", "ZymoID"]:
            widget_dict[col] = form.multiselect(
                label=col,
                options=df[col].unique().tolist(),
                default=df[col].unique().tolist(),
                key=str(col),
            )
            query += f"`{col}`  in {widget_dict[col]} & "

    submit = form.form_submit_button("Filter Data")
    df_filtered = df.query(query[:-2])
    filters.markdown(
        f"##### No. of samples before filtering: {len(df['SampleID'].unique())}"
    )
    filters.markdown(
        f"##### No. of samples after filtering: {len(df_filtered['SampleID'].unique())} ({len(df_filtered['SampleID'].unique()) * 100 / len(df['SampleID'].unique())}% of all samples)"
    )
    filters.markdown("""---""")


def st_sidebar_top_bacteria_slider():
    # Slider to select number of top bacteria to show
    global top_val, df_top, font_size
    top_val = st.sidebar.slider(
        label="Select number of bacteria to show",
        min_value=1,
        max_value=len(df_filtered["OTU"].unique()),
        value=10,
        key="top_val",
    )
    df_top = df_filtered[df_filtered["OTU"].isin(sorter[:top_val])]
    global df_top_corr, df1

    df1 = df_filtered.pivot(index=meta_columns, columns=["OTU"], values="RA")
    df_top_corr = df_top.pivot(index=meta_columns, columns=["OTU"], values="RA").corr(
        method="spearman"
    )
    df_top_corr = df_top_corr.mask(np.tril(np.ones(df_top_corr.shape), -1).astype(bool))
    font_size = st.sidebar.slider(
        label="Font Size", min_value=1, max_value=25, value=12, key="font_size"
    )


def st_main_top_bacteria_plot():
    # A box plot showing the relative abundance of the top bacteria
    top_bacteria_boxplot = st.container()
    top_bacteria_boxplot.subheader(f"Top {top_val} bacteria (ordered by mean)")
    top_bacteria_boxplot.text(
        "Change number of top bacteria by using the slider in the sidebar"
    )
    fig_all = px.box(
        df_top,
        y="OTU",
        x="RA",
        height=900,
    )  # template='plotly_white')
    fig_all.update_layout(
        font=dict(
            size=font_size,
        )
    )
    fig_all.update_traces(
        boxmean=True,
        orientation="h",
    )
    fig_all.update_yaxes(
        autorange="reversed",
        dtick=1,
    )
    fig_all.update_xaxes(title="Relative Abundance")
    top_bacteria_boxplot.plotly_chart(fig_all, use_container_width=True)
    top_bacteria_boxplot.markdown("""---""")


def st_main_correlation_heatmap_between_top_bac():
    # Show heatmap of top correlations

    top_bac_correlation_heatmap = st.container()
    corr_plot = px.imshow(
        df_top_corr,
        text_auto=".2f",
        template="plotly_white",
        color_continuous_scale="RdBu",
    )
    corr_plot.update_layout(autosize=True, height=900)
    corr_plot.update_xaxes(showticklabels=False, showgrid=False)
    corr_plot.update_yaxes(showticklabels=False, showgrid=False)
    corr_plot.update_layout(
        font=dict(
            size=font_size,
        )
    )
    top_bac_correlation_heatmap.subheader(
        f"Correlation Matrix of the top {top_val} bacteria"
    )
    top_bac_correlation_heatmap.text("Hover over the plot to see OTU names")
    top_bac_correlation_heatmap.text(
        "Change number of shown bacteria by using the slider in the sidebar"
    )
    top_bac_correlation_heatmap.plotly_chart(corr_plot, use_container_width=True)
    top_bac_correlation_heatmap.markdown("""---""")


def st_main_top_correlations_plot():
    # Show bar chart of main correlations
    top_correlation_plot = st.container()
    df_top_corr.index.name = "OTU1"
    df_top_corr_new = (
        df_top_corr.reset_index()
        .melt(id_vars="OTU1", value_vars=df_top_corr.columns)
        .sort_values("value", ascending=False)
        .dropna()
    )
    df_top_corr_new["value"] = df_top_corr_new["value"].round(3)
    df_top_corr_new = df_top_corr_new[df_top_corr_new["value"] < 1]
    df_top_corr_new["x"] = df_top_corr_new["OTU1"].astype(str) + df_top_corr_new[
        "OTU"
    ].astype(str)
    df_top_corr_new["abs"] = abs(df_top_corr_new["value"])
    df_top_corr_new["Correlation"] = np.where(
        df_top_corr_new["value"] < 0, "Negative", "Positive"
    )
    top_correlation_plot.subheader(f"Top Correlations")
    top_correlation_plot.text("Hover over the plot to see OTU names")
    top_correlation_plot.text(
        "Change number of correlations shown by using the slider in the sidebar"
    )
    top_corr = top_correlation_plot.slider(
        label="Number of top correlations to show",
        min_value=1,
        max_value=len(df_filtered["OTU"].unique()),
        value=50,
    )
    df_top_corr_new = (
        df_top_corr_new.sort_values("abs", ascending=False)
        .head(top_corr)
        .sort_values("value", ascending=False)
    )
    text_on_plot = top_correlation_plot.checkbox("Show text on plot", value=False)
    top_corr_plot = px.bar(
        df_top_corr_new,
        x="x",
        y="value",
        hover_data=["OTU", "OTU1"],
        color="Correlation",
        text_auto=text_on_plot,
    )
    top_corr_plot.update_traces(
        textfont_size=12,
        textangle=0,
        textposition="outside",
    )  # cliponaxis=False)
    top_corr_plot.update_xaxes(showticklabels=False, title="OTU Pair")
    top_corr_plot.update_yaxes(title="Correlation")
    top_corr_plot.update_layout(
        font=dict(
            size=font_size,
        )
    )
    top_correlation_plot.plotly_chart(top_corr_plot, use_container_width=True)
    top_correlation_plot.markdown("""---""")


def st_main_bacteria_to_compare_selection():
    # Selection to bacteria to compare
    global x, y
    compare_bac_selection = st.container()
    compare_bac_selection.subheader("Select Bacteria to Compare")
    compare_bac_selection.text("Dropdowns are sorted by mean relative abundance")
    x = compare_bac_selection.selectbox(
        label="Bacteria 1", options=sorter[:], key="bac1"
    )
    y = compare_bac_selection.selectbox(
        label="Bacteria 2", options=sorter[:], index=1, key="bac2"
    )
    compare_bac_selection.markdown("""---""")
    global df2_piv
    df2_piv = df1.reset_index()
    df2_piv["ratio"] = df2_piv[x] / df2_piv[y]
    df2_piv["ratio"]=df2_piv["ratio"].replace(np.inf, np.nan)
    
    
    
    
    

def print_corr(df, x, y, param, where):
    # Function to print correlations
    corr = spearmanr(df[[x, y]])
    if corr[1] < 0.05:
        return where.success(
            f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]}"
        )
    elif corr[1] > 0.05:
        return where.error(
            f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]}"
        )
    else:
        return where.info(
            f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]}"
        )


def st_main_ra_plot_of_selected_bacteria():
    # stacked bar chart of RA of the two selected bacteria

    ra_of_selected_bacteria = st.container()

    ra_of_selected_bacteria.subheader(
        f"Relative abundance of {x.split(';')[-1]} and {y.split(';')[-1]}"
    )
    split_by_donor = ra_of_selected_bacteria.checkbox(
        "Split plot by Donor", key="split_by_donor"
    )

    if split_by_donor:
        split = "DonorName"
    else:
        split = None

    fig1 = px.bar(df2_piv, x="SampleID", y=[x, y], facet_col=split)
    fig1.update_xaxes(matches=None)
    fig1.update_layout(showlegend=False)
    fig1.update_layout(
        font=dict(
            size=font_size,
        )
    )
    ra_of_selected_bacteria.plotly_chart(fig1, use_container_width=True)
    ra_of_selected_bacteria.markdown("""---""")


def st_main_correlation_scatter_between_selected_baceria():
    # Scatter plot of RA of the two selected bacteria
    correlation_between_selected_bacteria = st.container()
    correlation_between_selected_bacteria.subheader(
        "Spearman correlations between two selected bacteria"
    )
    correlation_between_selected_bacteria.text(
        f"Correlation between {x.split(';')[-1]} and {y.split(';')[-1]} - overall and for each color group"
    )
    color = correlation_between_selected_bacteria.selectbox(
        label="Color By",
        options=[None] + meta_columns,
        index=2,
        key="color_by_corr_bacteria_plot",
    )

    col1, col2 = correlation_between_selected_bacteria.columns(2)
    corr_expander = correlation_between_selected_bacteria.expander(
        label="Spearman Correlation (Click to Expand)", expanded=False
    )
    fig2_trend = px.scatter(
        df2_piv,
        x=x,
        y=y,
        color=color,
        hover_data=meta_columns,
        trendline="ols",
        trendline_scope="overall",
    )
    fig2_trend.layout.xaxis.title = fig2_trend.layout.xaxis.title["text"].split(";")[-1]
    fig2_trend.layout.yaxis.title = fig2_trend.layout.yaxis.title["text"].split(";")[-1]
    fig2_trend.update_layout(
        font=dict(
            size=font_size,
        )
    )
    col1.plotly_chart(fig2_trend, use_container_width=True)
    print_corr(df2_piv, x, y, "Overall", corr_expander)
    fig2_trend_each = px.scatter(
        df2_piv,
        x=x,
        y=y,
        color=color,
        hover_data=meta_columns,
        trendline="ols",
    )
    fig2_trend_each.layout.xaxis.title = fig2_trend_each.layout.xaxis.title[
        "text"
    ].split(";")[-1]
    fig2_trend_each.layout.yaxis.title = fig2_trend_each.layout.yaxis.title[
        "text"
    ].split(";")[-1]
    fig2_trend_each.update_layout(
        font=dict(
            size=font_size,
        )
    )
    col2.plotly_chart(fig2_trend_each, use_container_width=True)

    for donor in df2_piv[color].unique():
        print_corr(df2_piv[df2_piv[color].isin([donor])], x, y, donor, corr_expander)

    correlation_between_selected_bacteria.markdown("""---""")


def st_main_ratio_between_selected_bacteria_boxplot():
    # Plot showing ratio between the two selected bacteria
    ratio_between_selected_bacteria = st.container()
    ratio_between_selected_bacteria.subheader(
        f"Ratio between {x.split(';')[-1]}:{y.split(';')[-1]}"
    )
    widget1, widget2 = ratio_between_selected_bacteria.columns(2)
    split_by = widget1.selectbox(
        label="Group By", options=[None] + meta_columns, index=0, key="split_by"
    )
    color_by2 = widget2.selectbox(
        label="Color By",
        options=[None] + meta_columns,
        index=0,
        key="color_by_ratios_plot",
    )
    box_or_strip = widget2.selectbox(
        label="Box or Strip Plot", options=["Box Plot", "Strip Plot"], index=0
    )
    if box_or_strip == "Box Plot":
        fig3 = px.box(
            df2_piv, x=split_by, y="ratio", hover_data=meta_columns, color=color_by2,
        )
        fig3.update_traces(boxmean=True)
    elif box_or_strip == "Strip Plot":
        fig3 = px.strip(
            df2_piv, x=split_by, y="ratio", hover_data=meta_columns, color=color_by2
        )
    fig3.update_xaxes(matches=None, autorange=True)
    fig3.update_layout(boxmode="group", boxgap=0)
    fig3.layout.yaxis.title = f"{x.split(';')[-1]}:{y.split(';')[-1]} ratio"
    fig3.update_layout(
        font=dict(
            size=font_size,
        )
    )
    ratio_between_selected_bacteria.plotly_chart(fig3, use_container_width=True)
    if split_by:

        values_dict={split_category:df2_piv[df2_piv[split_by]==split_category]['ratio'].dropna().tolist() for split_category in df2_piv[split_by].unique()}
        
        perms=list(combinations(values_dict.keys(),2))

        group1=[]
        group2=[]
        mean_group1=[]
        mean_group2=[]
        stat_list=[]
        p_list=[]
        statistic=ratio_between_selected_bacteria.selectbox('Choose Statistic Test',options=["Mann Whitney U",'Kruskal Wallis'])
        for perm in perms:
            group1.append(perm[0])
            group2.append(perm[1])
            mean_group1.append(np.mean(values_dict[perm[0]]))
            mean_group2.append(np.mean(values_dict[perm[1]]))
            if statistic=='Mann Whitney U':
                stat,p=mannwhitneyu(values_dict[perm[0]],values_dict[perm[1]])
            elif statistic=='Wilcoxon':
                stat,p=wilcoxon(values_dict[perm[0]],values_dict[perm[1]])
            elif statistic=='Kruskal Wallis':
                stat,p=kruskal(values_dict[perm[0]],values_dict[perm[1]])
            
            stat_list.append(stat)
            p_list.append(p)
            
        stat_df=pd.DataFrame.from_dict({'Group1':group1,'Group2':group2,'Mean_Group1':mean_group1,'Mean_Group2':mean_group2,'Stat':stat_list,'P-Value':p_list})            
        stat_df['Significant']=stat_df['P-Value']<=0.05
        ratio_between_selected_bacteria.markdown("Groups are selected by 'Group By' dropdown above the plot")
        ratio_between_selected_bacteria.write(stat_df,use_container_width=True)
            
        
        
        
    
    ratio_between_selected_bacteria.markdown("""---""")


def st_main_correlation_scatter_between_bacteria_and_metadata_parameter():
    # Scatter plot between one selected bacteria and any column in the metadata
    correlation_to_metadata_scatter = st.container()
    correlation_to_metadata_scatter.subheader(
        f"Correlate any bacteria with any metadata column"
    )
    if "SampleDay" in meta_columns:
        loc = meta_columns.index("SampleDay")
    else:
        loc = 0
    if "DonorName" in meta_columns:
        loc2 = meta_columns.index("DonorName")
    else:
        loc2 = 2
    bacteria_picker = correlation_to_metadata_scatter.selectbox(
        label="Pick Bacteria",
        options=sorter[:],
        index=0,
        key="bac_to_correlate_to_meta",
    )
    correlate_to, color, marker = correlation_to_metadata_scatter.columns(3)
    plot_type, _2, marker_size = correlation_to_metadata_scatter.columns(3)
    
    correlate_to_selection = correlate_to.selectbox(
        label="Correlate to",
        options=meta_columns,
        index=loc,
        key="correleta_to_what_meta_column",
    )  # 7
    plot_type_selction = correlate_to.selectbox(
        "Select Plot Type", options=["Dot", "Box"]
    )
    color_by_picker = color.selectbox(
        label="Color by(for the Heatmap this is the Y Axis Parameter)",
        options=[None] + meta_columns,
        index=loc2,
        key="color_by_what_meta_column",
    )
    marker_picker = marker.selectbox(
        label="Marker",
        options=[None] + meta_columns,
        index=0,
        key="marker_by_what_meta_column",
    )
    marker_size = marker.slider(
        label="Marker Size", min_value=5, max_value=30, value=8, key="marker_size"
    )
    try:
        df2_piv[correlate_to_selection] = df2_piv[correlate_to_selection].astype(float)
    except:
        pass
    df2_piv[bacteria_picker + "_no_zero"] = df2_piv[bacteria_picker].replace(0, np.nan)
    df2_piv[bacteria_picker + "_log2"] = np.log2(df2_piv[bacteria_picker + "_no_zero"])
    df2_piv[bacteria_picker + "_log10"] = np.log10(
        df2_piv[bacteria_picker + "_no_zero"]
    )
    if color_by_picker is None:
        color_seq = ["#808080", "#5A5A5A"]
    else:
        color_seq = px.colors.qualitative.Plotly
    if plot_type_selction == "Dot":
        plot = px.scatter(
            df2_piv.sort_values(by=correlate_to_selection),
            x=correlate_to_selection,
            y=bacteria_picker,
            hover_data=meta_columns,
            color=color_by_picker,
            symbol=marker_picker,
            color_discrete_sequence=color_seq,
            symbol_sequence=[
                "circle",
                "square",
                "diamond",
                "cross",
                "x",
                "triangle-up",
                "triangle-down",
                "pentagon",
                "bowtie",
            ],
        )
    elif plot_type_selction == "Box":
        plot = px.box(
            df2_piv.sort_values(by=correlate_to_selection),
            x=correlate_to_selection,
            y=bacteria_picker,
            hover_data=meta_columns,
            color=color_by_picker,
            # symbol=marker_picker,
            color_discrete_sequence=color_seq,
        )
    plot.update_layout(
        font=dict(
            size=font_size,
        ),
    )
    plot.layout.yaxis.title = plot.layout.yaxis.title["text"].split(";")[-1]
    plot.update_layout(plot_bgcolor="white")
    plot.update_traces(marker_size=marker_size)
    correlation_to_metadata_scatter.plotly_chart(plot, use_container_width=True)


    transformation = correlation_to_metadata_scatter.selectbox(
        "Choose Transformation for Heatmap",
        options=["None", "Log2", "Log10"],
        key="heatmap_trasformation",
    )
    grayscale = correlation_to_metadata_scatter.checkbox(
        "Greyscale", value=False, key="grayscale_heatmap"
    )
    boolean = correlation_to_metadata_scatter.checkbox(
        "Yes/No Heatmap", value=False, key="boolean_heatmap"
    )
    agg_func='mean'
    if grayscale:
        color_seq_heatmap = px.colors.sequential.gray_r
    else:
        color_seq_heatmap = px.colors.sequential.Blues
    if transformation == "None":
        val = bacteria_picker
    elif transformation == "Log2":
        val = bacteria_picker + "_log2"
    elif transformation == "Log10":
        val = bacteria_picker + "_log10"
    
    # fig.layout.coloraxis.colorbar.title = "log2RA"
    
    df_heatmap = df2_piv.pivot_table(
        columns=correlate_to_selection,
        index=color_by_picker,
        values=val,
        aggfunc=agg_func,
    )
    if boolean:   
        # df2_piv[val]=df2_piv[val].astype(bool)
        # agg_func=np.any
        df_heatmap=df_heatmap.replace(np.nan,0).astype(bool).astype(int)
        # st.write(df_heatmap.astype(str))
        
    df_heatmap.columns = df_heatmap.columns.astype(str)
    df_heatmap.index = df_heatmap.index.astype(str)
    fig1 = px.imshow(
        df_heatmap,
        color_continuous_scale=color_seq_heatmap,
        aspect="auto",
        title=bacteria_picker.split(";")[-1],
    )
    fig1.layout.coloraxis.colorbar.tickformat='.2E'
    fig1.update_layout(
        plot_bgcolor="white",
        autosize=True,
        font=dict(
            size=font_size,
        )
    )
    
    correlation_to_metadata_scatter.plotly_chart(fig1, use_container_width=True)


def st_main_correlation_scatter_between_ratio_and_metadata_parameter():
    # Scatter plot between one selected bacteria and any column in the metadata
    ratio_correlation_to_metadata_scatter = st.container()
    ratio_correlation_to_metadata_scatter.subheader(
        f"Correlate the Ratio With Any Metadata Parameter"
    )
    ratio_correlation_to_metadata_scatter.text(
        "Uses the above selected bacteria for the ratios"
    )
    # bacteria_picker=ratio_correlation_to_metadata_scatter.selectbox(label="Pick Bacteria",options=sorter[:],index=0)
    correlate_to1, color1 = ratio_correlation_to_metadata_scatter.columns(2)
    correlate_to_selection1 = correlate_to1.selectbox(
        label="Correlate to ",
        options=meta_columns,
        index=0,
        key="correlate_ratio_to_what_meta_column",
    )  # 7
    color_by_picker1 = color1.selectbox(
        label="Color by ",
        options=[None] + meta_columns,
        index=2,
        key="color_by_ratio_correlation_plot",
    )
    try:
        df2_piv[correlate_to_selection1] = df2_piv[correlate_to_selection1].astype(
            float
        )
    except:
        pass
    plot = px.scatter(
        df2_piv.sort_values(by=correlate_to_selection1),
        x=correlate_to_selection1,
        y="ratio",
        hover_data=meta_columns,
        color=color_by_picker1,
    )
    plot.update_layout(
        font=dict(
            size=font_size,
        )
    )
    plot.layout.yaxis.title = plot.layout.yaxis.title["text"].split(";")[-1]
    ratio_correlation_to_metadata_scatter.plotly_chart(plot, use_container_width=True)
    ratio_correlation_to_metadata_scatter.markdown("""---""")


def save_and_upload_settings():
    from json import dumps, loads
    from time import strftime

    global save_and_use_settings
    save_and_use_settings = st.sidebar.expander(
        "Save Current Settings Or Upload Saved Settings"
    )
    st.sidebar.markdown("---")
    settings_to_download = {
        k: v
        for k, v in st.session_state.items()
        if "button" not in k and "file_uploader" not in k and "FormSubmitter" not in k
    }
    custom_filename = save_and_use_settings.text_input(
        label="Choose name for the settings file",
        placeholder="Leave blank to use current date and time as the file name.",
        value="",
    )
    add_date_to_name = save_and_use_settings.checkbox(
        "Add date and time to filename", value=True
    )
    if add_date_to_name:
        timeanddate = strftime("%Y%m%d-%H%M%S")
    else:
        timeanddate = ""
    settings_filename = timestr = (
        str(timeanddate) + " " + str(custom_filename) + str(".json")
    )
    save_and_use_settings.download_button(
        label="Save Current Settings as a File",
        data=dumps(settings_to_download, default=str),
        file_name=settings_filename,
    )
    save_and_use_settings.markdown("---")

    upload_settings_widget = save_and_use_settings.file_uploader(
        label="Upload Previously Saved Settings File",
        type=["json"],
        accept_multiple_files=False,
    )

    if upload_settings_widget:
        with upload_settings_widget.getvalue() as f:
            uploaded_settings = loads(f)
            failed = []
            succeeded = []
            button_apply_uploaded_settings = save_and_use_settings.button(
                "Apply Settings",
                on_click=apply_uploaded_settings,
                args=(uploaded_settings,),
            )


def apply_uploaded_settings(json_settings):
    failed = []
    succeeded = []
    for k, v in json_settings.items():
        try:
            st.session_state[k] = v
            succeeded.append(k)
        except:
            failed.append(k)
    save_and_use_settings.success(
        f"Successfully uploaded {str(len(succeeded))} out of {str(len(succeeded) + len(failed))} settings"
    )
    if len(failed) > 0:
        save_and_use_settings.error(
            f"Failed to upload the following settings: {failed}"
        )


def main():
    # Main part of function
    st_main_header()
    # st_sidebar_pick_file()
    st_sidebar_upload_file()
    if upload_data_widget or sample_data:
        with st.spinner("Wait for it..."):
            st_main_raw_data()
            options_1 = [
                "Most Abundant Bacteria",
                "RA of All Samples",
                "Correlation Matrix of the top bacteria",
                "Top Correlations",
            ]
            options_2 = [
                "Relative Abundance of Two Selected Bacteria",
                "Correlations Between Two Selected Bacteria",
                "Ratios Between Two Selected Bacteria",
                "Correlate the Ratio With Any Metadata Paramterer",
                "Correlate Any Bacteria With Any Metadata Column",
            ]

            save_and_upload_settings()
            st_sidebar_sort_samples()
            st_sidebar_data_filters()
            st_sidebar_top_bacteria_slider()

            plot_to_show = st.radio(
                "Choose Which Plot To Show",
                options=options_1 + options_2,
                key="which_plot",
            )
            if plot_to_show in options_1:
                if plot_to_show == "Most Abundant Bacteria":
                    st_main_top_bacteria_plot()
                elif plot_to_show == "RA of All Samples":
                    show_ra_of_all()
                elif plot_to_show == "Correlation Matrix of the top bacteria":
                    st_main_correlation_heatmap_between_top_bac()
                elif plot_to_show == "Top Correlations":
                    st_main_top_correlations_plot()
            elif plot_to_show in options_2:
                st_main_bacteria_to_compare_selection()
                if plot_to_show == "Relative Abundance of Two Selected Bacteria":
                    st_main_ra_plot_of_selected_bacteria()
                elif plot_to_show == "Correlations Between Two Selected Bacteria":
                    st_main_correlation_scatter_between_selected_baceria()
                elif plot_to_show == "Ratios Between Two Selected Bacteria":
                    st_main_ratio_between_selected_bacteria_boxplot()
                elif plot_to_show == "Correlate the Ratio With Any Metadata Paramterer":
                    st_main_correlation_scatter_between_ratio_and_metadata_parameter()
                elif plot_to_show == "Correlate Any Bacteria With Any Metadata Column":
                    st_main_correlation_scatter_between_bacteria_and_metadata_parameter()
            # elif plot_to_show in options_3:


# %%
if __name__ == "__main__":
    # if check_password():
        main()
