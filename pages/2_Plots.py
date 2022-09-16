# %%
import streamlit as st
import plotly.express as px
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu,kruskal,wilcoxon
from itertools import combinations,combinations
from src.functions import update_session_state

from src.data_functions import sort_samples

update_session_state(update_all=True)

# %%
# st.write(st.session_state
# 
def general_plot_settings():
    plot_settings=st.sidebar.expander("Plot Settings")
    
    font_size = plot_settings.slider(
        label="Font Size", min_value=1, max_value=30, value=12, key="font_size"
    )
    sorter=plot_settings.radio("Sort by",options=['Mean','Median'],key='sorter_radio',horizontal=True)
    st.session_state['sorter']=st.session_state['sorter_mean'] if sorter=="Mean" else st.session_state['sorter_median']
    
def show_df():
    st.write(st.session_state['filtered_df'].astype(str))



def show_ra_of_all():
    ra_all = st.container()
    ra_all.subheader("Relative Abundance of All Samples")
    ra_all.text("For many samples (>30) this might be very slow")
    ra_all.text(
        f"You currently have {len(st.session_state['filtered_df']['SampleID'].unique())} samples after applying filters with {len(st.session_state['filtered_df']['OTU'].unique())} unique OTUs"
    )
    load_plot = ra_all.checkbox("Load Plot", value=False)
    ra_all_plot = ra_all.container()
    if load_plot:
        ra_all_legend = ra_all_plot.checkbox("Show Legend on Plot", value=False)
        with st.spinner("Loading plot"):
            fig = px.bar(
                st.session_state['filtered_df'],
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

    ra_all.markdown("""---""")


def top_bacteria_plot():
    # A box plot showing the relative abundance of the top bacteria
    
    
    top_bacteria_boxplot = st.container()
    
    top_bacteria_boxplot.subheader(f"Most abundant bacteria")
    top_val = top_bacteria_boxplot.slider(
        label="Select number of bacteria to show",
        min_value=1,
        max_value=len(st.session_state['filtered_df']["OTU"].unique()),
        value=10,
        key="top_val",
    )
    df_top = st.session_state['filtered_df'][st.session_state['filtered_df']["OTU"].isin(st.session_state['sorter'][:top_val])]
    with st.spinner("Loading"):
        fig_all = px.box(
            df_top,
            y="OTU",
            x="RA",
            height=900,
        )  # template='plotly_white')
        fig_all.update_layout(
            font=dict(
                size=st.session_state['font_size'],
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
    
    
def correlation_heatmap_between_top_bac():
    # Show heatmap of top correlations
    top_bac_correlation_heatmap = st.container()

    top_val = top_bac_correlation_heatmap.slider(
        label="Select number of bacteria to correlate (sorted by mean)",
        min_value=1,
        max_value=len(st.session_state['filtered_df']["OTU"].unique()),
        value=10,
        key="top_val1",
    )
    
    df_top = st.session_state['filtered_df'][st.session_state['filtered_df']["OTU"].isin(st.session_state['sorter'][:top_val])]
    df_top_corr = df_top.pivot(index=st.session_state['meta_columns'], columns=["OTU"], values="RA").corr(
        method="spearman"
    )
    df_top_corr = df_top_corr.mask(np.tril(np.ones(df_top_corr.shape), -1).astype(bool))
    

    selection_corr=top_bac_correlation_heatmap.radio("Select Plot",options=['Correlation heatmap','Top Correlations'])
    
    if selection_corr=='Correlation heatmap':
        top_bac_correlation_heatmap.subheader(
            f"Correlation Matrix of the most abundant bacteria"
        )
        top_bac_correlation_heatmap.text("Hover over the plot to see OTU names")
        top_bac_correlation_heatmap.text(
            "Change number of shown bacteria by using the slider in the sidebar"
        )
        with st.spinner("Loading"):
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
                    size=st.session_state['font_size'],
                )
            )

            top_bac_correlation_heatmap.plotly_chart(corr_plot, use_container_width=True)
    elif selection_corr=='Top Correlations':
        st_main_top_correlations_plot(df_top_corr)
    # top_bac_correlation_heatmap.markdown("""---""")


def st_main_top_correlations_plot(df_top_corr):
    # Show bar chart of main correlations
    top_correlation_plot_container = st.container()
    
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
    top_correlation_plot_container.subheader(f"Top Correlations")
    top_correlation_plot_container.text("Hover over the plot to see OTU names")
    top_correlation_plot_container.text(
        "Change number of correlations shown by using the slider in the sidebar"
    )
    top_corr = top_correlation_plot_container.slider(
        label="Number of top correlations to show",
        min_value=1,
        max_value=len(st.session_state['filtered_df']["OTU"].unique()),
        value=5,
        key='number_of_top_corr'
    )
    df_top_corr_new = (
        df_top_corr_new.sort_values("abs", ascending=False)
        .head(top_corr)
        .sort_values("value", ascending=False)
    )
    text_on_plot = top_correlation_plot_container.checkbox("Show text on plot", value=False)
    top_corr_plot = px.bar(
        df_top_corr_new,
        x="x",
        y="value",
        hover_data=["OTU", "OTU1"],
        color="Correlation",
        text_auto=text_on_plot,
        # orientation='h'
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
            size=st.session_state['font_size']
        )
    )
    top_correlation_plot_container.plotly_chart(top_corr_plot, use_container_width=True)
    top_correlation_plot_container.markdown("""---""")
    
    

def ra_plot_of_selected_bacteria(df,x,y):
    # stacked bar chart of RA of the two selected bacteria

    ra_of_selected_bacteria = st.container()

    ra_of_selected_bacteria.subheader(
        f"Relative abundance of {x.split(';')[-1]} and {y.split(';')[-1]}"
    )
    split=ra_of_selected_bacteria.selectbox("Split By",options=[None]+st.session_state['meta_columns'])
    
    fig1 = px.bar(df, x="SampleID", y=[x, y], facet_col=split)
    fig1.update_xaxes(matches=None)
    fig1.update_layout(showlegend=False)
    fig1.update_layout(
        font=dict(
            size=st.session_state['font_size'],
        )
    )
    fig1.for_each_annotation(lambda a: a.update(text=a.text.replace(str(split)+"=", "")))
    ra_of_selected_bacteria.plotly_chart(fig1, use_container_width=True)
    ra_of_selected_bacteria.markdown("""---""")
    
    
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

def correlation_scatter_between_selected_baceria(df,x,y):
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
        options=[None] + st.session_state['meta_columns'],
        index=2,
        key="color_by_corr_bacteria_plot",
    )

    col1, col2 = correlation_between_selected_bacteria.columns(2)
    corr_expander = correlation_between_selected_bacteria.expander(
        label="Spearman Correlation (Click to Expand)", expanded=False
    )
    fig2_trend = px.scatter(
        df,
        x=x,
        y=y,
        color=color,
        hover_data=st.session_state['meta_columns'],
        trendline="ols",
        trendline_scope="overall",
    )
    fig2_trend.layout.xaxis.title = fig2_trend.layout.xaxis.title["text"].split(";")[-1]
    fig2_trend.layout.yaxis.title = fig2_trend.layout.yaxis.title["text"].split(";")[-1]
    fig2_trend.update_layout(
        font=dict(
            size=st.session_state['font_size'],
        )
    )
    col1.plotly_chart(fig2_trend, use_container_width=True)
    print_corr(df, x, y, "Overall", corr_expander)
    fig2_trend_each = px.scatter(
        df,
        x=x,
        y=y,
        color=color,
        hover_data=st.session_state['meta_columns'],
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
            size=st.session_state['font_size'],
        )
    )
    col2.plotly_chart(fig2_trend_each, use_container_width=True)

    for donor in df[color].unique():
        print_corr(df[df[color].isin([donor])], x, y, donor, corr_expander)

    correlation_between_selected_bacteria.markdown("""---""")

def correlate_two_bacteria():
    
    compare_bac_selection = st.container()
    compare_bac_selection.subheader("Select Bacteria to Compare")
    compare_bac_selection.text("Dropdowns are sorted by mean relative abundance")
    with compare_bac_selection.form("Select Bacteria to Compare"):
        x = st.selectbox(
            label="Bacteria 1", options=st.session_state['sorter'][:], key="bac1"
        )
        y = st.selectbox(
            label="Bacteria 2", options=st.session_state['sorter'][:], index=1, key="bac2"
        )
        st.form_submit_button("Apply")
    compare_bac_selection.markdown("""---""")
    # global df2_piv
    df1 = st.session_state['filtered_df'].pivot(index=st.session_state['meta_columns'], columns=["OTU"], values="RA")
    df2_piv = df1.reset_index()
    df2_piv["ratio"] = df2_piv[x] / df2_piv[y]
    df2_piv["ratio"]=df2_piv["ratio"].replace(np.inf, np.nan)
    
    cor_two_radio=st.radio("Choose plot",options=['RA of selected bacteria','Correlation of selected bacteria',])
    
    if cor_two_radio=='RA of selected bacteria':
        ra_plot_of_selected_bacteria(df2_piv,x,y)
    elif cor_two_radio=='Correlation of selected bacteria':
        correlation_scatter_between_selected_baceria(df2_piv,x,y)
        

def ratio_between_selected_bacteria_boxplot(df2_piv,x,y):
    # Plot showing ratio between the two selected bacteria
    ratio_between_selected_bacteria = st.container()
    ratio_between_selected_bacteria.subheader(
        f"Ratio between {x.split(';')[-1]}:{y.split(';')[-1]}"
    )
    widget1, widget2 = ratio_between_selected_bacteria.columns(2)
    split_by = widget1.selectbox(
        label="Group By", options=[None] + st.session_state['meta_columns'], index=0, key="split_by"
    )
    color_by2 = widget2.selectbox(
        label="Color By",
        options=[None] + st.session_state['meta_columns'],
        index=0,
        key="color_by_ratios_plot",
    )
    box_or_strip = widget2.selectbox(
        label="Box or Strip Plot", options=["Box Plot", "Strip Plot"], index=0
    )
    if box_or_strip == "Box Plot":
        fig3 = px.box(
            df2_piv, x=split_by, y="ratio", hover_data=st.session_state['meta_columns'], color=color_by2,
        )
        fig3.update_traces(boxmean=True)
    elif box_or_strip == "Strip Plot":
        fig3 = px.strip(
            df2_piv, x=split_by, y="ratio", hover_data=st.session_state['meta_columns'], color=color_by2
        )
    fig3.update_xaxes(matches=None, autorange=True)
    fig3.update_layout(boxmode="group", boxgap=0)
    fig3.layout.yaxis.title = f"{x.split(';')[-1]}:{y.split(';')[-1]} ratio"
    fig3.update_layout(
        font=dict(
            size=st.session_state['font_size']
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
        median_group1=[]
        median_group2=[]
        stat_list=[]
        p_list=[]
        statistic=ratio_between_selected_bacteria.selectbox('Choose Statistic Test',options=["Mann Whitney U",'Kruskal Wallis'])
        for perm in perms:
            group1.append(perm[0])
            group2.append(perm[1])
            mean_group1.append(np.mean(values_dict[perm[0]]))
            mean_group2.append(np.mean(values_dict[perm[1]]))
            median_group1.append(np.median(values_dict[perm[0]]))
            median_group2.append(np.median(values_dict[perm[1]]))
            if statistic=='Mann Whitney U':
                stat,p=mannwhitneyu(values_dict[perm[0]],values_dict[perm[1]])
            elif statistic=='Wilcoxon':
                stat,p=wilcoxon(values_dict[perm[0]],values_dict[perm[1]])
            elif statistic=='Kruskal Wallis':
                stat,p=kruskal(values_dict[perm[0]],values_dict[perm[1]])
            
            stat_list.append(stat)
            p_list.append(p)
            
        stat_df=pd.DataFrame.from_dict({'Group1':group1,'Group2':group2,'Mean_Group1':mean_group1,'Mean_Group2':mean_group2,'Median_Group1':median_group1,'Median_Group2':median_group2,'Stat':stat_list,'P-Value':p_list})            
        stat_df['Significant']=stat_df['P-Value']<=0.05
        # stat_df['Which is heigher']=str(split_by)+": "
        ratio_between_selected_bacteria.markdown("Groups are selected by 'Group By' dropdown above the plot")
        ratio_between_selected_bacteria.write(stat_df,use_container_width=True)
            
        
        
        
    
    ratio_between_selected_bacteria.markdown("""---""")


def correlation_scatter_between_bacteria_and_metadata_parameter():
    # Scatter plot between one selected bacteria and any column in the metadata
    correlation_to_metadata_scatter = st.container()
    df1 = st.session_state['filtered_df'].pivot(index=st.session_state['meta_columns'], columns=["OTU"], values="RA")

    df2_piv = df1.reset_index()
    correlation_to_metadata_scatter.subheader(
        f"Correlate any bacteria with any metadata column"
    )
    if "SampleDay" in st.session_state['meta_columns']:
        loc = st.session_state['meta_columns'].index("SampleDay")
    else:
        loc = 0
    if "DonorName" in st.session_state['meta_columns']:
        loc2 = st.session_state['meta_columns'].index("DonorName")
    else:
        loc2 = 2
    bacteria_picker = correlation_to_metadata_scatter.selectbox(
        label="Pick Bacteria",
        options=st.session_state['sorter'][:],
        index=0,
        key="bac_to_correlate_to_meta",
    )
    correlate_to, color, marker = correlation_to_metadata_scatter.columns(3)
    plot_type, _2, marker_size = correlation_to_metadata_scatter.columns(3)
    
    correlate_to_selection = correlate_to.selectbox(
        label="Correlate to",
        options=st.session_state['meta_columns'],
        index=loc,
        key="correleta_to_what_meta_column",
    )  # 7
    plot_type_selction = correlate_to.selectbox(
        "Select Plot Type", options=["Dot", "Box"]
    )
    color_by_picker = color.selectbox(
        label="Color by(for the Heatmap this is the Y Axis Parameter)",
        options=[None] + st.session_state['meta_columns'],
        index=loc2,
        key="color_by_what_meta_column",
    )
    marker_picker = marker.selectbox(
        label="Marker",
        options=[None] + st.session_state['meta_columns'],
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
            hover_data=st.session_state['meta_columns'],
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
            hover_data=st.session_state['meta_columns'],
            color=color_by_picker,
            # symbol=marker_picker,
            color_discrete_sequence=color_seq,
        )
    plot.update_layout(
        font=dict(
            size=st.session_state['font_size']
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
    change=correlation_to_metadata_scatter.checkbox("calculate change")
    if change:
        if boolean:
            correlation_to_metadata_scatter.error("Disable YES/NO heatmap")
            st.stop()
        ref=correlation_to_metadata_scatter.selectbox("select reference",options=df_heatmap.columns.tolist())
        df_heatmap_change=df_heatmap.subtract(df_heatmap[ref],axis='index')
        
    df=df_heatmap_change if change else df_heatmap
    fig1 = px.imshow(
        df,
        x=df.columns,
        y=df.index,
        color_continuous_scale= px.colors.sequential.RdBu if change else color_seq_heatmap,
        color_continuous_midpoint=0 if change else None,
        aspect="auto",
        title=bacteria_picker.split(";")[-1],
    )
    fig1.layout.coloraxis.colorbar.tickformat='.2E'
    fig1.update_layout(
        plot_bgcolor="white",
        autosize=True,
        font=dict(
            size=st.session_state['font_size']
        )
    )

    correlation_to_metadata_scatter.plotly_chart(fig1, use_container_width=True)


def correlation_scatter_between_ratio_and_metadata_parameter(df2_piv):
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
        options=st.session_state['meta_columns'],
        index=0,
        key="correlate_ratio_to_what_meta_column",
    )  # 7
    color_by_picker1 = color1.selectbox(
        label="Color by ",
        options=[None] + st.session_state['meta_columns'],
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
        hover_data=st.session_state['meta_columns'],
        color=color_by_picker1,
    )
    plot.update_layout(
        font=dict(
            size=st.session_state['font_size']
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
        
def bacteria_ratios():
    ratios_bac_selection = st.container()
    ratios_bac_selection.subheader("Select Bacteria to Compare")
    with ratios_bac_selection.form("Select Bacteria to Compare"):
        x = st.selectbox(
            label="Bacteria 1", options=st.session_state['sorter'][:], key="bac1.1"
        )
        y = st.selectbox(
            label="Bacteria 2", options=st.session_state['sorter'][:], index=1, key="bac2.1"
        )
        st.form_submit_button("Apply")
        
        
    df1 = st.session_state['filtered_df'].pivot(index=st.session_state['meta_columns'], columns=["OTU"], values="RA")
    df2_piv = df1.reset_index()
    df2_piv["ratio"] = df2_piv[x] / df2_piv[y]
    df2_piv["ratio"]=df2_piv["ratio"].replace(np.inf, np.nan)
    cor_ratio_radio=st.radio("Choose plot",options=['Ratios Boxplot','Correlation of ratio to metadata',])
    
    if cor_ratio_radio=='Ratios Boxplot':
        # ra_plot_of_selected_bacteria(df2_piv,x,y)
        ratio_between_selected_bacteria_boxplot(df2_piv,x,y)
    elif cor_ratio_radio=='Correlation of ratio to metadata':
        # correlation_scatter_between_selected_baceria(df2_piv,x,y)
        correlation_scatter_between_ratio_and_metadata_parameter(df2_piv)
    
# %%

menu={
    'DataFrame':show_df,
    'Relative Abundance':show_ra_of_all,
    'Most Abundant Bacteria':top_bacteria_plot,
    'Correlation between most abundant bacteria':correlation_heatmap_between_top_bac,
    'Correlation between two bacteria':correlate_two_bacteria,
    'Ratios between two bacteria':bacteria_ratios,
    'Correlate Bacterium to Metadata':correlation_scatter_between_bacteria_and_metadata_parameter
    
}
def main():

    if 'filtered_df' in st.session_state:
        menu_radio=st.sidebar.radio("Menu",options=menu.keys())
        general_plot_settings()
        menu[menu_radio]()
        
if __name__=='__main__':
        main()