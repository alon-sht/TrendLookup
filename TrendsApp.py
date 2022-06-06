# %%
import os
import pandas as pd
import plotly.express as px
import numpy as np
from scipy.stats import spearmanr
import streamlit as st
from PIL import Image
from src.check_password import check_password

st.set_page_config(layout="wide",page_title="TrendAnalysis",page_icon=Image.open("fav.ico"))
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
    st.title("Trend Analysis")
    
    
def st_sidebar_pick_file():
    #Sidebar dropdown to pick file to look at 
    global level,df
    st.sidebar.image("Mybiotics_LOGO - Large.png",width=350)
    st.sidebar.markdown("""---""")
    file_select=st.sidebar.container()
    file_options={  'Level 3 (Class) - Stool Only':"all_donors L3.csv",
                    'Level 4 (Order) - Stool Only':"all_donors L4.csv",
                    'Level 5 (Family) - Stool Only':"all_donors L5.csv",
                    'Level 6 (Genus) - Stool Only':"all_donors L6.csv",
                    # 'Level 7 (Species) - Stool Only':"all_donors L7.csv",
                    'Level 3 (Class) - All Samples':'all_samples L3.csv',
                    'Level 4 (Order) - All Samples':'all_samples L4.csv',
                    'Level 5 (Family) - All Samples':'all_samples L5.csv',
                    }
    
    level=file_select.selectbox(label='Pick file to use', options=list(file_options.keys()), index=0)#'6 - All Samples'
    file_select.markdown("""---""")
    wd=''
    try:
        file=file_options[level]
        df=pd.read_csv(os.path.join(wd,file))
    except:
        st.warning("This file doesn't currently work")
    
    
    
def st_main_raw_data():
    # Show raw data table in main container
    global meta_columns,sorter,data
    
    raw_data=st.container()
    meta_columns=[x for x in df.columns.tolist() if x not in ['OTU','level_1','RA','Notes']]
    df[meta_columns]=df[meta_columns].astype(str)
    sorter_mean=df.groupby("OTU").agg({'RA':'mean'}).sort_values("RA",ascending=False).index.tolist()
    sorter_median=df.groupby("OTU").agg({'RA':'mean'}).sort_values("RA",ascending=False).index.tolist()
    sorter_choose=st.sidebar.radio("Sort samples by Mean or Median",options=["Mean","Median"],index=0)#,horizontal=True)
    st.sidebar.markdown("""---""")
    if sorter_choose=="Mean":
        sorter=sorter_mean
    elif sorter_choose=="Median":
        sorter=sorter_median
    sorterIndex = dict(zip(sorter, range(len(sorter))))
    df['ind']=df['OTU'].map(sorterIndex)
    df.sort_values(by='ind',inplace=True)
    raw_data.subheader("Raw Data")
    raw_data.write(df.astype(str),use_container_width=True)
    raw_data.markdown("""---""")
    # data=st.container()
    
    
    
def st_sidebar_data_filters():
    # Show data filters in the sidebar inside an expander
    # Filters are updated with the press of a button (and not automatically)
    # Also number of samples before and after filtering is shown
    
    global df_filtered
    
    filters=st.sidebar.container()
    filters.subheader('Data Filters')
    filter_widgets=filters.expander("Filter Widgets (Click to Expand). After selecting filters click the button at the bottom.")
    query=f""
    widget_dict={}
    form=filter_widgets.form("form1")
    for col in meta_columns:
        if col not in ['SampleID','Replicates','ReplicateGroup','ZymoID']:
            widget_dict[col]=form.multiselect(label=col,options=df[col].unique().tolist(),default=df[col].unique().tolist())
            query+=f"`{col}`  in {widget_dict[col]} & "
    
    submit=form.form_submit_button("Filter Data")
    df_filtered=df.query(query[:-2])
    filters.markdown(f"##### No. of samples before filtering: {len(df['SampleID'].unique())}")
    filters.markdown(f"##### No. of samples after filtering: {len(df_filtered['SampleID'].unique())} ({len(df_filtered['SampleID'].unique())*100/len(df['SampleID'].unique())}% of all samples)")
    filters.markdown("""---""")



def st_sidebar_top_bacteria_slider():
    #Slider to select number of top bacteria to show
    global top_val, df_top
    top_val=st.sidebar.slider(label='Select number of bacteria to show',min_value=1,max_value=len(df_filtered['OTU'].unique()),value=10)
    df_top=df_filtered[df_filtered['OTU'].isin(sorter[:top_val])]
    

def st_main_top_bacteria_plot():
    # A box plot showing the relative abundance of the top bacteria 
    top_bacteria_boxplot=st.container()
    top_bacteria_boxplot.subheader(f"Top {top_val} bacteria (ordered by mean)")
    top_bacteria_boxplot.text("Change number of top bacteria by using the slider in the sidebar")
    fig_all=px.box(df_top,y='OTU',x='RA',height=900,)#template='plotly_white')
    fig_all.update_traces(boxmean=True,orientation='h',)
    fig_all.update_yaxes(autorange="reversed",dtick=1,)
    fig_all.update_xaxes(title='Relative Abundance')
    top_bacteria_boxplot.plotly_chart(fig_all,use_container_width=True)
    top_bacteria_boxplot.markdown("""---""")
    
    
def st_main_correlation_heatmap_between_top_bac():
    # Show heatmap of top correlations
    global df_top_corr, df1
    top_bac_correlation_heatmap=st.container()
    df1=df_filtered.pivot(index=meta_columns, columns=["OTU"], values="RA")
    df_top_corr=df_top.pivot(index=meta_columns, columns=["OTU"], values="RA").corr(method='spearman')
    df_top_corr=df_top_corr.mask(np.tril(np.ones(df_top_corr.shape),-1).astype(bool))
    
    corr_plot=px.imshow(df_top_corr,text_auto=".2f",template='plotly_white',color_continuous_scale='RdBu')
    corr_plot.update_layout(autosize=True,height=900)
    corr_plot.update_xaxes(showticklabels=False,showgrid=False)
    corr_plot.update_yaxes(showticklabels=False,showgrid=False)
    top_bac_correlation_heatmap.subheader(f"Correlation Matrix of the top {top_val} bacteria")
    top_bac_correlation_heatmap.text("Hover over the plot to see OTU names")
    top_bac_correlation_heatmap.text("Change number of shown bacteria by using the slider in the sidebar")
    top_bac_correlation_heatmap.plotly_chart(corr_plot,use_container_width=True)
    top_bac_correlation_heatmap.markdown("""---""")
    
    
    
def st_main_top_correlations_plot():
    #Show bar chart of main correlations
    top_correlation_plot=st.container()
    df_top_corr.index.name='OTU1'
    df_top_corr_new=df_top_corr.reset_index().melt(id_vars='OTU1',value_vars=df_top_corr.columns).sort_values('value',ascending=False).dropna()
    df_top_corr_new['value']=df_top_corr_new['value'].round(3)
    df_top_corr_new=df_top_corr_new[df_top_corr_new['value']<1]
    df_top_corr_new['x']=df_top_corr_new['OTU1'].astype(str)+df_top_corr_new['OTU'].astype(str)
    df_top_corr_new['abs']=abs(df_top_corr_new['value'])
    df_top_corr_new["Correlation"] = np.where(df_top_corr_new["value"]<0, 'Negative', 'Positive')
    top_corr=st.sidebar.slider(label='Number of top correlations to show',min_value=1,max_value=len(df_filtered['OTU'].unique()),value=50)
    top_correlation_plot.subheader(f"{top_corr} Top Correlations")
    top_correlation_plot.text("Hover over the plot to see OTU names")
    top_correlation_plot.text("Change number of correlations shown by using the slider in the sidebar")
    df_top_corr_new=df_top_corr_new.sort_values('abs',ascending=False).head(top_corr).sort_values('value',ascending=False)
    top_corr_plot=px.bar(df_top_corr_new,x='x',y='value',hover_data=['OTU','OTU1'],color='Correlation',text_auto=True)
    top_corr_plot.update_traces(textfont_size=12, textangle=0, textposition="outside", )#cliponaxis=False)
    top_corr_plot.update_xaxes(showticklabels=False,title='OTU Pair')
    top_corr_plot.update_yaxes(title='Correlation')
    top_correlation_plot.plotly_chart(top_corr_plot,use_container_width=True)
    top_correlation_plot.markdown("""---""")
    
def st_main_bacteria_to_compare_selection():
    # Selection to bacteria to compare
    global x,y
    compare_bac_selection=st.container()
    compare_bac_selection.subheader("Select Bacteria to Compare")
    compare_bac_selection.text("Dropdowns are sorted by mean relative abundance")
    x=compare_bac_selection.selectbox(label="Bacteria 1",options=sorter[:])
    y=compare_bac_selection.selectbox(label="Bacteria 2",options=sorter[:],index=1)
    compare_bac_selection.markdown("""---""")
    
def print_corr(df,x,y,param,where):
    # Function to print correlations
    corr=spearmanr(df[[x,y]])
    if corr[1]<0.05:
        return where.success(f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]}")
    elif corr[1]>0.05:
        return where.error(f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]}")
    else:
        return where.info(f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]}")
    
    
def st_main_ra_plot_of_selected_bacteria():
    # stacked bar chart of RA of the two selected bacteria
    
    global df2_piv
    ra_of_selected_bacteria=st.container()
    df2_piv=df1.reset_index()
    df2_piv["ratio"]=df2_piv[x]/df2_piv[y]

    
    ra_of_selected_bacteria.subheader("Relative abundance of selected bacteria")
    split_by_donor=ra_of_selected_bacteria.checkbox("Split plot by Donor")

    if split_by_donor:
        split="DonorName"
    else:
        split=None

    fig1=px.bar(df2_piv,x="SampleID",y=[x,y],facet_col=split)
    fig1.update_xaxes(matches=None)
    fig1.update_yaxes(title="Relative Abundance")
    fig1.update_layout(showlegend=False)
    ra_of_selected_bacteria.plotly_chart(fig1,use_container_width=True)
    ra_of_selected_bacteria.markdown("""---""")


def st_main_correlation_scatter_between_selected_baceria():
    #Scatter plot of RA of the two selected bacteria    
    correlation_between_selected_bacteria=st.container()
    correlation_between_selected_bacteria.subheader("Spearman correlations between two selected bacteria")
    correlation_between_selected_bacteria.text(f"Correlation between {x.split(';')[-1]} and {y.split(';')[-1]} - overall and for each color group")
    color=correlation_between_selected_bacteria.selectbox(label="Color By",options=[None]+meta_columns,index=2)

    col1,col2=correlation_between_selected_bacteria.columns(2)
    corr_expander=correlation_between_selected_bacteria.expander(label='Spearman Correlation (Click to Expand)',expanded=False)
    fig2_trend=px.scatter(df2_piv,x=x,y=y,color=color,hover_data=meta_columns,trendline="ols",trendline_scope="overall")
    fig2_trend.layout.xaxis.title=fig2_trend.layout.xaxis.title['text'].split(';')[-1]
    fig2_trend.layout.yaxis.title=fig2_trend.layout.yaxis.title['text'].split(';')[-1]
    col1.plotly_chart(fig2_trend,use_container_width=True)    
    print_corr(df2_piv,x,y,'Overall',corr_expander)
    fig2_trend_each=px.scatter(df2_piv,x=x,y=y,color=color,hover_data=meta_columns,trendline="ols",)
    fig2_trend_each.layout.xaxis.title=fig2_trend_each.layout.xaxis.title['text'].split(';')[-1]
    fig2_trend_each.layout.yaxis.title=fig2_trend_each.layout.yaxis.title['text'].split(';')[-1]
    col2.plotly_chart(fig2_trend_each,use_container_width=True)        
        
    for donor in df2_piv[color].unique():
        print_corr(df2_piv[df2_piv[color].isin([donor])],x,y,donor,corr_expander)
    
    
    correlation_between_selected_bacteria.markdown("""---""")
    
    
    
def st_main_ratio_between_selected_bacteria_boxplot():
    # Plot showing ratio between the two selected bacteria
    ratio_between_selected_bacteria=st.container()
    ratio_between_selected_bacteria.subheader(f"Ratio between {x.split(';')[-1]}:{y.split(';')[-1]}")
    widget1,widget2=ratio_between_selected_bacteria.columns(2)
    split_by=widget1.selectbox(label="Group By",options=[None]+meta_columns,index=0)
    color_by2=widget2.selectbox(label="Color By",options=[None]+meta_columns,index=0)
    box_or_strip=widget2.selectbox(label="Box or Strip Plot",options=['Box Plot','Strip Plot'],index=0)
    if box_or_strip=="Box Plot":
        fig3=px.box(df2_piv,x=split_by,y="ratio",hover_data=meta_columns,color=color_by2)
    elif box_or_strip=="Strip Plot":
        fig3=px.strip(df2_piv,x=split_by,y="ratio",hover_data=meta_columns,color=color_by2)
    fig3.update_xaxes(matches=None,autorange=True)
    fig3.update_layout(boxmode='group',boxgap=0)
    fig3.layout.yaxis.title=f"{x.split(';')[-1]}:{y.split(';')[-1]} ratio"

    ratio_between_selected_bacteria.plotly_chart(fig3,use_container_width=True)
    ratio_between_selected_bacteria.markdown("""---""")
    
    
    
    
    
def st_main_correlation_scatter_between_bacteria_and_metadata_parameter():
    # Scatter plot between one selected bacteria and any column in the metadata 
    correlation_to_metadata_scatter=st.container()
    correlation_to_metadata_scatter.subheader(f"Correlate any bacteria with any metadata column")
    bacteria_picker=correlation_to_metadata_scatter.selectbox(label="Pick Bacteria",options=sorter[:],index=0)
    correlate_to,color=correlation_to_metadata_scatter.columns(2)
    correlate_to_selection=correlate_to.selectbox(label="Correlate to",options=meta_columns,index=0)#7
    color_by_picker=color.selectbox(label="Color by",options=[None]+meta_columns,index=2)
    try:
        df2_piv[correlate_to_selection]=df2_piv[correlate_to_selection].acorrelation_to_metadata_scatterype(float)
    except:
        pass
    plot=px.scatter(df2_piv.sort_values(by=correlate_to_selection),x=correlate_to_selection,y=bacteria_picker,hover_data=meta_columns,color=color_by_picker)
    plot.layout.yaxis.title=plot.layout.yaxis.title['text'].split(';')[-1]
    correlation_to_metadata_scatter.plotly_chart(plot,use_container_width=True)


    
def main():
    #Main part of function
    st_main_header()
    st_sidebar_pick_file()
    st_main_raw_data()
    st_sidebar_data_filters()
    st_sidebar_top_bacteria_slider()
    st_main_top_bacteria_plot()
    st_main_correlation_heatmap_between_top_bac()
    st_main_top_correlations_plot()
    st_main_bacteria_to_compare_selection()
    st_main_ra_plot_of_selected_bacteria()
    st_main_correlation_scatter_between_selected_baceria()
    st_main_ratio_between_selected_bacteria_boxplot()
    st_main_correlation_scatter_between_bacteria_and_metadata_parameter()
    
    
    
    

# %%
if __name__=='__main__':
    if check_password():
        main()  