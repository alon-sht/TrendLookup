# %%
import pandas as pd
import os
import plotly.express as px
import numpy as np
from scipy.stats import spearmanr
import streamlit as st
from PIL import Image
im = Image.open("fav.ico")

st.set_page_config(layout="wide",page_title="TrendAnalysis",page_icon=im)
hide_streamlit_style = """
            <style>
            
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True) 

def check_password():
    """Returns `True` if the user had the correct password."""

    def password_entered():
        """Checks whether a password entered by the user is correct."""
        if st.session_state["password"] == st.secrets["password"]:
            st.session_state["password_correct"] = True
            del st.session_state["password"]  # don't store password
        else:
            st.session_state["password_correct"] = False

    if "password_correct" not in st.session_state:
        # First run, show input for password.
        st.text_input(
            "Password", type="password", on_change=password_entered, key="password"
        )
        return False
    elif not st.session_state["password_correct"]:
        # Password not correct, show input + error.
        st.text_input(
            "Password", type="password", on_change=password_entered, key="password"
        )
        st.error("😕 Password incorrect")
        return False
    else:
        # Password correct.
        return True

if check_password():

    st.title("Trend Analysis")
    main=st.container()
    data=st.container()
    file_options={  'Level 3 (Class) - Stool Only':"all_donors L3.csv",
                    'Level 4 (Order) - Stool Only':"all_donors L4.csv",
                    'Level 5 (Family) - Stool Only':"all_donors L5.csv",
                    'Level 6 (Genus) - Stool Only':"all_donors L6.csv",
                    'Level 7 (Species) - Stool Only':"all_donors L7.csv",
                    'Level 3 (Class) - All Samples':'all_samples L3.csv',
                    'Level 4 (Order) - All Samples':'all_samples L4.csv',
                    'Level 5 (Family) - All Samples':'all_samples L5.csv',}
    level=main.radio(label='Pick file to use', options=list(file_options.keys()), index=0)#'6 - All Samples'

    wd=""
    try:
        file=file_options[level]
    except:
        st.warning("This file doesn't currently work")
    # elif level=='6 - All Samples':
    #     file='all_samples L6.csv'
    df=pd.read_csv(os.path.join(wd,file))#.rename(columns={"#OTU ID":"OTU",
                                          #          "0":"RA"})#[['OTU','RA',"SampleID","DonorName","ExperimentName"]]
    
    meta_columns=[x for x in df.columns.tolist() if x not in ['OTU','level_1','RA','Notes']]
    df[meta_columns]=df[meta_columns].astype(str)
    sorter=df.groupby("OTU").agg({'RA':'mean'}).sort_values("RA",ascending=False).index.tolist()
    sorterIndex = dict(zip(sorter, range(len(sorter))))
    df['rank']=df['OTU'].map(sorterIndex)
    df=df.sort_values('rank')
    data.subheader("Raw Data")
    data.write(df.astype(str),use_container_width=True)
    data.markdown("""---""")
    data.write(f"### No. of samples before filtering: {len(df['SampleID'].unique())}")
    data.markdown("""---""")
    filters=data.container()
    filters.title('Data Filters')
    query=f""
    widget_dict={}
    for col in meta_columns:
        if col not in ['SampleID','Replicates','ReplicateGroup','ZymoID']:
            widget_dict[col]=filters.multiselect(label=col,options=df[col].unique().tolist(),default=df[col].unique().tolist())
            query+=f"`{col}`  in {widget_dict[col]} & "
    df_filtered=df.query(query[:-2])
    filters.markdown("""---""")
    filters.write(f"### No. of samples after filtration {len(df_filtered['SampleID'].unique())} ({len(df_filtered['SampleID'].unique())*100/len(df['SampleID'].unique())}% of all samples)")
    filters.markdown("""---""")
    top_val=data.slider(label='Select number of bacteria to show',min_value=1,max_value=len(df_filtered['OTU'].unique()),value=10)
    df_top=df_filtered[df_filtered['OTU'].isin(sorter[:top_val])]
    data.subheader(f"Top {top_val} bacteria (ordered by mean)")
    fig_all=px.box(df_top,y='OTU',x='RA',height=900)
    fig_all.update_traces(boxmean=True,orientation='h')
    fig_all.update_yaxes(autorange="reversed",dtick=1)
    
    data.plotly_chart(fig_all,use_container_width=True)
    df1=df_filtered.pivot(index=meta_columns, columns=["OTU"], values="RA")
    df_top_corr=df_top.pivot(index=meta_columns, columns=["OTU"], values="RA").corr(method='spearman')
    df_top_corr=df_top_corr.mask(np.tril(np.ones(df_top_corr.shape),-1).astype(np.bool))
    
    corr_plot=px.imshow(df_top_corr,text_auto=".2f",template='plotly_white',color_continuous_scale='RdBu')
    corr_plot.update_layout(autosize=True,height=900)
    corr_plot.update_xaxes(showticklabels=False,showgrid=False)
    corr_plot.update_yaxes(showticklabels=False,showgrid=False)
    data.subheader("Correlation Matrix")
    data.text("Hover over the plot to see OTU names")
    data.plotly_chart(corr_plot,use_container_width=True)
    df_top_corr.index.name='OTU1'
    df_top_corr_new=df_top_corr.reset_index().melt(id_vars='OTU1',value_vars=df_top_corr.columns).sort_values('value',ascending=False).dropna()
    df_top_corr_new['value']=df_top_corr_new['value'].round(3)
    df_top_corr_new=df_top_corr_new[df_top_corr_new['value']<1]
    df_top_corr_new['x']=df_top_corr_new['OTU1'].astype(str)+df_top_corr_new['OTU'].astype(str)
    df_top_corr_new['abs']=abs(df_top_corr_new['value'])
    df_top_corr_new["Correlation"] = np.where(df_top_corr_new["value"]<0, 'Negative', 'Positive')
    data.subheader("Top Correlations")
    data.text("Hover over the plot to see OTU names")
    top_corr=data.slider(label='Number of top correlations to show',min_value=1,max_value=len(df_filtered['OTU'].unique()),value=50)
    df_top_corr_new=df_top_corr_new.sort_values('abs',ascending=False).head(top_corr).sort_values('value',ascending=False)
    top_corr_plot=px.bar(df_top_corr_new,x='x',y='value',hover_data=['OTU','OTU1'],color='Correlation')
    top_corr_plot.update_xaxes(showticklabels=False)
    data.plotly_chart(top_corr_plot,use_container_width=True)
    
    
    
    
    
    st.markdown("""---""")
    st.subheader("Select Bacteria to Compare")
    x=st.selectbox(label="Bacteria 1",options=sorter[:])
    y=st.selectbox(label="Bacteria 2",options=sorter[:],index=1)
    config = {"responsive": True, "toImageButtonOptions": {"format": "svg", }}
    df2_piv=df1.reset_index()
    df2_piv["ratio"]=df2_piv[x]/df2_piv[y]
    title=f"{x.split(';')[-1]} to {y.split(';')[-1]} ratio per donor"
    def print_corr(df,x,y,param,where):
        corr=spearmanr(df[[x,y]])
        if corr[1]<0.05:
            return where.success(f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]}")
        elif corr[1]>0.05:
            return where.error(f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]}")
        else:
            return where.info(f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]}")
    st.markdown("""---""")
    st.subheader("Relative abundance of selected bacteria")
    split_by_donor=st.checkbox("Split plot by Donor")

    if split_by_donor:
        split="DonorName"
    else:
        split=None

    fig1=px.bar(df2_piv,x="SampleID",y=[x,y],facet_col=split)
    fig1.update_xaxes(matches=None)
    fig1.update_layout(showlegend=False)
    st.plotly_chart(fig1,use_container_width=True)
    st.markdown("""---""")
    st.subheader("Correlations")
    color=st.selectbox(label="Color By",options=[None]+meta_columns,index=2)

    col1,col2=st.columns(2)
    corr_expander=st.expander(label='Spearman Correlation',expanded=False)
    fig2_trend=px.scatter(df2_piv,x=x,y=y,color=color,hover_data=meta_columns,trendline="ols",trendline_scope="overall")
    fig2_trend.layout.xaxis.title=fig2_trend.layout.xaxis.title['text'].split(';')[-1]
    fig2_trend.layout.yaxis.title=fig2_trend.layout.yaxis.title['text'].split(';')[-1]
    col1.plotly_chart(fig2_trend,use_container_width=True)    
    print_corr(df2_piv,x,y,'All Donor',corr_expander)
    fig2_trend_each=px.scatter(df2_piv,x=x,y=y,color=color,hover_data=meta_columns,trendline="ols",)
    fig2_trend_each.layout.xaxis.title=fig2_trend_each.layout.xaxis.title['text'].split(';')[-1]
    fig2_trend_each.layout.yaxis.title=fig2_trend_each.layout.yaxis.title['text'].split(';')[-1]

    col2.plotly_chart(fig2_trend_each,use_container_width=True)        
    for donor in df2_piv[color].unique():
        print_corr(df2_piv[df2_piv[color].isin([donor])],x,y,donor,corr_expander)
        
    st.subheader(f"{x.split(';')[-1]}:{y.split(';')[-1]} ratio")
    widget1,widget2=st.columns(2)
    split_by=widget1.selectbox(label="Group By",options=[None]+meta_columns,index=0)
    color_by2=widget2.selectbox(label="Color By",options=[None]+meta_columns,index=0)
    fig3=px.box(df2_piv,x=split_by,y="ratio",hover_data=meta_columns,color=color_by2,)
    fig3.update_xaxes(matches=None,autorange=True)
    fig3.update_layout(boxmode='group',boxgap=0)
    fig3.layout.yaxis.title=f"{x.split(';')[-1]}:{y.split(';')[-1]} ratio"

    st.plotly_chart(fig3,use_container_width=True)
    
    bacteria_picker=st.selectbox(label="Pick Bacteria",options=sorter[:],index=0)
    correlate_to,color=st.columns(2)
    correlate_to_selection=correlate_to.selectbox(label="Correlate to",options=meta_columns,index=0)#7
    color_by_picker=color.selectbox(label="Color by",options=[None]+meta_columns,index=2)
    try:
        df2_piv[correlate_to_selection]=df2_piv[correlate_to_selection].astype(float)
    except:
        pass
    xx=px.scatter(df2_piv.sort_values(by=correlate_to_selection),x=correlate_to_selection,y=bacteria_picker,hover_data=meta_columns,color=color_by_picker)
    #trendline="ols",trendline_scope="overall")
    # fig4_trend=px.scatter(df2_piv,x=correlate_to,y=bac_picker,)#hover_data=meta_columns)#,trendline="ols",trendline_scope="overall")
    # fig4_trend.layout.xaxis.title=fig4_trend.layout.xaxis.title['text'].split(';')[-1]
    xx.layout.yaxis.title=xx.layout.yaxis.title['text'].split(';')[-1]
    st.plotly_chart(xx,use_container_width=True)  
    

# %%
