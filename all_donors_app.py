# %%
import pandas as pd
import os
import plotly.express as px

from scipy.stats import spearmanr
import streamlit as st
st.set_page_config(layout="wide")
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
    st.title("Correlations")
    main=st.container()
    data=st.container()


    wd=""
    file="all_donors L6.csv"
    df=pd.read_csv(os.path.join(wd,file)).rename(columns={"#OTU ID":"OTU",
                                                    "0":"RA"})[['OTU','RA','ID','Donor','Experiment']]
    sorter=df.groupby("OTU").agg({'RA':'mean'}).sort_values("RA",ascending=False).index.tolist()
    sorterIndex = dict(zip(sorter, range(len(sorter))))
    df['rank']=df['OTU'].map(sorterIndex)
    df=df.sort_values('rank')
    data.subheader("Raw Data")
    data.write(df.astype(str),use_container_width=True)
    data.markdown("""---""")
    top_val=data.slider(label='Select number of bacteria to show',min_value=1,max_value=len(df['OTU'].unique()),value=10)
    df_top=df[df['OTU'].isin(sorter[:top_val])]
    data.subheader(f"Top {top_val} bacteria (ordered by mean)")
    fig_all=px.box(df_top,y='OTU',x='RA',height=900)
    fig_all.update_traces(boxmean=True,orientation='h')
    fig_all.update_yaxes(autorange="reversed",dtick=1)
    data.plotly_chart(fig_all,use_container_width=True)

    st.markdown("""---""")
    st.subheader("Select Bacteria to Compare")
    x=st.selectbox(label="Bacteria 1",options=sorter[:top_val])
    y=st.selectbox(label="Bacteria 2",options=sorter[:top_val],index=1)
    config = {"responsive": True, "toImageButtonOptions": {"format": "svg", }}
    df2_piv=df.pivot(index=["ID", "Donor", "Experiment"], columns=["OTU"], values="RA").reset_index()
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
        split='Donor'
    else:
        split=None

    fig1=px.bar(df2_piv,x="ID",y=[x,y],facet_col=split)
    fig1.update_xaxes(matches=None)
    fig1.update_layout(showlegend=False)
    st.plotly_chart(fig1,use_container_width=True)
    st.markdown("""---""")
    st.subheader("Correlations")
    color=st.selectbox(label="Color By",options=[None,'ID','Donor','Experiment'],index=2)

    col1,col2=st.columns(2)
    corr_expander=st.expander(label='Spearman Correlation',expanded=True)
    fig2_trend=px.scatter(df2_piv,x=x,y=y,color=color,hover_data=["ID","Experiment"],trendline="ols",trendline_scope="overall")
    fig2_trend.layout.xaxis.title=fig2_trend.layout.xaxis.title['text'].split(';')[-1]
    fig2_trend.layout.yaxis.title=fig2_trend.layout.yaxis.title['text'].split(';')[-1]
    col1.plotly_chart(fig2_trend,use_container_width=True)    
    print_corr(df2_piv,x,y,'All Donor',corr_expander)
    fig2_trend_each=px.scatter(df2_piv,x=x,y=y,color=color,hover_data=["ID","Experiment"],trendline="ols",)
    fig2_trend_each.layout.xaxis.title=fig2_trend_each.layout.xaxis.title['text'].split(';')[-1]
    fig2_trend_each.layout.yaxis.title=fig2_trend_each.layout.yaxis.title['text'].split(';')[-1]

    col2.plotly_chart(fig2_trend_each,use_container_width=True)        
    for donor in df2_piv['Donor'].unique():
        print_corr(df2_piv[df2_piv['Donor'].isin([donor])],x,y,donor,corr_expander)
        
    st.subheader(f"{x.split(';')[-1]}:{y.split(';')[-1]} ratio")
    widget1,widget2=st.columns(2)
    split_by=widget1.selectbox(label="Group By",options=[None,'ID','Donor','Experiment'],index=2)
    color_by2=widget2.selectbox(label="Color By",options=[None,'ID','Donor','Experiment'],index=0)
    fig3=px.box(df2_piv,x=split_by,y="ratio",hover_data=["ID","Experiment"],color=color_by2)
    fig3.update_xaxes(matches=None,autorange=True)
    fig3.update_layout(boxmode='group',boxgap=0)
    fig3.layout.yaxis.title=f"{x.split(';')[-1]}:{y.split(';')[-1]} ratio"

    st.plotly_chart(fig3,use_container_width=True)
