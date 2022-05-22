# %%
import pandas as pd
import panel as pn
import os
import plotly.express as px
import itertools
pn.extension("plotly",
             sizing_mode="stretch_width",)
from scipy.stats import spearmanr
import streamlit as st
st.set_page_config(layout="wide")
hide_streamlit_style = """
            <style>
            
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True) 


st.title("Correlations")
main=st.container()
data=st.container()

# %%
# wd=f"C:\\Users\\Alon Shtrikman\\OneDrive - Mybiotix"
# wd=f"C:\\Users\\alons\\OneDrive - Mybiotix"
wd=""
# C:\Users\alons\OneDrive - Mybiotix
file="all_donors L6.csv"
# %%
df=pd.read_csv(os.path.join(wd,file)).rename(columns={"#OTU ID":"OTU",
                                                    #   "level_1":"ID",
                                                      "0":"RA"})[['OTU','RA','ID','Donor','Experiment']]
sorter=df.groupby("OTU").agg({'RA':'mean'}).sort_values("RA",ascending=False).index.tolist()
sorterIndex = dict(zip(sorter, range(len(sorter))))
df['rank']=df['OTU'].map(sorterIndex)
df=df.sort_values('rank')
# %%
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

# df.head()
# %%
# fig=px.bar(df,x="ID",y="0",color=")
# fig.update_layout(showlegend=False)
# %%

st.markdown("""---""")
# comb_list=list(itertools.combinations(top15, 2))
st.subheader("Select Bacteria to Compare")
x=st.selectbox(label="Bacteria 1",options=sorter[:top_val])
y=st.selectbox(label="Bacteria 2",options=sorter[:top_val],index=1)


# %%

    # x=top15[0]
    # y=top15[2]
    
    


# %%

config = {"responsive": True, "toImageButtonOptions": {"format": "svg", }}
# df2=df[df["].isin([x,y])]
df2_piv=df.pivot(index=["ID", "Donor", "Experiment"], columns=["OTU"], values="RA").reset_index()
df2_piv["ratio"]=df2_piv[x]/df2_piv[y]
title=f"{x.split(';')[-1]} to {y.split(';')[-1]} ratio per donor"
def print_corr(df,x,y,param):
    corr=spearmanr(df[[x,y]])
    # alert=pn.pane.Alert(f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]")
    if corr[1]<0.05:
        return st.sucess(f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]}")
    elif corr[1]>0.05:
        return st.error(f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]}")
    else:
        return st.info(f"{param} Correlation: Spearman R: {corr[0]}, P-Value: {corr[1]}")
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
st.plotly_chart(fig1,use_container_width=True)#config=config)
st.markdown("""---""")
st.subheader("Correlations")
color=st.selectbox(label="Color By",options=[None,'ID','Donor','Experiment'],index=2)
# fig2=px.scatter(df2_piv,x=x,y=y,color=color,title=title,hover_data=["ID","Experiment"],)#trendline="ols",trendline_scope="overall")
# st.plotly_chart(fig2,use_container_width=True)
col1,col2=st.columns(2)
fig2_trend=px.scatter(df2_piv,x=x,y=y,color=color,title=title,hover_data=["ID","Experiment"],trendline="ols",trendline_scope="overall")
fig2_trend.layout.xaxis.title=fig2_trend.layout.xaxis.title['text'].split(';')[-1]
fig2_trend.layout.yaxis.title=fig2_trend.layout.yaxis.title['text'].split(';')[-1]
col1.plotly_chart(fig2_trend,use_container_width=True)        
fig2_trend_each=px.scatter(df2_piv,x=x,y=y,color=color,title=title,hover_data=["ID","Experiment"],trendline="ols",)#marginal_x="box",marginal_y="box")#trendline_scope="overall")
fig2_trend_each.layout.xaxis.title=fig2_trend_each.layout.xaxis.title['text'].split(';')[-1]
fig2_trend_each.layout.yaxis.title=fig2_trend_each.layout.yaxis.title['text'].split(';')[-1]
col2.plotly_chart(fig2_trend_each,use_container_width=True)        
fig3=px.box(df2_piv,x="Donor",y="ratio",title=title,hover_data=["ID","Experiment"])#facet_col="Donor")
st.plotly_chart(fig3,use_container_width=True)

def plots(x,y):
        title=f"{x.split(';')[-1]} to {y.split(';')[-1]} ratio per donor"
        print(title)
        df2_piv["ratio"]=df2_piv[x]/df2_piv[y]
        fig1=px.bar(df2_piv,x="ID",y=[x,y])
        fig1.update_layout(showlegend=False)

        fig2=px.scatter(df2_piv,x=x,y=y,color="Donor",title=title,hover_data=["ID","Experiment"],)#trendline="ols",trendline_scope="overall")

        col_all=pn.Column()
        col=pn.Column()
        col_all.append(print_corr(df2_piv,x,y,"All Donors"))
        for donor in df2_piv.Donor.unique().tolist():
            df_temp=df2_piv[df2_piv["Donor"].isin([donor])]
            col.append(print_corr(df_temp,x,y,donor))
        fig2_trend=px.scatter(df2_piv,x=x,y=y,color="Donor",title=title,hover_data=["ID","Experiment"],trendline="ols",trendline_scope="overall")
        
        fig2_trend_each=px.scatter(df2_piv,x=x,y=y,color="Donor",title=title,hover_data=["ID","Experiment"],trendline="ols",)#marginal_x="box",marginal_y="box")#trendline_scope="overall")

        fig3=px.box(df2_piv,x="Donor",y="ratio",title=title,hover_data=["ID","Experiment"])#facet_col="Donor")

        fig3.update_xaxes(matches=None)

        fig4=px.strip(df2_piv,x="Donor",y="ratio",title=title,color="Experiment",hover_data=["ID","Experiment"])
        fig4.update_xaxes(matches=None)

        config = {"responsive": True, "toImageButtonOptions": {"format": "svg", }}
        return fig1#, fig2, fig3, fig4
        # return title, pn.Column(pn.panel(fig1,config=config),
        #                         pn.layout.Divider(),
        #                         pn.panel(fig2,config=config),
        #                         pn.layout.Divider(),
                                
        #                         pn.panel(fig2_trend,config=config),
        #                         pn.Row(col_all),
        #                         pn.layout.Divider(),
                                
        #                         pn.panel(fig2_trend_each,config=config),
        #                         pn.Row(col),
        #                         pn.layout.Divider(),
        #                         pn.panel(fig3,config=config),
        #                         pn.layout.Divider(),
        #                         pn.panel(fig4,config=config),sizing_mode="stretch_width")
        

# for (x,y) in comb_list:
#     title,plots_col= plots(x,y)
#     plots_col.save(title+".html")
# %%
