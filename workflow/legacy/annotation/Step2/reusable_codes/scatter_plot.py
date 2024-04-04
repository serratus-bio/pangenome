import pandas as pd
from pathlib import Path
from sklearn.linear_model import LinearRegression
import numpy as np

import matplotlib.pyplot as plt
from typing import List
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.axes import Axes
# import numpy as np

'''
(entry name for cov19:SARS-CoV-2||MN908947)
(accession for rdrp: PF00680)
'''

def parse_diamond(matchfile:str)->pd.DataFrame:
    '''
    matchfile:
        output of diamond blastp
        with '--header', 'verbose', 
        and 'tsv' suffix
    '''
    head_lines=open(matchfile,'r').readlines()[2].strip().split(': ')[1]
    diamond_aligns=pd.read_csv(matchfile,skiprows=3,names=head_lines.split(', '),delim_whitespace=True)
    return diamond_aligns

def get_identities():
    '''
    from diamond tsv
    to df of ['domain','virus','identity']
    '''
    o_dict={}
    _={}
    _['domain'],_['virus'],_['identity']=[],[],[]
    for i in Path('cov19-hits').iterdir():
        if i.suffix=='.tsv':
            stem=i.stem.replace('-match','')
            diamond=parse_diamond(i.absolute())
            o_dict[stem]={i:j for i,j in zip(diamond['Subject ID'],diamond['Percentage of identical matches'])}
            for idx,s in diamond.iterrows():
                _['domain'].append(stem)
                _['virus'].append(s['Subject ID'])
                _['identity'].append(s['Percentage of identical matches'])
    identities=pd.DataFrame(_)
    # identities.columns
    return o_dict,identities

def get_annotations():
    '''
    from domains.csv
    to unqiue accession-annotation .csv file
    
    then read in the csv
    '''
    domains=pd.read_csv('nido-domains.csv')
    o_dict={}
    o_dict['accession'],o_dict['annotation']=[],[]
    for accession in domains['domain_accession'].unique():
        annot=domains[domains['domain_accession']==accession].iloc[0]['domain_annotation']
        o_dict['accession'].append(accession)
        o_dict['annotation'].append(annot)
        
    o_df=pd.DataFrame(o_dict)
    o_df.to_csv('accession-annotation.csv',index=False)
    #--- ---#
    annotations=pd.read_csv('accession-annotation.csv',index_col='accession')
    return annotations
    

def cal_lr(x,y,refx=np.linspace(1,100,100)):
    lr=LinearRegression()
    x_=np.array(x).reshape(-1,1)
    y_=np.array(y).reshape(-1,1)
    reg=lr.fit(X=x_,y=y_)
    r2=reg.score(x_,y_)
    coef=reg.coef_[0][0]
    intercept=reg.intercept_[0]
    
    predy=reg.predict(refx.reshape(-1,1))
    sel=(predy>0) & (predy<100)
    
    reg_x=refx[sel.reshape(-1)].reshape(-1)
    reg_y=predy[sel].reshape(-1)
    
    return (
        coef,intercept,r2,
        reg_x,reg_y,
    )
    
def plot_scatter():
    '''
    cleansing needed
    '''
    o_dict=get_identities()[0]
    annotations=get_annotations()
    
    d:dict
    refx=refy=np.linspace(1,100,100)
    per_page=2
    with PdfPages('tmp1.pdf') as pdf:
        count=0
        fig,axs=plt.subplots(per_page,1,figsize=(10,10.5*per_page))
        axs:List[Axes]
        
        for accession,d in o_dict.items():
            if accession!='PF00680':
                x,y=[],[]
                for k in d.keys():
                    if k !='SARS-CoV-2||MN908947':
                        r=o_dict['PF00680'].get(k,None)
                        if r is not None:
                            x.append(r)
                            y.append(d[k])
                (coef,intercept,
                r2,reg_x,reg_y)=cal_lr(x+[100.],y+[100.])
                #for better regression
                axs[count].scatter(x,y)
                axs[count].plot(refx,refy,'--',color='grey')
                axs[count].plot(reg_x,reg_y,'--',color='red')
                axs[count].set_xlim(0,100)
                axs[count].set_ylim(0,100)
                annot=annotations['annotation'].loc[accession]
                axs[count].set_title(f'{accession}:{annot}')
                s=(f'annot: {annot}\n'
                f'coef:{coef:.1f}\n'
                f'R^2:{r2:.2f}\n'
                f'freq:{len(x)}\n'
                f'intercept:{intercept:.1f}'
                )
                axs[count].text(x=5,y=80,s=s,fontsize=12)
                axs[count].set_aspect(1,'box')
                count+=1
                if count==per_page:
                    count=0
                    plt.tight_layout()
                    pdf.savefig(fig)
                    plt.close(fig)
                    fig,axs=plt.subplots(per_page,1,figsize=(10,10.5*per_page))
                    axs:List[Axes]
        if count!=0:
            while count <per_page:
                axs[count].set_axis_off()
                count+=1
            pdf.savefig(fig)
        plt.close(fig)
                    # break