from dna_features_viewer import GraphicFeature, GraphicRecord
from foreground import get_foreground 
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib.patches as patches
from matplotlib import colormaps
from typing import List
from matplotlib.axes import Axes
import seaborn as sns
import pandas as pd

sequence_length=42000
ref=21000

domains=pd.read_csv('nido-domains.csv')
sample_size=len(domains['genome_name'].unique())

coloring = {
 'orf1a': ['bCoV_NSP1', 'CoV_NSP2_N', 'CoV_NSP2_C',
           'bCoV_NSP3_N', 'Macro', 'bCoV_SUD_M', 'bCoV_SUD_C', 
           'CoV_peptidase', 'bCoV_NAR',  'CoV_NSP3_C', 
           'CoV_NSP4_N', 'CoV_NSP4_C', 'Peptidase_C30', 
           'CoV_NSP6', 'CoV_NSP7', 
           'CoV_NSP8', 'CoV_NSP9', 'CoV_NSP10'],
 'rdrp_helicase': [ 'CoV_RPol_N', 'RdRP_1', 'Viral_helicase1'  ],
 'orf1b':         [ 'CoV_Methyltr_1',  'CoV_NSP15_N', 'CoV_NSP15_M',
                    'CoV_NSP15_C', 'CoV_Methyltr_2'],
 '3p_orfs': [       'bCoV_S1_N', 'bCoV_S1_RBD', 'CoV_S1_C',
                    'CoV_S2',  'bCoV_viroporin', 'CoV_E', 'CoV_M', 
                    'bCoV_NS6', 'bCoV_NS7A', 'bCoV_NS7B', 'bCoV_NS8', 
                    'bCoV_lipid_BD', 'CoV_nucleocap','bCoV_Orf14']
}
inv_coloring  = dict([(y,x) for x in coloring for y in coloring[x]])
dpalette = dict()
all_names = set([x for x in inv_coloring])
for region in coloring:
    #region_color = {'orf1a': 'BuPu', 'rdrp_helicase': 'Greens', 'orf1b': 'GnBu', '3p_orfs': 'Greys'}[region]
    #palette = sns.color_palette(region_color, len(coloring[region]) + offset_color)
    region_color, additional_length, offset_color = {'orf1a': ('cubehelix',20,10), 'rdrp_helicase': ('Greens',10,5), 'orf1b': ('BrBG',30,30), '3p_orfs': ('Greys',10,7)}[region]
    palette = sns.color_palette(region_color, len(coloring[region]) + additional_length)
    if region == 'orf1a': palette = palette[::-1]
    for i,name in enumerate(coloring[region]):
        dpalette[name] = palette[offset_color+i] 

# some manual palette fixes
dpalette['Macro'] = "#d42065"
dpalette['Peptidase_C30'] = "#c2da29"
dpalette['RdRP_1'] = "#23d32a"

nido_exclusives=[]
for i in domains['domain_annotation'].apply(lambda x:x.split(':')[0]).unique():
    if i not in dpalette:
        nido_exclusives.append(i)
nidopalette = sns.color_palette('YlGnBu', len(nido_exclusives))

def hide_axes(ax,genome_start_pos,genome_end_pos):
    #hack: since I don't~~~~~ didn't know how to erase axes around coordinates, draw a white box over them.
    # is that a ugly hack? yes. does it render well? also yes
    ax.add_patch(
     patches.Rectangle(
        (-100, -0.25),
        genome_start_pos,
        0.5,
        fill=True,      # remove background
        facecolor='white'
     ) ) 
    ax.add_patch(
     patches.Rectangle(
        (genome_end_pos, -0.25),
        sequence_length+100,
        0.5,
        fill=True,      # remove background
        facecolor='white'
     ) ) 
    
import pandas as pd
sort_rdrps=pd.read_csv('nido-rdrp-sort.list')

plt.close('all')

c=0
fig,ax=plt.subplots(len(sort_rdrps),1,sharex=True,figsize=(60,len(sort_rdrps)*2))
ax:List[Axes]
for i,j in zip(sort_rdrps['id'],sort_rdrps['match']):
    idf=domains[domains['genome_name']==i]
    # for i,idf in domains.groupby('genome_name'):
        # sequence_length=idf.iloc[0]['genome_length']
    if (idf['domain_accession']=='PF00680').any():
        rdrp_start=idf[idf['domain_accession']=='PF00680'].iloc[0]['start']
        rdrp_hmm_start=idf[idf['domain_accession']=='PF00680'].iloc[0]['hmmStart']
        del_start=ref-(rdrp_start-rdrp_hmm_start*3)
    else:
        del_start=ref-idf.iloc[0]['genome_length']/2
        print(f'no rdrp warning:{i}')
    
    features=[]
    for _,d in idf.iterrows():
        strand=+1 if d['strand']=='SENSE' else -1
        a_=d['domain_annotation'].split(':')[0]
        if  a_ in dpalette:
            color=dpalette[a_]
        elif a_ in nido_exclusives:
            color=nidopalette[nido_exclusives.index(a_)]
        else:
            color='red'
        # color= 'green' if d['domain_accession']=='PF00680' else 'grey'
        features.append(GraphicFeature(start=d['start']+del_start, end=d['end']+del_start, strand=strand, color=color,
                    label=d['domain_accession']))
    record = GraphicRecord(sequence_length=sequence_length, features=features)
    record.plot(ax=ax[c],figure_width=20)
    ax[c].set_facecolor("white")
    ax[c].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # remove ticks
    ax[c].annotate(i+f'|{j}' ,(-1000,-0), xycoords='data', annotation_clip=False, ha='left',fontsize=25) 
    hide_axes(ax[c],0+del_start,idf.iloc[0]['genome_length']+del_start)
    c+=1
    if c>=sample_size:
        break
fig.tight_layout()
# fig.set_dpi(400)
fig.savefig('xxxx-sort.pdf')
plt.close('all')