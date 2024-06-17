#%%
from workflow.fetchdata import vmr_parse as vp 
from glob import glob
import re
import pandas as pd

import os
import shutil
import pickle as pkl
from tqdm import  tqdm

#%%
patten=r'[^a-zA-Z0-9]'
def curate_gbid(s:str):
    s=s.replace(' ','')
    s=re.sub(r'\s*\(.*?\)\s*','',s)
    s=s.split(':')[-1]
    return s
gbids={curate_gbid(i.split('|')[-1].replace('.pkl','')):i for i in glob('data/CoreData/genbank_meta/*') }

vmr=pd.read_csv('data/VMR_MSL38_v2.csv')
r_vmr=vmr[vmr['Genome composition'].apply(lambda x:'RNA' in x)].fillna('Null')
r_vmr['true_access']=r_vmr['Virus GENBANK accession'].apply(vp.get_genbank_id)
r_vmr['true_name']=r_vmr['Virus name abbreviation(s)'].apply(vp.get_correct_name)

#%%
def fetch_seq(f:str)->str:
    gb:dict=pkl.load(open(f,'rb'))
    seq:str=gb['GBSet']['GBSeq']['GBSeq_sequence']
    return seq

for name,access in tqdm(zip(r_vmr['true_name'],r_vmr['true_access'])):
    for k,v in access.items():
        if v=='Null':
            continue
        k='' if k=='_' else k
        segname=f'{name}|{k}|{v}'
        prev_fasta=f'data/CoreData/genome_fasta/{segname}:genome.fasta'
        new_fasta=f'data/CoreData/curated_genome_fasta/{segname}:genome.fasta'
        if os.path.exists(prev_fasta):
            # pass
            if not os.path.exists(new_fasta):
                shutil.copy(prev_fasta,new_fasta)
            # print(segname)
        else:
            if not os.path.exists(new_fasta):
                metafile=gbids[v]
                seq=fetch_seq(metafile)
                with open(new_fasta,'w') as f:
                    f.write(f'> {segname}\n{seq}\n')